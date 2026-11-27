# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *              Yunior C. Fonseca Reyna (cfonseca@cnb.csic.es)
# *
# * MRC Laboratory of Molecular Biology (MRC-LMB)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
import os
from typing import Tuple, Set, List

import numpy as np
import math
import scipy.fft
import xml.etree.ElementTree as eTree

from pwem import logger
from pwem.convert import transformations
from pwem.convert.transformations import translation_from_matrix, euler_from_matrix
from tomo.constants import TR_RELION
from tomo.objects import CTFTomoSeries, TiltSeries
from warp import TILTIMAGES_FOLDER, FRAMESERIES_FOLDER

""" This code is adapted from https://github.com/dtegunov/tom_deconv
    and https://github.com/Heng-Z/IsoNet

    Original code is from:

    Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
    Journal of Structural Biology, 149 (2005), 227-234.

    Copyright (c) 2004-2007
    TOM toolbox for Electron Tomography
    Max-Planck-Institute of Biochemistry
    Dept. Molecular Structural Biology
    82152 Martinsried, Germany
    http://www.biochem.mpg.de/tom
"""


def tom_ctf1d(length, pixelsize, voltage, cs, defocus,
              amplitude, phaseshift, bfactor):
    ny = 1 / pixelsize
    lambda1 = 12.2643247 / np.sqrt(voltage * (1.0 + voltage * 0.978466e-6)) * 1e-10
    lambda2 = lambda1 * 2

    k2 = (np.arange(length, dtype=np.float32) / (2 * length) * ny) ** 2
    term1 = lambda1 ** 3 * cs * k2 ** 2

    w = np.pi / 2 * (term1 + lambda2 * defocus * k2) - phaseshift

    acurve = np.cos(w) * amplitude
    pcurve = -np.sqrt(1 - amplitude ** 2) * np.sin(w)
    bfactor = np.exp(-bfactor * k2 * 0.25)
    ctf = (pcurve + acurve) * bfactor

    return ctf


def tom_deconv(vol, angpix=1.0, voltage=300.0, cs=2.7, defocus=3, snrfalloff=1.1,
               deconvstrength=1, highpassnyquist=0.02, phaseflipped=False,
               phaseshift=0, ncpu=1, gpu=False, gpuid=0):
    """
    :param vol: tomogram volume (or 2D image)
    :param angpix: angstrom per pixel
    :param voltage: voltage in kilovolts
    :param cs: cs in millimeters
    :param defocus: defocus in micrometers, positive = underfocus
    :param snrfalloff: how fast does SNR fall off, i.e. higher values will downweight high frequencies; values like 1.0 or 1.2 seem reasonable
    :param deconvstrength: how much will the signal be deconvoluted overall, i.e. a global scale for SNR; exponential scale: 1.0 is SNR = 1000 at zero frequency, 0.67 is SNR = 100, and so on
    :param highpassnyquist: fraction of Nyquist frequency to be cut off on the lower end (since it will be boosted the most)
    :param phaseflipped: whether the data are already phase-flipped
    :param phaseshift: CTF phase shift in degrees (e.g. from a phase plate)
    :param ncpu: number of CPUs for FFT
    :param gpu: use GPU instead
    :param gpuid: gpu ID

    Example:
    deconv = tom_deconv(mytomo, 3.42, 300, 2.7, 6, 1.1, 1, 0.02, False, 0);

    """

    # Precompute highpass filter
    data = np.arange(0, 1 + 1 / 2047, 1 / 2047, dtype=np.float32)
    highpass = np.minimum(1, data / highpassnyquist) * np.pi
    highpass = 1 - np.cos(highpass)

    snr = np.exp(-data * snrfalloff * 100 / angpix) * (10 ** (3 * deconvstrength)) * highpass + 1e-6

    # Precompute some constants
    angpix *= 1e-10
    voltage *= 1000
    cs *= 1e-3
    defocus = -defocus * 1e-6
    phaseshift = phaseshift / 180 * np.pi

    ctf = tom_ctf1d(2048, angpix, voltage, cs, defocus, 0.07, phaseshift, 0)
    if phaseflipped:
        ctf = np.abs(ctf)
    wiener = ctf / (ctf * ctf + 1.0 / snr)

    s1, f1 = -math.floor(vol.shape[0] / 2), -math.floor(vol.shape[0] / 2) + vol.shape[0] - 1
    m1 = np.arange(s1, f1 + 1, dtype=np.int32)
    s2, f2 = -math.floor(vol.shape[1] / 2), -math.floor(vol.shape[1] / 2) + vol.shape[1] - 1
    m2 = np.arange(s2, f2 + 1, dtype=np.int32)

    if vol.ndim == 3:
        s3, f3 = -math.floor(vol.shape[2] / 2), -math.floor(vol.shape[2] / 2) + vol.shape[2] - 1
        m3 = np.arange(s3, f3 + 1, dtype=np.int32)
        x, y, z = np.meshgrid(m1, m2, m3, indexing='ij')
        x = np.divide(x, abs(s1))
        y = np.divide(y, abs(s2))
        z = np.divide(z, max(1, abs(s3)))
        r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    else:
        x, y = np.meshgrid(m1, m2, indexing='ij')
        x = np.divide(x, abs(s1))
        y = np.divide(y, abs(s2))
        r = np.sqrt(x ** 2 + y ** 2)

    r = np.minimum(1, r)
    r = np.fft.ifftshift(r)
    ramp = np.interp(r, data, wiener).astype(np.float32)
    vol = vol.astype(np.float32)

    if gpu:
        import cupy as cp
        with cp.cuda.Device(gpuid):
            deconv = cp.real(cp.fft.ifftn(cp.fft.fftn(cp.asarray(vol)) * cp.asarray(ramp)))
            deconv = cp.asnumpy(deconv).astype(np.float32)
    else:
        deconv = np.real(scipy.fft.ifftn(scipy.fft.fftn(vol, overwrite_x=True, workers=ncpu) * ramp,
                                         overwrite_x=True, workers=ncpu))

    return deconv


def _getDefocusHand(xmlPath) -> int:
    """Read the CTF handedness value from the .xml file generated by
    WarpTools ts_ctf. Code:
    - AreAnglesInverted = 'False' means _rlnTomoHand = -1.
    - AreAnglesInverted = 'True' means _rlnTomoHand = 1.
    """
    tree = eTree.parse(xmlPath)
    root = tree.getroot()
    areAnglesInverted = root.get('AreAnglesInverted')
    return 1 if areAnglesInverted == 'True' else -1


def _parseCtfGrid(root, gridName, scale=1.0):
    grid = root.find(gridName)
    if grid is None:
        return {}

    return {
        int(node.get("Z")): float(node.get("Value")) * scale
        for node in grid.findall("Node")
    }


def parseCtfXMLFile(xmlPath):
    """Parse Warp tilt-series CTF data."""
    tree = eTree.parse(xmlPath)
    root = tree.getroot()

    ctfParams = root.find('CTF')
    ctfData = {
        param.get('Name'): param.get('Value')
        for param in ctfParams.findall('Param')
    }

    gridCtf = root.find('GridCTF')
    gridCtfData = {
        "Width": gridCtf.get("Width") if gridCtf is not None else None,
        "Height": gridCtf.get("Height") if gridCtf is not None else None,
        "Depth": gridCtf.get("Depth") if gridCtf is not None else None,
        "Nodes": _parseCtfGrid(root, 'GridCTF', 1e4),
        "DeltaNodes": _parseCtfGrid(root, 'GridCTFDefocusDelta', 1e4),
        "AngleNodes": _parseCtfGrid(root, 'GridCTFDefocusAngle')
    }

    return ctfData, gridCtfData


def extractGlobalResolution(xmlPath):
    globalResolution = ''

    with open(xmlPath, "r", encoding="utf-8") as f:
        xml_data = f.read().lstrip()

    root = eTree.fromstring(xml_data)

    for param in root.findall('Param'):
        if param.attrib.get('Name') == 'GlobalResolution':
            globalResolution = f"{float(param.attrib.get('Value')):.2f}"
            break
    return globalResolution


def extractAxisOffsets(xmlPath):
    """
    Extracts AxisOffsetX and AxisOffsetY values from an XML file.
    Args:
        xmlPath (str): Path to the XML file.
    Returns:
        tuple: (axisOffsetX, axisOffsetY), both as lists of floats.
    """
    tree = eTree.parse(xmlPath)
    root = tree.getroot()

    # Get text content of the XML elements
    axisOffsetXText = root.findtext('AxisOffsetX')
    axisOffsetYText = root.findtext('AxisOffsetY')

    # Convert text to lists of floats
    axisOffsetX = [float(value) for value in axisOffsetXText.strip().split()]
    axisOffsetY = [float(value) for value in axisOffsetYText.strip().split()]

    return axisOffsetX, axisOffsetY


def updateCtFXMLFile(xmlPath, ctfTomoSeries, ts):
    tree = eTree.parse(xmlPath)
    root = tree.getroot()

    ctfByAcqOrder = {
        ctf.getAcquisitionOrder(): ctf
        for ctf in ctfTomoSeries.iterItems(iterate=False)
        if ctf.isEnabled()
    }

    tiltImages = sorted(
        (ti for ti in ts.iterItems(iterate=False) if ti.isEnabled()),
        key=lambda ti: ti.getTiltAngle()
    )

    defocusValues = []
    defocusDeltaValues = []
    defocusAngleValues = []

    for ti in tiltImages:
        acqOrder = ti.getAcquisitionOrder()
        ctf = ctfByAcqOrder.get(acqOrder)

        if ctf is None:
            raise ValueError(
                f'No CTF found for acquisition order {acqOrder} '
                f'in tilt-series {ts.getTsId()}'
            )

        defocusU = ctf.getDefocusU() / 1e4
        defocusV = ctf.getDefocusV() / 1e4

        defocusValues.append((defocusU + defocusV) / 2)
        defocusDeltaValues.append((defocusU - defocusV) / 2)
        defocusAngleValues.append(ctf.getDefocusAngle())

    def _updateGrid(gridName, values):
        grid = root.find(gridName)
        if grid is None:
            return

        nodes = sorted(
            grid.findall('Node'),
            key=lambda node: int(node.get('Z'))
        )

        if len(nodes) != len(values):
            raise ValueError(
                f'{gridName} contains {len(nodes)} nodes but '
                f'{len(values)} CTF values were provided'
            )

        for node, value in zip(nodes, values):
            node.set('Value', str(value))

    # Global values are fallbacks; Warp uses the per-tilt grids below.
    for param in root.find('CTF').findall('Param'):
        name = param.get('Name')

        if name == 'Defocus':
            param.set('Value', str(np.mean(defocusValues)))
        elif name == 'DefocusDelta':
            param.set('Value', str(np.mean(defocusDeltaValues)))

    _updateGrid('GridCTF', defocusValues)
    _updateGrid('GridCTFDefocusDelta', defocusDeltaValues)
    _updateGrid('GridCTFDefocusAngle', defocusAngleValues)

    tree.write(xmlPath, encoding='utf-8', xml_declaration=True)


def getWarpTiltImages(ts: TiltSeries, presentAcqOrders: List[int]):
    """Tilt images of ``ts`` in the order Warp indexes them, sorted by ascending tilt
    angle. Only tilts that are enabled and whose acquisition order is in
    ``presentAcqOrders`` are kept. That list is meant to come from
    tomo.utils.getPresentAcqOrdersTsAndCtf, so each kept tilt is also enabled in the
    CTF series and has a CTF to write. Every per-series file Warp matches by tilt index
    (.tomostar rows, CTF XML grid nodes, IMOD .xf/.tlt lines) must be written from
    this same list."""
    return sorted(
        (ti for ti in ts.iterItems(iterate=False) if ti.isEnabled() and ti.getAcquisitionOrder() in presentAcqOrders),
        key=lambda ti: ti.getTiltAngle()
    )


def writeCtfXml(xmlPath: str,
                ctfTomoSeries: CTFTomoSeries,
                ts: TiltSeries,
                presentAcqOrders: List[int],
                apix: float,
                voltage: int,
                cs: float,
                amplitude: float,
                imageDimPix: Tuple[int, int],
                areAnglesInverted: bool= False) -> None:
    """Build a Warp per-series CTF XML directly from Scipion's (already estimated)
    CTF metadata, replacing the WarpTools ts_ctf estimation step.

    Only the elements that ts_export_particles / ts_import_alignments / M actually
    read are written: the AreAnglesInverted flag, the <CTF> parameter block and the
    four per-tilt CTF grids (defocus mean, delta, angle, phase). Everything else Warp
    default-initialises, reads from the .tomostar, or fills in during
    ts_import_alignments (which re-serialises this file to the full schema).

    Conventions (must stay in lockstep with the rest of the plugin):
      * Defocus is stored in MICRONS (Scipion keeps Angstrom -> /1e4), split into the
        symmetric mean (GridCTF) and the astigmatism half-difference (GridCTFDefocusDelta),
        matching updateCtFXMLFile.
      * Grid Z-nodes follow getWarpTiltImages (enabled tilts, ascending tilt angle), the
        same list used to write the .tomostar, so grid index z lines up with Warp's
        internal tilt index.
      * CTFs are matched to tilts by acquisition order, like updateCtFXMLFile.
    """
    ctfByAcqOrder = {
        ctfAcqOrder: ctf
        for ctf in ctfTomoSeries.iterItems(iterate=False)
        if ctf.isEnabled() and (ctfAcqOrder := ctf.getAcquisitionOrder()) in presentAcqOrders
    }

    tiltImages = getWarpTiltImages(ts, presentAcqOrders)

    defocusValues = []
    defocusDeltaValues = []
    defocusAngleValues = []
    for ti in tiltImages:
        # Intersect the ts and the ctf via the acquisition order, as there may
        # be excluded views
        acqOrder = ti.getAcquisitionOrder()
        ctf = ctfByAcqOrder[acqOrder]
        defocusU = ctf.getDefocusU() / 1e4
        defocusV = ctf.getDefocusV() / 1e4
        defocusValues.append((defocusU + defocusV) / 2)
        defocusDeltaValues.append((defocusU - defocusV) / 2)
        defocusAngleValues.append(ctf.getDefocusAngle())

    nTilts = len(tiltImages)
    phaseValues = [0.0] * nTilts
    nx, ny = imageDimPix

    root = eTree.Element('TiltSeries', {
        'AreAnglesInverted': 'True' if areAnglesInverted else 'False',
        'ImageDimensionsAngstrom': f'{nx * apix}, {ny * apix}',
    })
    ctfEl = eTree.SubElement(root, 'CTF')

    def _param(name, value):
        eTree.SubElement(ctfEl, 'Param', {'Name': name, 'Value': str(value)})

    _param('PixelSize', apix)
    _param('Voltage', voltage)
    _param('Cs', cs)
    _param('Amplitude', amplitude)
    # Global Defocus/DefocusDelta are fallbacks; Warp drives CTF from the grids below.
    _param('Defocus', np.mean(defocusValues))
    _param('DefocusDelta', np.mean(defocusDeltaValues))
    _param('DefocusAngle', defocusAngleValues[nTilts // 2])
    _param('PhaseShift', 0)

    def _grid(name, values):
        grid = eTree.SubElement(root, name, {
            'Width': '1', 'Height': '1', 'Depth': str(nTilts),
            'MarginX': '0', 'MarginY': '0', 'MarginZ': '0',
        })
        for z, value in enumerate(values):
            eTree.SubElement(grid, 'Node',
                             {'X': '0', 'Y': '0', 'Z': str(z), 'Value': str(value)})

    # Element name is GridCTF (its property name in Warp is GridCTFDefocus).
    _grid('GridCTF', defocusValues)
    _grid('GridCTFDefocusDelta', defocusDeltaValues)
    _grid('GridCTFDefocusAngle', defocusAngleValues)
    _grid('GridCTFPhase', phaseValues)

    eTree.ElementTree(root).write(xmlPath, encoding='utf-8', xml_declaration=True)


def tomoStarGenerate(tsId, tiValues, otputFolder, isTiltSeries, perTs=False):
    """Generate the .tomostar files from TS"""
    _fileName = os.path.abspath(otputFolder) + '/%s.tomostar' % tsId
    _file = open(_fileName, 'a+')
    header = """
data_

loop_
_wrpMovieName #1
_wrpAngleTilt #2
_wrpAxisAngle #3
_wrpDose #4
_wrpAverageIntensity #5  
_wrpMaskedFraction #6
"""
    _file.write(header)
    sortedTiValues = sorted(tiValues)
    imagesFolder = TILTIMAGES_FOLDER if isTiltSeries else FRAMESERIES_FOLDER
    for acqOrder in sortedTiValues:
        value = tiValues[acqOrder]
        tiPath = '../%s/' % imagesFolder + value[0] if not perTs else '../../%s/' % imagesFolder + value[0]
        angleTilt = value[1]
        axisAngle = value[2]
        shiftX = f"{value[3]:.6f}"
        shiftY = f"{value[4]:.6f}"
        dose = value[5]
        averageIntensity = 3.721
        maskedFraction = value[7]
        _file.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (
            tiPath, angleTilt, axisAngle, dose, averageIntensity, maskedFraction))

    _file.close()


def genTransformMatrix(rot, tilt, psi):
    angles = (float(rot), float(tilt), float(psi))
    radAngles = -np.deg2rad(angles)
    M = transformations.euler_matrix(radAngles[0], radAngles[1], radAngles[2], 'szyz')

    # These 3 lines are the ones for "invert" flag.
    M[0, 3] = 0
    M[1, 3] = 0
    M[2, 3] = 0
    M = np.linalg.inv(M)

    return M


def getTransformInfoFromCoordOrSubtomo(coord, samplingRate, mirrorZ=False):
    """RELION Euler angles (rot, tilt, psi; degrees) and shifts (Å) of a coordinate or
    subtomogram.

    :param mirrorZ: express the orientation in the Z-mirrored tomogram frame (z -> -z),
        to go with Z-mirrored positions. A mirror can't be applied to a rotation alone,
        so the pose is conjugated by the mirror F = diag(1, 1, -1): R' = F R F. That is
        the same particle pose seen in the mirrored tomogram, provided the reference is
        mirrored too (the map refined from these particles is the enantiomer of the
        one in the original frame). In Euler terms the tilt changes sign, which is
        (rot + 180, tilt, psi + 180) in the usual range; the Z shift changes sign too.
    """
    M = coord.getMatrix(convention=TR_RELION)
    if mirrorZ:
        F = np.diag([1., 1., -1., 1.])
        M = F @ M @ F
    shifts = translation_from_matrix(M)
    # These 2 lines below were done when inverting, which is now what we always do.
    shifts = -shifts
    M = np.linalg.inv(M)
    angles = -np.rad2deg(euler_from_matrix(M, axes='szyz'))
    shifts *= samplingRate

    return angles, shifts


def modifyOptFileMultiTable(filePath, columnToModify, modifierFunc):
    """
        Modify a given column in all data blocks in a .star file that contain a loop_ table.
        Keeps untouched blocks without loop_.
        """
    tempPath = filePath + '.tmp'

    with open(filePath, 'r') as fin, open(tempPath, 'w') as fout:

        for line in fin:
            stripped = line.strip()
            if stripped.startswith('_'):
                if stripped.startswith(columnToModify):
                    parts = line.strip().split()
                    try:
                        oldValue = parts[-1]
                        parts[-1] = str(modifierFunc(oldValue))
                        fout.write(' '.join(parts) + '\n')
                    except IndexError:
                        fout.write(line)  # Malformed line; write as-is
                    continue
            fout.write(line)

    os.replace(tempPath, filePath)


def modifyStarFileMultiTable(filePath, columnToModify, modifierFunc):
    """
    Modify a given column in all data blocks in a .star file that contain a loop_ table.
    Keeps untouched blocks without loop_.
    """
    tempPath = filePath + '.tmp'

    with open(filePath, 'r') as fin, open(tempPath, 'w') as fout:
        inLoop = False
        headers = []
        colIndex = -1

        for line in fin:
            stripped = line.strip()

            # Handle start of new data block
            if stripped.startswith('data_'):
                inLoop = False
                headers = []
                colIndex = -1
                fout.write(line)
                continue

            # Start of loop
            if stripped.startswith('loop_'):
                inLoop = True
                headers = []
                colIndex = -1
                fout.write(line)
                continue

            # Header lines
            if inLoop and stripped.startswith('_'):
                headers.append(stripped)
                fout.write(line)
                if stripped.startswith(columnToModify):
                    colIndex = len(headers) - 1  # 0-based
                continue

            # Data line in a loop
            if inLoop and stripped and not stripped.startswith('_'):
                if colIndex == -1:
                    fout.write(line)  # Column not in this block
                    continue

                parts = line.strip().split()
                try:
                    oldValue = parts[colIndex]
                    parts[colIndex] = str(modifierFunc(oldValue))
                    fout.write(' '.join(parts) + '\n')
                except IndexError:
                    fout.write(line)  # Malformed line; write as-is
                continue

            # Any other line (outside of loop)
            fout.write(line)

    os.replace(tempPath, filePath)


# --------------------------- Warp tilt-series geometry ---------------------------------------
# Python port of the image-position part of Warp's TiltSeries.GetPositionInAllTilts
# (WarpLib/TiltSeries/TiltSeries.cs), for the geometry stored in a <tsId>.xml after
# ts_import_alignments. Warp's deformation models (GridVolumeWarpX/Y/Z, GridMovementX/Y)
# are not ported: they stay at zero until Warp refines the tilt series (e.g. in M), and
# readWarpTsGeometry reports whether they are set.
WARP_DEFORMATION_GRIDS = ('GridVolumeWarpX', 'GridVolumeWarpY', 'GridVolumeWarpZ',
                          'GridMovementX', 'GridMovementY')


def readWarpTsGeometry(xmlPath: str) -> dict:
    """Read the per-tilt geometry of a Warp tilt-series XML, in Warp's tilt order
    (the .tomostar rows). Angles in degrees, offsets and dimensions in Å.
    'hasDeformation' is True if any deformation grid has a non-zero node."""
    root = eTree.parse(xmlPath).getroot()

    def _floats(tag):
        el = root.find(tag)
        if el is None or not (el.text or '').strip():
            return np.array([])
        return np.array([float(v) for v in el.text.split()])

    hasDeformation = any(
        float(node.get('Value', 0)) != 0
        for gridName in WARP_DEFORMATION_GRIDS
        for grid in root.findall(gridName)
        for node in grid.findall('Node')
    )
    imageDims = np.array([float(v) for v in root.get('ImageDimensionsAngstrom').split(',')])
    return {
        'angles': _floats('Angles'),
        'axisAngles': _floats('AxisAngle'),
        'offsetX': _floats('AxisOffsetX'),
        'offsetY': _floats('AxisOffsetY'),
        'levelAngleX': float(root.get('LevelAngleX', 0)),
        'levelAngleY': float(root.get('LevelAngleY', 0)),
        'imageDimsAngst': imageDims,
        'hasDeformation': hasDeformation,
    }


def _warpEuler(rot: float, tilt: float, psi: float) -> np.ndarray:
    """Warp's Matrix3.Euler (angles in radians), WarpLib/Tools/Matrix3.cs."""
    ca, cb, cg = math.cos(rot), math.cos(tilt), math.cos(psi)
    sa, sb, sg = math.sin(rot), math.sin(tilt), math.sin(psi)
    cc, cs, sc, ss = cb * ca, cb * sa, sb * ca, sb * sa
    return np.array([[cg * cc - sg * sa, cg * cs + sg * ca, -cg * sb],
                     [-sg * cc - cg * sa, -sg * cs + cg * ca, sg * sb],
                     [sc, ss, cb]])


def _warpRotateX(angle: float) -> np.ndarray:
    """Warp's Matrix3.RotateX (radians). Its constructor takes the elements column by
    column: (1, 0, 0, 0, c, s, 0, -s, c) is [[1, 0, 0], [0, c, -s], [0, s, c]]."""
    c, s = math.cos(angle), math.sin(angle)
    return np.array([[1., 0., 0.], [0., c, -s], [0., s, c]])


def projectWarpPositions(geom: dict, coordsAngst: np.ndarray,
                         volumeDimsAngst: np.ndarray) -> np.ndarray:
    """Image positions (Å, from the corner of each tilt image) of 3D positions given
    in Å from the corner of Warp's volume, as Warp extracts them. Returns an array of
    shape (nParticles, nTilts, 2). Same model as TiltSeries.GetPositionInAllTilts:
    centre on the volume, rotate by Euler(0, angle + LevelAngleY, -axisAngle) *
    RotateX(LevelAngleX), add the per-tilt axis offset and the image centre. The
    defocus handedness (AreAnglesInverted) only changes the depth (defocus), not X/Y."""
    coords = np.atleast_2d(np.asarray(coordsAngst, dtype=float))
    centered = coords - np.asarray(volumeDimsAngst, dtype=float) / 2
    imageCenter = geom['imageDimsAngst'] / 2
    rotX = _warpRotateX(math.radians(geom['levelAngleX']))
    nTilts = len(geom['angles'])
    result = np.empty((len(coords), nTilts, 2))
    for t in range(nTilts):
        rotation = _warpEuler(0.,
                              math.radians(geom['angles'][t] + geom['levelAngleY']),
                              math.radians(-geom['axisAngles'][t])) @ rotX
        transformed = centered @ rotation.T
        result[:, t, 0] = transformed[:, 0] + geom['offsetX'][t] + imageCenter[0]
        result[:, t, 1] = transformed[:, 1] + geom['offsetY'][t] + imageCenter[1]
    return result


def writeWarpTiltImages(ts: TiltSeries, tiltImages: List, imagesDir: str, averagesDir: str) -> List[str]:
    """Write each tilt of ``tiltImages`` (getWarpTiltImages) as a single-image .mrc in
    ``imagesDir`` and link it from ``averagesDir``, and return the file names, in the
    same order.

    Warp opens each tilt as a "movie" at <tomostar dir>/<wrpMovieName>: it reads the
    header of that file (size, data hash), then the image itself from the sibling
    average/ folder (TiltSeries.LoadMovieData / LoadMovieSizes), which must hold a
    single layer. The names are <tsId>_<stack index>.mrc, so the same tilt list
    always gives the same names."""
    # Imported here: pwem's image readers aren't needed by the rest of the module
    from pwem.emlib.image.image_readers import ImageStack, ImageReadersRegistry
    from pyworkflow.utils import createLink
    tsId = ts.getTsId()
    # The pixel size goes to the MRC header: Warp takes the physical image size from it.
    properties = {'sr': ts.getSamplingRate()}
    imageNames = []
    for ti in tiltImages:
        imageName = getWarpTiltImageName(tsId, ti)
        tiltImage = ImageStack(properties=properties)
        tiltImage.append(ImageReadersRegistry.open(f'{ti.getIndex()}@{ti.getFileName()}'))
        ImageReadersRegistry.write(tiltImage, os.path.join(imagesDir, imageName))
        # createLink resolves the source against the current directory, not against
        # the link's folder, so it must get the image's real path; it then writes the
        # link as relative (../<imageName>).
        createLink(os.path.join(imagesDir, imageName), os.path.join(averagesDir, imageName))
        imageNames.append(imageName)
    return imageNames


def getWarpTiltImageName(tsId: str, ti) -> str:
    """File name of a tilt image written by writeWarpTiltImages."""
    return f'{tsId}_{ti.getIndex():03d}.mrc'


def readWarpTsMoviePaths(xmlPath: str) -> List[str]:
    """Per-tilt image paths (MoviePath) of a Warp tilt-series XML, in Warp's tilt order."""
    node = eTree.parse(xmlPath).getroot().find('MoviePath')
    if node is None or not node.text:
        return []
    return [line.strip() for line in node.text.split('\n') if line.strip()]


def copyWarpTsXml(srcXml: str, dstXml: str, dataDirectory: str):
    """Copy a Warp tilt-series XML, pointing it to another data (.tomostar) folder.

    DataDirectory is an absolute path, and the per-tilt image paths (MoviePath) are
    relative to it, so this is the only change needed to move a series to another
    project with the same layout (<dataDirectory>/../tiltimages)."""
    tree = eTree.parse(srcXml)
    tree.getroot().set('DataDirectory', os.path.abspath(dataDirectory))
    tree.write(dstXml, encoding='utf-8', xml_declaration=True)


def _getWarpSettingsImportParam(root, name: str):
    node = root.find(f"./Import/Param[@Name='{name}']")
    if node is None:
        raise ValueError(f'Warp settings have no Import parameter {name}')
    return node


def copyWarpSettings(srcSettings: str, dstSettings: str, dataFolder: str, processingFolder: str):
    """Copy a Warp .settings file, pointing it to other data (.tomostar) and processing
    (per-series XML) folders. Every other value is kept, in particular the pixel size
    and the volume size (Tomo/Dimensions*), which must stay those the particles were
    placed with."""
    tree = eTree.parse(srcSettings)
    root = tree.getroot()
    _getWarpSettingsImportParam(root, 'DataFolder').set('Value', os.path.abspath(dataFolder))
    _getWarpSettingsImportParam(root, 'ProcessingFolder').set('Value', os.path.abspath(processingFolder))
    tree.write(dstSettings, encoding='utf-8', xml_declaration=True)


def readWarpSettingsTomoDims(settingsFile: str) -> Tuple[int, int, int]:
    """Volume size (Tomo/DimensionsX, Y, Z) of a Warp .settings file, in pixels of the
    settings' pixel size."""
    root = eTree.parse(settingsFile).getroot()
    dims = []
    for axis in 'XYZ':
        node = root.find(f"./Tomo/Param[@Name='Dimensions{axis}']")
        if node is None:
            raise ValueError(f'Warp settings {settingsFile} have no Tomo Dimensions{axis}')
        dims.append(int(float(node.get('Value'))))
    return tuple(dims)
