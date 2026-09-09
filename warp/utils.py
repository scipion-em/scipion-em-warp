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

import numpy as np
import math
import scipy.fft
import xml.etree.ElementTree as eTree

from pwem.convert import transformations
from pwem.convert.transformations import translation_from_matrix, euler_from_matrix
from tomo.constants import TR_RELION
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
    term1 = lambda1**3 * cs * k2**2

    w = np.pi/2 * (term1 + lambda2 * defocus * k2) - phaseshift

    acurve = np.cos(w) * amplitude
    pcurve = -np.sqrt(1 - amplitude**2) * np.sin(w)
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
    data = np.arange(0, 1+1/2047, 1/2047, dtype=np.float32)
    highpass = np.minimum(1, data/highpassnyquist) * np.pi
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
    wiener = ctf / (ctf*ctf + 1.0/snr)

    s1, f1 = -math.floor(vol.shape[0]/2), -math.floor(vol.shape[0]/2) + vol.shape[0] - 1
    m1 = np.arange(s1, f1 + 1, dtype=np.int32)
    s2, f2 = -math.floor(vol.shape[1]/2), -math.floor(vol.shape[1]/2) + vol.shape[1] - 1
    m2 = np.arange(s2, f2 + 1, dtype=np.int32)

    if vol.ndim == 3:
        s3, f3 = -math.floor(vol.shape[2]/2), -math.floor(vol.shape[2]/2) + vol.shape[2] - 1
        m3 = np.arange(s3, f3 + 1, dtype=np.int32)
        x, y, z = np.meshgrid(m1, m2, m3, indexing='ij')
        x = np.divide(x, abs(s1))
        y = np.divide(y, abs(s2))
        z = np.divide(z, max(1, abs(s3)))
        r = np.sqrt(x**2 + y**2 + z**2)
    else:
        x, y = np.meshgrid(m1, m2, indexing='ij')
        x = np.divide(x, abs(s1))
        y = np.divide(y, abs(s2))
        r = np.sqrt(x**2 + y**2)

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


def parseTsDefocusHandFile(xmlPath):
    """Parse the .xml file generated by ts_defocus_hand command containing the ctf_hand value"""
    tree = eTree.parse(xmlPath)
    root = tree.getroot()
    ctfParams = root.find('CTF')


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


def getTransformInfoFromCoordOrSubtomo(coord, samplingRate):
    M = coord.getMatrix(convention=TR_RELION)
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
