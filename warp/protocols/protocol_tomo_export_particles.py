# ******************************************************************************
# *
# * Authors:     Yunior C. Fonseca Reyna
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
# ******************************************************************************

import os
from enum import Enum

from emtable import Table

from pyworkflow import BETA
import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from pyworkflow.object import Integer
from pyworkflow.utils import Message
from reliontomo.convert import readSetOfPseudoSubtomograms
from reliontomo.objects import createSetOfRelionPSubtomograms, RelionSetOfPseudoSubtomograms
from tomo.constants import BOTTOM_LEFT_CORNER

from warp.constants import *
from warp.protocols.protocol_base import ProtWarpBase
from warp.utils import updateCtFXMLFile, getTransformInfoFromCoordOrSubtomo


class outputObjects(Enum):
    relionParticles = RelionSetOfPseudoSubtomograms


class ProtWarpExportParticles(ProtWarpBase):
    """
    Export particles as 3D volumes or 2D image series
    More info:
       https://warpem.github.io/warp/user_guide/warptools/quick_start_warptools_tilt_series/#export-particles
    """

    _label = 'Export particles'
    _devStatus = BETA

    def __init__(self, **kwargs):
        ProtWarpBase.__init__(self, **kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('coordinates', params.PointerParam,
                      pointerClass='SetOfCoordinates3D',
                      label="3D Coordinates",
                      important=True,
                      allowsNull=False,
                      help='3D coordinates')
        form.addParam('inputSet', params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      label="Tilt series",
                      important=True,
                      help='Tilt series with alignment (non interpolated) used in the tomograms reconstruction.')
        form.addParam('inputSetOfCtfTomoSeries',
                      params.PointerParam,
                      label="Input CTF estimation",
                      pointerClass='SetOfCTFTomoSeries',
                      help='Select the CTF estimation for the set '
                           'of tilt-series.')

        form.addSection(label='Reconstruct')
        form.addParam('output_angpix', params.FloatParam,
                      label='Output pixel size',
                      important=True,
                      default=None,
                      help='Pixel size at which to export particles')

        form.addParam('box', params.IntParam,
                      label='Output box size',
                      important=True,
                      default=None,
                      help='Output has this many pixels/voxels on each side')

        form.addParam('diameter', params.IntParam,
                      label='Particles diameter (Å)',
                      important=True,
                      default=None,
                      help='Particle diameter in angstroms')

        form.addParam('writeStacks', params.EnumParam,
                      label='Export type',
                      default=0,
                      choices=['2D', '3D'],
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='If set to 2D, this program will write output particles as 2d image series centered on '
                           'the particle (particle series). If set to 3D, this program will write output '
                           'particles as 3d images (subtomograms)')

    def _insertAllSteps(self):
        self._insertFunctionStep(self.prepareDataStep, needsGPU=True)
        self._insertFunctionStep(self.exportParticlesStep, needsGPU=True)
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)
        self._insertFunctionStep(self.cleanIntermediateResults, needsGPU=False)


    def prepareDataStep(self):
        inputTs = self.inputSet.get()
        coordSet = self.coordinates.get()

        tsSr = inputTs.getSamplingRate()
        tomoSr = coordSet.getSamplingRate()
        tomo_dim = coordSet.getPrecedents().getDim()

        scaleFactor = tomoSr / tsSr
        self.tomo_thickness = Integer(round(tomo_dim[2] * scaleFactor))
        self.x_dimension = Integer(round(tomo_dim[0] * scaleFactor))
        self.y_dimension = Integer(round(tomo_dim[1] * scaleFactor))

        # Create the warp settings file and retrieve ctf and alignments values
        self.info(">>> Creating the warp settings file and retrieve ctf and alignments values...")
        self.createTiltSeriesSetting(None)
        for ts in inputTs.iterItems(iterate=False):
            self.tsDataPrepare(ts)
            self.createImodFiles(ts)
        self.tsCtfEstimation()
        self.updateCTFValues()
        self.tsImportAligments()

        sRate = coordSet.getSamplingRate()
        outPath = self._getExtraPath(MATCHING_FOLDER)
        pwutils.makePath(outPath)

        currentTomoId = None
        particlesTable = None
        f = None

        # Generate .star file per tomogram
        self.info(">>> Generate .star file per tomogram...")
        for coord in coordSet.iterCoordinates(orderBy='_tomoId'):
            tomoId = coord.getTomoId()

            if tomoId != currentTomoId:
                if f is not None:
                    particlesTable.writeStar(f, tableName='')
                    f.close()

                currentTomoId = tomoId
                tomoDim = coord.getVolume().getDim()
                filePath = os.path.join(outPath, f'{currentTomoId}.star')
                f = open(filePath, 'w')
                particlesTable = Table(columns=tomoStarFields)

            angles, _ = getTransformInfoFromCoordOrSubtomo(coord, sRate)
            particlesTable.addRow(
                coord.getX(BOTTOM_LEFT_CORNER) / tomoDim[0],
                coord.getY(BOTTOM_LEFT_CORNER) / tomoDim[1],
                coord.getZ(BOTTOM_LEFT_CORNER) / tomoDim[2],
                angles[0],
                angles[1],
                angles[2],
                f'{currentTomoId}.tomostar',
                coord.getScore()
            )
        # Close the last open file
        if f is not None:
            particlesTable.writeStar(f, tableName='')
            f.close()

    def exportParticlesStep(self):
        self.info(">>> Exporting particles...")
        settingFile = self._getExtraPath(TILTSERIE_SETTINGS)
        matchinFolder = self._getExtraPath(MATCHING_FOLDER)
        output = self._getExtraPath(RELION_FOLDER)
        pwutils.makePath(output)
        argsDict = {
            "--settings": os.path.abspath(settingFile),
            "--input_directory": matchinFolder,
            "--input_pattern": "*.star",
            "--output_star": os.path.join(output, 'matching.star'),
            "--output_angpix": self.output_angpix.get(),
            "--box": self.box.get(),
            "--diameter": self.diameter.get(),
        }
        cmd = '--relative_output_paths --normalized_coords'
        if self.writeStacks.get() == 0:
            cmd += ' --2d'
        else:
            cmd += ' --3d'
        self.runProgram(argsDict, WARP_TOOLS, TS_EXPORT_PARTICLES, othersCmds=cmd)

    def createOutputStep(self):
        coordSet = self.coordinates.get()
        tsSet = self.inputSet.get()
        tsSRate = tsSet.getSamplingRate()
        boxSize = self.box.get()
        acq = tsSet.getAcquisition()
        relionFolder = self._getExtraPath(RELION_FOLDER)
        are2dStacks = self.writeStacks.get() == 0

        psubtomoSet = createSetOfRelionPSubtomograms(self._getPath(),
                                                     os.path.join(relionFolder, OPTIMISATION_SET_STAR),
                                                     coordSet,
                                                     template='pseudosubtomograms%s.sqlite',
                                                     tsSamplingRate=tsSRate,
                                                     relionBinning=self.output_angpix.get()/tsSet.getSamplingRate(),
                                                     boxSize=boxSize,
                                                     are2dStacks=are2dStacks,
                                                     acquisition=acq)

        self.modifyStarFileMultiTable(os.path.join(relionFolder, MATCHING_PARTICLES_STAR),
                                      '_rlnImageName', lambda v: self.normalizeParticlesPath(v))
        self.modifyStarFileMultiTable(os.path.join(relionFolder, MATCHING_TOMOGRAMS_STAR),
                                      '_rlnTomoTiltSeriesName', lambda v: self.normalizeTomogramsPath(v))
        # Fill the set with the generated particles
        readSetOfPseudoSubtomograms(psubtomoSet)
        outDict = {outputObjects.relionParticles.name: psubtomoSet}
        self._defineOutputs(**outDict)
        self._defineSourceRelation(self.coordinates, psubtomoSet)
        self._defineSourceRelation(self.inputSetOfCtfTomoSeries, psubtomoSet)
        self._defineSourceRelation(self.inputSet, psubtomoSet)

    def normalizeParticlesPath(self, particleValue):
        particleFile = particleValue
        if not os.path.exists(particleFile):
            particleFile = os.path.normpath(os.path.join(self._getExtraPath(TILTSERIES_FOLDER), particleValue))
        return particleFile

    def normalizeTomogramsPath(self, tomogramValue):
        tomogramFile = tomogramValue
        if not os.path.exists(tomogramFile):
            tomogramFile = os.path.join(self._getExtraPath(RELION_FOLDER), tomogramValue)
        return tomogramFile

    def tsCtfEstimation(self):
        """CTF estimation"""

        self.info(">>> Generating ctf estimation...")
        settingFile = self._getExtraPath(TILTSERIE_SETTINGS)
        argsDict = {
            "--settings": os.path.abspath(settingFile)
        }
        try:
            self.runProgram(argsDict, WARP_TOOLS, TS_CTF)
        except Exception:
            self.info(">>> Error generating ctf estimation file...")

    def updateCTFValues(self):
        for ts in self.inputSet.get().iterItems():
            ctfTomoSeries = self.inputSetOfCtfTomoSeries.get().getItem('_tsId', ts.getTsId())
            processingFolder = os.path.abspath(self._getExtraPath(TILTSERIES_FOLDER))
            defocusFilePath = os.path.join(processingFolder, ts.getTsId() + '.xml')
            updateCtFXMLFile(defocusFilePath, ctfTomoSeries)

    def tsImportAligments(self):
        self.info(">>> Starting import alignments...")
        processingFolder = os.path.abspath(self._getExtraPath(TILTSERIES_FOLDER))
        tiltstackFolder = os.path.join(processingFolder, 'tiltstack')
        angpix = self.inputSet.get().getSamplingRate()
        settingFile = self._getExtraPath(TILTSERIE_SETTINGS)
        argsDict = {
            "--settings": os.path.abspath(settingFile),
            '--alignments': os.path.abspath(tiltstackFolder),
            "--alignment_angpix": angpix,
        }
        cmd = ' '.join(['%s %s' % (k, v) for k, v in argsDict.items()])
        self.runJob(self.getPlugin().getProgram(WARP_TOOLS, TS_IMPORT_ALIGNMENTS), cmd, executable='/bin/bash')

    def createImodFiles(self, ts):
        processingFolder = os.path.abspath(self._getExtraPath(TILTSERIES_FOLDER))
        tiltstackFolder = os.path.join(processingFolder, 'tiltstack', ts.getTsId())
        pwutils.makePath(tiltstackFolder)
        ts.writeImodFiles(tiltstackFolder, delimiter=' ', factor=1)

    def modifyStarFileMultiTable(self, filePath, columnToModify, modifierFunc):
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

    def cleanIntermediateResults(self):
        self.info(">>> Cleaning intermediate results...")
        imagesFolder = self._getExtraPath(TILTIMAGES_FOLDER)
        pwutils.cleanPath(imagesFolder)




