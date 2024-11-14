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
from pwem.emlib.image.image_readers import ImageStack, ImageReadersRegistry
from pyworkflow import BETA
import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from pyworkflow.object import Set
from tomo.objects import SetOfTiltSeriesM, SetOfTiltSeries, TiltImage, TiltSeries
from tomo.protocols import ProtTomoBase
from .protocol_base import ProtWarpBase

from .. import (Plugin, CREATE_SETTINGS, FS_MOTION, FRAMESERIES_FOLDER, FRAMESERIES_SETTINGS, AVERAGE_FOLDER,
                OUTPUT_TILTSERIES)


class ProtWarpTSMotionCorr(ProtWarpBase, ProtTomoBase):
    """ This protocol wraps WarpTools programs.
        Align tilt-series movies
    """

    _label = 'align tilt-series movies'
    _devStatus = BETA
    evenOddCapable = True
    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputTSMovies', params.PointerParam, pointerClass=SetOfTiltSeriesM,
                      important=True,
                      label=pwutils.Message.LABEL_INPUT_MOVS,
                      help='Select a set of previously imported tilt series movies.')
        form.addSection('Alignment')
        self._defineAlignmentParams(form)

    def _defineAlignmentParams(self, form):
        form.addHidden(params.GPU_LIST, params.StringParam, default='0',
                       expertLevel=params.LEVEL_ADVANCED,
                       label="Choose GPU IDs",
                       help="Space-separated list of GPU IDs to use for processing. Default: all GPUs in the system."
                            " Warp can use multiple GPUs - in that case"
                            " set to i.e. *0 1 2*.")

        form.addParam('binFactor', params.IntParam, default=1,
                      label="Binning factor",
                      help="Binning factor, applied in Fourier "
                           "space when loading raw data. 1 = no binning, "
                           "2 = 2x2 binning, 4 = 4x4 binning, supports "
                           "non-integer values")

        line = form.addLine('Resolution to fit',
                            help='Resolution in Angstrom to consider in fit.')
        line.addParam('range_min', params.IntParam, default=500,
                      label='Min')
        line.addParam('range_max', params.IntParam, default=10,
                      label='Max')

        form.addParam('bfactor', params.IntParam, default=-500,
                      label="B-factor",
                      help="Downweight higher spatial frequencies using a "
                           "B-factor, in Angstrom^2")

        line = form.addLine('Motion model grid',
                            help="Resolution of the motion model grid in "
                                 "X, Y, and temporal dimensions, e.g. 5x5x40; "
                                 "0 = auto")
        line.addParam('x', params.IntParam, default=2, label='X')
        line.addParam('y', params.IntParam, default=2, label='Y')
        line.addParam('z', params.IntParam, default=1, label='Temporal')

        form.addParam('average_halves', params.BooleanParam,
                      default=False,
                      label='Do even and odd ?',
                      help='Export aligned averages of odd and even frames separately, e.g. for denoiser training')

        form.addSection(label="Gain and defects")
        form.addParam('gainSwap', params.EnumParam,
                      choices=['no swap', 'transpose X/Y'],
                      label="Transpose gain reference:",
                      default=0,
                      display=params.EnumParam.DISPLAY_COMBO)

        form.addParam('gainFlip', params.EnumParam,
                      choices=['no flip', 'flip X', 'flip Y'],
                      label="Flip gain reference:", default=0,
                      display=params.EnumParam.DISPLAY_COMBO)

    # --------------------------- STEPS functions -----------------------------
    def _insertAllSteps(self):
        inputTSMovies = self.inputTSMovies.get()
        self.samplingRate = inputTSMovies.getSamplingRate()
        self._insertFunctionStep(self.createSettingStep, needsGPU=False)
        self._insertFunctionStep(self.proccessMoviesStep,  needsGPU=True)

    def createSettingStep(self):
        """ Create a settings file. """
        tsMovies = self.inputTSMovies.get()
        firstTSMovie = tsMovies.getFirstItem()
        fileName, extension = os.path.splitext(firstTSMovie.getFirstItem().getFileName())
        folderData = os.path.abspath(os.path.dirname(fileName))
        processingFolder = os.path.abspath(self._getExtraPath(FRAMESERIES_FOLDER))
        exposure = tsMovies.getAcquisition().getDosePerFrame()
        gainPath = os.path.abspath(tsMovies.getGain()) if tsMovies.getGain() else None
        pwutils.makePath(processingFolder)
        argsDict = {
            "--folder_data": folderData,
            "--extension": "*%s" % extension,
            "--folder_processing": processingFolder,
            "--bin": self.getBinFactor(),
            "--angpix": self.samplingRate,
            "--exposure": exposure,
            "--output": os.path.abspath(self._getExtraPath(FRAMESERIES_SETTINGS)),
        }

        cmd = ' '.join(['%s %s' % (k, v) for k, v in argsDict.items()])
        if gainPath:
            cmd += " --gain_path %s" % gainPath
            if self.gainFlip.get() == 1:
                cmd += ' --gain_flip_x'
            elif self.gainFlip.get() == 2:
                cmd += ' --gain_flip_y'
            if self.gainSwap.get() == 1:
                cmd += ' --gain_transpose'

        self.runJob(Plugin.getProgram(CREATE_SETTINGS), cmd, executable='/bin/bash')

    def proccessMoviesStep(self) -> None:
        """Estimate motion in frame series, produce aligned averages and register the output"""
        inputTSMovies = self.inputTSMovies.get()
        for tsMovie in inputTSMovies.iterItems():
            tsId = tsMovie.getTsId()
            self.info(">>> Starting estimate motion for %s..." % tsId)
            # Prepare a list of absolute paths for the movies to process
            # Each movie name in micNamesList is converted to an absolute path and join them into a
            # single string separated by spaces (warp specification)
            warpMoviesNamesList = [os.path.abspath(tiName.getFileName()) for tiName in tsMovie.iterItems()]
            warpMoviesNamesList = " ".join(warpMoviesNamesList)
            outputProcessingFolder = os.path.abspath(os.path.join(self._getExtraPath(AVERAGE_FOLDER), tsId))
            argsDict = {
                "--settings": os.path.abspath(self._getExtraPath(FRAMESERIES_SETTINGS)),
                "--range_min": self.range_min.get(),
                "--range_max": self.range_max.get(),
                "--bfac": self.bfactor.get(),
                "--input_data": warpMoviesNamesList,
                "--output_processing": outputProcessingFolder
            }
            gpuList = self.getGpuList()
            if gpuList:
                argsDict['--device_list'] = ' '.join(map(str, gpuList))

            cmd = ' '.join(['%s %s' % (k, v) for k, v in argsDict.items()])
            cmd += ' --averages'

            if self.average_halves.get():
                cmd += ' --average_halves'

            if self.x.get() and self.y.get() and self.z.get():
                cmd += ' --grid %sx%sx%s' % (self.x.get(), self.y.get(), self.z.get())

            self.runJob(self.getPlugin().getProgram(FS_MOTION), cmd, executable='/bin/bash')
            self.createOutput(tsMovie)

        self._closeOutputSet()

    def createOutput(self, tsMovie):
        self.info(">>> Generating output for %s..." % tsMovie.getTsId())
        output = self.getOutputSetOfTS('TiltSeries')
        averageFolder = os.path.join(self._getExtraPath(AVERAGE_FOLDER), tsMovie.getTsId(), AVERAGE_FOLDER)

        tsId = tsMovie.getTsId()
        newTs = TiltSeries(tsId=tsId)
        output.append(newTs)

        tiOrderDict = {}
        properties = {"sr": tsMovie.getSamplingRate()}
        newStack = ImageStack(properties=properties)
        oddFileNames = ImageStack(properties=properties)
        evenFileNames = ImageStack(properties=properties)
        newBinaryName = os.path.join(averageFolder, tsId + '.mrcs')
        hasAverageHalves = self.average_halves.get()
        newOddBinaryName = os.path.join(averageFolder, 'odd', tsId + '_odd.mrcs')
        newEvenBinaryName = os.path.join(averageFolder, 'even', tsId + '_even.mrcs')

        for index, tiM in enumerate(tsMovie):
            fileName = os.path.splitext(os.path.basename(tiM.getFileName()))[0] + '.mrc'
            newTi = TiltImage(location=(index + 1, newBinaryName))
            newTi.copyInfo(tiM)
            newTi.setAcquisition(tiM.getAcquisition().clone())
            newTi.setSamplingRate(tiM.getSamplingRate() * self.binFactor.get())
            tiOrderDict[newTi.getTiltAngle()] = (newTi, fileName)
            if hasAverageHalves:
                newTi.setOddEven([newOddBinaryName, newEvenBinaryName])

        sortedTiltAngle = sorted(tiOrderDict.keys())

        for angle in sortedTiltAngle:
            fileName = tiOrderDict[angle][1]
            newStack.append(ImageReadersRegistry.open(os.path.join(averageFolder, fileName)))
            if hasAverageHalves:
                oddFileNames.append(ImageReadersRegistry.open(os.path.join(averageFolder, 'odd', fileName)))
                evenFileNames.append(ImageReadersRegistry.open(os.path.join(averageFolder, 'even', fileName)))

        ImageReadersRegistry.write(newStack, newBinaryName, isStack=True)
        if hasAverageHalves:
            ImageReadersRegistry.write(oddFileNames, newOddBinaryName, isStack=True)
            ImageReadersRegistry.write(evenFileNames, newEvenBinaryName, isStack=True)

        for index, angle in enumerate(sortedTiltAngle):
            ti = tiOrderDict[angle][0]
            ti.setIndex(index + 1)
            ti.setObjId(index + 1)
            newTs.append(ti)

        output.update(newTs)
        output.write()
        self._store(output)

    def _summary(self):
        summary = []
        if self.hasAttribute(OUTPUT_TILTSERIES):
            summary.append(f"Aligned tiltseries: {self.TiltSeries.getSize()} of {self.inputTSMovies.get().getSize()}\n")
        else:
            summary.append("Outputs are not ready yet.")
        return summary

    def getOutputSetOfTS(self, outputSetName):
        outputSetOfTiltSeries = getattr(self, outputSetName, None)

        if outputSetOfTiltSeries:
            outputSetOfTiltSeries.enableAppend()
        else:
            outputSetOfTiltSeries = SetOfTiltSeries.create(self._getPath(), template='tiltseries')
            tsMovieSet = self.inputTSMovies.get()
            outputSetOfTiltSeries.setSamplingRate(tsMovieSet.getSamplingRate() * self.binFactor.get())
            outputSetOfTiltSeries.setAcquisition(tsMovieSet.getAcquisition())
            outputSetOfTiltSeries.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(**{outputSetName: outputSetOfTiltSeries})
            self._defineSourceRelation(outputSetOfTiltSeries, tsMovieSet)

        return outputSetOfTiltSeries

    def getBinFactor(self):
        import math
        return math.floor(math.log2(self.binFactor.get()))


















