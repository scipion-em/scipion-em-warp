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
import time

from pwem.emlib.image.image_readers import ImageStack, ImageReadersRegistry, logger
from pyworkflow import BETA
import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from pyworkflow.object import Set, Float, Boolean
from tomo.objects import (SetOfTiltSeriesM, SetOfTiltSeries, TiltImage,
                          TiltSeries, SetOfCTFTomoSeries, CTFTomoSeries,
                          CTFTomo)
from tomo.protocols import ProtTomoBase

from warp import Plugin
from warp.protocols.protocol_base import ProtWarpBase
from warp.constants import *
from warp.utils import parseCtfXMLFile


class ProtWarpTSMotionCorr(ProtWarpBase, ProtTomoBase):
    """ This protocol wraps WarpTools programs.
        Estimate motion in frame series, produce aligned averages, estimate CTF
    """

    _label = 'tilt-series motion and ctf estimation'
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

        form.addParam('binFactor', params.FloatParam, default=1,
                      label="Binning factor",
                      help="Binning factor, applied in Fourier "
                           "space when loading raw data. 1 = no binning, "
                           "2 = 2x2 binning, 4 = 4x4 binning, supports "
                           "non-integer values")

        line = form.addLine('Resolution to fit',
                            help='Resolution in Angstrom to consider in fit.')
        line.addParam('m_range_min', params.FloatParam, default=500,
                      label='Min', help='Minimum resolution in Angstrom to consider in fit')
        line.addParam('m_range_max', params.FloatParam, default=10,
                      label='Max', help='Maximun resolution in Angstrom to consider in fit')

        form.addParam('bfactor', params.FloatParam, default=-500,
                      label="B-factor",
                      help="Downweight higher spatial frequencies using a "
                           "B-factor, in Angstrom^2")

        line = form.addLine('Motion model grid',
                            help="Resolution of the motion model grid in X, Y, and temporal dimensions, "
                                 "separated by 'x': e.g. 5x5x40; empty = auto")
        line.addParam('x', params.IntParam, default=None,
                      allowsNull=True,
                      label='X')
        line.addParam('y', params.IntParam,
                      default=None,
                      allowsNull=True,
                      label='Y')
        line.addParam('z', params.IntParam, default=None,
                      allowsNull=True,
                      label='Temporal')

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

        form.addSection(label="CTF")
        form.addParam('window', params.IntParam, default=512,
                      label='Windows', help='Patch size for CTF estimation in binned pixels')

        line = form.addLine('Resolution (Å)',
                            help='Resolution in Angstrom to consider in fit.')

        line.addParam('range_min', params.FloatParam, default=30,
                      label='Min', help='Lowest (worst) resolution in Angstrom to consider in fit')

        line.addParam('range_max', params.FloatParam, default=4,
                      label="Max",
                      help="Highest (best) resolution in Angstrom to consider in fit")

        line = form.addLine('Defocus search range (um)',
                            help='Defocus values in um to explore during fitting (positive = underfocus). '
                                 'The units are microns!!')
        line.addParam('defocus_min', params.FloatParam, default=0.5,
                      label='Min', help='Minimum defocus value in um to explore during fitting (positive = underfocus)')
        line.addParam('defocus_max', params.FloatParam, default=5,
                      label='Max', help='Maximum defocus value in um to explore during fitting (positive = underfocus)')

        line = form.addLine('Defocus model grid',
                            help="Resolution of the defocus model grid in X, Y, and temporal dimensions, " 
                                 "separated by x: e.g. 5x5x40; empty = auto; Z > 1 is purely experimental")

        line.addParam('c_x', params.IntParam, default=None,
                      allowsNull=True, label='X')
        line.addParam('c_y', params.IntParam, default=None,
                      allowsNull=True, label='Y')
        line.addParam('c_z', params.IntParam, default=None, allowsNull=True,
                      label='Temporal')

        form.addParam('fit_phase', params.BooleanParam, default=False,
                      label='Fit phase', help='Fit the phase shift of a phase plate')

        form.addParam('use_sum', params.BooleanParam, default=False,
                      label='Use the movie average',
                      help='Use the movie average spectrum instead of the average of individual '
                           'frames spectra. Can help in the absence of an energy filter, or when signal is low')

        form.addParam('handedness', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Check the handedness ?',
                      help='Checking defocus handedness across a dataset ')

    # --------------------------- STEPS functions -----------------------------
    def _insertAllSteps(self):
        self.averageCorrelation = Float()
        inputTSMovies = self.inputTSMovies.get()
        self.samplingRate = inputTSMovies.getSamplingRate()
        self._insertFunctionStep(self.createFrameSeriesSettingStep, needsGPU=False)
        self._insertFunctionStep(self.createTiltSeriesSettingStep, needsGPU=False)
        self._insertFunctionStep(self.dataPrepare, inputTSMovies, needsGPU=False)
        self._insertFunctionStep(self.proccessMoviesStep,  needsGPU=True)
        self._insertFunctionStep(self.tsCtfEstimationStep, needsGPU=True)
        if self.handedness.get():
            self._insertFunctionStep(self.tsDefocusHandStep, needsGPU=True)
        self._insertFunctionStep(self.deleteIntermediateOutputsStep, needsGPU=False)

    def createFrameSeriesSettingStep(self):
        """ Create a settings file. """
        self.info(">>> Starting frame series settings creation...")
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
            "--extension": "'*%s'" % extension,
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

    def createTiltSeriesSettingStep(self):
        self.info(">>> Starting tilt-series settings creation...")
        objSet = self.inputTSMovies.get()
        sr = objSet.getSamplingRate()
        exposure = objSet.getAcquisition().getDosePerFrame()
        processingFolder = os.path.abspath(self._getExtraPath(TILTSERIES_FOLDER))
        pwutils.makePath(processingFolder)
        argsDict = {
            "--folder_data": os.path.abspath(self._getExtraPath(TOMOSTAR_FOLDER)),
            "--extension": "*.tomostar",
            "--folder_processing": processingFolder,
            "--bin": self.getBinFactor(),
            '--angpix': sr,
            "--output": os.path.abspath(self._getExtraPath(TILTSERIE_SETTINGS))
        }

        if exposure is not None:
            argsDict['--exposure'] = exposure

        if hasattr(self, 'tomo_thickness'):
            z = self.tomo_thickness.get()
            x = self.x_dimension.get() or objSet.getDimensions()[0]
            y = self.y_dimension.get() or objSet.getDimensions()[1]

            argsDict['--tomo_dimensions'] = f'{x}x{y}x{z}'

        cmd = ' '.join(['%s %s' % (k, v) for k, v in argsDict.items()])

        self.runJob(Plugin.getProgram(CREATE_SETTINGS), cmd, executable='/bin/bash')

    def tsDefocusHandStep(self):
        """Defocus handedness"""
        self.info(">>> Starting defocus handedness...")
        argsDict = {
            "--settings": os.path.abspath(self._getExtraPath(TILTSERIE_SETTINGS)),
        }
        cmd = ' '.join(['%s %s' % (k, v) for k, v in argsDict.items()])
        cmd += ' --check'
        self.runJob(self.getPlugin().getProgram(TS_DEFOCUS_HAND), cmd, executable='/bin/bash')
        self.createOutputDefocusHand()

    def proccessMoviesStep(self) -> None:
        """Estimate motion in frame series, produce aligned averages and register the output"""
        inputTSMovies = self.inputTSMovies.get()
        for tsMovie in inputTSMovies.iterItems():
            tsId = tsMovie.getTsId()
            warpMoviesNamesList = [os.path.abspath(tiName.getFileName()) for tiName in tsMovie.iterItems()]
            warpMoviesNamesList = " ".join(warpMoviesNamesList)
            self.info(">>> Starting estimate motion for %s..." % tsId)
            self.fsMotionAndCTF(tsMovie, warpMoviesNamesList)
            self.createOutputTS(tsMovie)

        self._closeOutputSet()

    def fsMotionAndCTF(self, tsMovie, warpMoviesNamesList):
        # Prepare a list of absolute paths for the movies to process
        # Each movie name in micNamesList is converted to an absolute path and join them into a
        # single string separated by spaces (warp specification)
        self.info(">>> Starting align motion process...")
        inputTSAdquisition = tsMovie.getFirstItem().getAcquisition()
        outputProcessingFolder = os.path.abspath(os.path.join(self._getExtraPath(FRAMESERIES_FOLDER)))
        argsDict = {
            "--settings": os.path.abspath(self._getExtraPath(FRAMESERIES_SETTINGS)),
            "--m_range_min": self.m_range_min.get(),
            "--m_range_max": self.m_range_max.get(),
            "--m_bfac": self.bfactor.get(),
            '--c_window': self.window.get(),
            '--c_range_min': self.range_min.get(),
            '--c_range_max': self.range_max.get(),
            '--c_defocus_min': self.defocus_min.get(),
            '--c_defocus_max': self.defocus_max.get(),
            "--c_voltage": int(inputTSAdquisition.getVoltage()),
            "--c_cs": inputTSAdquisition.getSphericalAberration(),
            "--c_amplitude": inputTSAdquisition.getAmplitudeContrast(),
            "--input_data": warpMoviesNamesList,
            "--output_processing": outputProcessingFolder
        }
        gpuList = self.getGpuList()
        if gpuList:
            argsDict['--device_list'] = ' '.join(map(str, gpuList))

        cmd = ' '.join(['%s %s' % (k, v) for k, v in argsDict.items()])
        cmd += ' --out_averages'

        if self.average_halves.get():
            cmd += ' --out_average_halves'

        if self.x.get() and self.y.get() and self.z.get():
            cmd += ' --m_grid %sx%sx%s' % (self.x.get(), self.y.get(), self.z.get())

        if self.c_x.get() and self.c_y.get() and self.c_z.get():
            cmd += ' --c_grid %sx%sx%s' % (self.c_x.get(), self.c_y.get(), self.c_z.get())

        if self.fit_phase.get():
            cmd += ' --c_fit_phase'
        if self.use_sum.get():
            cmd += ' --c_use_sum'

        self.runJob(self.getPlugin().getProgram(FS_MOTION_AND_CTF), cmd, executable='/bin/bash')

    def tsCtfEstimationStep(self):
        """CTF estimation"""
        self.info(">>> Starting ctf estimation...")
        inputTSAdquisition = self.inputTSMovies.get().getFirstItem().getAcquisition()
        argsDict = {
            "--settings": os.path.abspath(self._getExtraPath(TILTSERIE_SETTINGS)),
            "--window": self.window.get(),
            "--range_low": self.range_min.get(),
            "--range_high": self.range_max.get(),
            "--defocus_min": self.defocus_min.get(),
            "--defocus_max": self.defocus_max.get(),
            "--voltage": int(inputTSAdquisition.getVoltage()),
            "--cs": inputTSAdquisition.getSphericalAberration(),
            "--amplitude": inputTSAdquisition.getAmplitudeContrast(),
        }
        self.runProgram(argsDict, TS_CTF)
        self.createOutputCTF()

    def createOutputTS(self, tsMovie):
        self.info(">>> Generating output for %s..." % tsMovie.getTsId())
        outputTS = self.getOutputSetOfTS(OUTPUT_TILTSERIES)
        averageFolder = os.path.join(self._getExtraPath(FRAMESERIES_FOLDER), AVERAGE_FOLDER)

        tsId = tsMovie.getTsId()
        newTs = TiltSeries(tsId=tsId)
        outputTS.append(newTs)

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

        outputTS.update(newTs)
        outputTS.write()
        self._store(outputTS)

    def deleteIntermediateOutputsStep(self):
        try:
            averageFolder = os.path.join(self._getExtraPath(FRAMESERIES_FOLDER), AVERAGE_FOLDER)
            if not os.path.exists(averageFolder):
                logger.info(f"The directory {averageFolder} does not exist.")
                return
            for filename in os.listdir(averageFolder):
                if filename.endswith(".mrc") or filename.endswith(".json"):
                    file_path = os.path.join(averageFolder, filename)
                    os.remove(file_path)

            logger.info("All .mrc files have been deleted.")

        except Exception as e:
            logger.error(f"An error occurred: {e}")

    def createOutputCTF(self):
        self.info(">>> Generating outputs...")
        processingFolder = os.path.abspath(self._getExtraPath(TILTSERIES_FOLDER))
        tsSet = self.TiltSeries
        if tsSet:
            for ts in tsSet.iterItems():
                if ts.isEnabled():
                    tsId = ts.getTsId()
                    outputSetOfCTFTomoSeries = self.getOutputSetOfCTFTomoSeries(OUTPUT_CTF_SERIE)

                    # CTF outputs
                    newCTFTomoSeries = CTFTomoSeries(tsId=tsId)
                    newCTFTomoSeries.copyInfo(ts)
                    newCTFTomoSeries.setTiltSeries(ts)
                    outputSetOfCTFTomoSeries.append(newCTFTomoSeries)
                    defocusFilePath = os.path.join(processingFolder, ts.getTsId() + '.xml')
                    ctfData, gridCtfData = parseCtfXMLFile(defocusFilePath)
                    defocusDelta = float(ctfData['DefocusDelta']) * 1e4
                    defocusAngle = float(ctfData['DefocusAngle'])

                    index = 0
                    for ti in ts.iterItems():
                        if ti.isEnabled():
                            newCTFTomo = CTFTomo()
                            newCTFTomo.setAcquisitionOrder(ti.getAcquisitionOrder())
                            newCTFTomo.setIndex(index)
                            newCTFTomo.setObjId(index)
                            defocusU = 0
                            defocusV = 0
                            if index in gridCtfData["Nodes"]:
                                defocusU = gridCtfData["Nodes"][index] + defocusDelta
                                defocusV = gridCtfData["Nodes"][index] - defocusAngle
                            newCTFTomo.setDefocusU(defocusU)
                            newCTFTomo.setDefocusV(defocusV)
                            newCTFTomo.setDefocusAngle(defocusAngle)
                            newCTFTomo.setResolution(0)
                            newCTFTomo.setFitQuality(0)
                            newCTFTomo.standardize()
                            newCTFTomoSeries.append(newCTFTomo)
                            index += 1

                    outputSetOfCTFTomoSeries.update(newCTFTomoSeries)
                    outputSetOfCTFTomoSeries.write()
                    self._store(outputSetOfCTFTomoSeries)
            self._closeOutputSet()

    def createOutputDefocusHand(self):
        # Registering the output
        stdoutFile = os.path.abspath(os.path.join(self.getPath(), 'logs', 'run.stdout'))
        with open(stdoutFile, 'r', encoding='utf-8') as file:
            lines = file.readlines()

        for line in reversed(lines):
            if 'Average correlation:' in line:
                self.averageCorrelation.set(float(line.split()[-1]))
                outputAverage = Boolean(True)
                if self.averageCorrelation.get() < 0:
                    outputAverage = Boolean(False)

                self._defineOutputs(**{OUTPUT_HANDEDNESS: outputAverage})
                break
        self._store(self.averageCorrelation)

    def _summary(self):
        summary = []
        tilseriesSize = 0
        ctfSize = 0
        if self.hasAttribute(OUTPUT_TILTSERIES):
            tilseriesSize = self.TiltSeries.getSize()
        summary.append(f"Aligned tiltseries: {tilseriesSize} of {self.inputTSMovies.get().getSize()}")

        if self.hasAttribute(OUTPUT_CTF_SERIE):
            ctfSize = self.CTFTomoSeries.getSize()
        summary.append(f"CTF estimated: {ctfSize} of {self.inputTSMovies.get().getSize()}")

        if self.handedness.get():
            if self.hasAttribute('averageCorrelation') and self.averageCorrelation.get():
                # text = " (The average correlation is positive, which means that the defocus handedness should be set to '%s')"
                # flip = 'no flip'
                # if self.averageCorrelation.get() < 0:
                #     flip = 'flip'
                text = 'Warp convention is inverted related to ours (IMOD, Relion,...)'
                summary.append(f"Handedness: {self.averageCorrelation}  {text}")
            else:
                summary.append('Handedness: Not ready')

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

    def getOutputSetOfCTFTomoSeries(self, outputSetName):
        outputSetOfCTFTomoSeries = getattr(self, outputSetName, None)

        if outputSetOfCTFTomoSeries:
            outputSetOfCTFTomoSeries.enableAppend()
        else:
            outputSetOfCTFTomoSeries = SetOfCTFTomoSeries.create(self._getPath(),
                                                                 template='CTFmodels%s.sqlite')
            tsSet = self.inputTSMovies.get()
            outputSetOfCTFTomoSeries.setSetOfTiltSeries(self.inputTSMovies)
            outputSetOfCTFTomoSeries.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(**{outputSetName: outputSetOfCTFTomoSeries})
            self._defineCtfRelation(outputSetOfCTFTomoSeries, tsSet)

        return outputSetOfCTFTomoSeries

    def getBinFactor(self):
        import math
        return math.floor(math.log2(self.binFactor.get()))
