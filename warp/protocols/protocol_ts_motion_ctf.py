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
from warp.protocols.protocol_base import ProtWarpBase, ProtTSMovieAlignBase
from warp.constants import *
from warp.utils import parseCtfXMLFile, tomoStarGenerate


class ProtWarpTSMotionCorr(ProtTomoBase, ProtTSMovieAlignBase):
    """ This protocol wraps WarpTools programs.
        Estimate motion in frame series, produce aligned averages, estimate CTF
    """

    _label = 'tilt-series motion and ctf estimation'
    _devStatus = BETA
    evenOddCapable = True

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.averageCorrelation = Float()

    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputTSMovies', params.PointerParam, pointerClass=SetOfTiltSeriesM,
                      important=True,
                      label=pwutils.Message.LABEL_INPUT_MOVS,
                      help='Select a set of previously imported tilt series movies.')
        form.addSection('Alignment')
        self._defineAlignmentParams(form)
        ProtTSMovieAlignBase._defineStreamingParams(self, form)

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

        form.addSection("EER")
        form.addParam('EERtext', params.LabelParam,
                      label="These options are ignored for non-EER movies.")
        form.addParam('eer_ngroups', params.IntParam, default=16,
                      allowsNull=True,
                      label='EER fractionation',
                      help="Number of groups to combine raw EER frames into, i.e. number of 'virtual' "
                           "frames in resulting stack; use negative value to specify the number of "
                           "frames per virtual frame instead")
        form.addParam('eer_groupexposure', params.FloatParam, default=None,
                      allowsNull=True,
                      label='EER group exposure',
                      help="As an alternative to EER fractionation, fractionate the frames so that a group will "
                           "have this exposure in e-/A^2; this overrides EER fractionation"
                           "\nFractionate such that each fraction "
                           "has about 0.5 to 1.25 e/A2.")
        form.addSection(label="CTF")

        form.addParam('estimateCTF', params.BooleanParam, default=True,
                      label='Estimate the CTF ?',
                      help='Estimate the CTF')

        form.addParam('window', params.IntParam, default=512,
                      condition='estimateCTF',
                      label='Windows', help='Patch size for CTF estimation in binned pixels')

        line = form.addLine('Resolution (Å)',
                            condition='estimateCTF',
                            help='Resolution in Angstrom to consider in fit.')

        line.addParam('range_min', params.FloatParam, default=30,
                      condition='estimateCTF',
                      label='Min', help='Lowest (worst) resolution in Angstrom to consider in fit')

        line.addParam('range_max', params.FloatParam, default=4,
                      condition='estimateCTF',
                      label="Max",
                      help="Highest (best) resolution in Angstrom to consider in fit")

        line = form.addLine('Defocus search range (um)',
                            condition='estimateCTF',
                            help='Defocus values in um to explore during fitting (positive = underfocus). '
                                 'The units are microns!!')
        line.addParam('defocus_min', params.FloatParam, default=0.5,
                      condition='estimateCTF',
                      label='Min', help='Minimum defocus value in um to explore during fitting (positive = underfocus)')
        line.addParam('defocus_max', params.FloatParam, default=5,
                      condition='estimateCTF',
                      label='Max', help='Maximum defocus value in um to explore during fitting (positive = underfocus)')

        line = form.addLine('Defocus model grid',
                            condition='estimateCTF',
                            help="Resolution of the defocus model grid in X, Y, and temporal dimensions, " 
                                 "separated by x: e.g. 5x5x40; empty = auto; Z > 1 is purely experimental")

        line.addParam('c_x', params.IntParam, default=None,
                      condition='estimateCTF',
                      allowsNull=True, label='X')
        line.addParam('c_y', params.IntParam, default=None,
                      condition='estimateCTF',
                      allowsNull=True, label='Y')
        line.addParam('c_z', params.IntParam, default=None, allowsNull=True,
                      condition='estimateCTF',
                      label='Temporal')

        form.addParam('fit_phase', params.BooleanParam, default=False,
                      condition='estimateCTF',
                      label='Fit phase', help='Fit the phase shift of a phase plate')

        form.addParam('use_sum', params.BooleanParam, default=False,
                      condition='estimateCTF',
                      label='Use the movie average',
                      help='Use the movie average spectrum instead of the average of individual '
                           'frames spectra. Can help in the absence of an energy filter, or when signal is low')

        form.addParam('handedness', params.BooleanParam, default=False,
                      condition='estimateCTF',
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Check the handedness ?',
                      help='Checking defocus handedness across a dataset ')
        form.addParallelSection(threads=2, mpi=0)

    # --------------------------- STEPS functions -----------------------------

    def insertInitialSteps(self):
        self.samplingRate = self.getInputTSMovies().getSamplingRate()
        createSettingStep = self._insertFunctionStep(self.createFrameSeriesSettingStep,
                                                     prerequisites=[], needsGPU=False)
        return [createSettingStep]

    def createFrameSeriesSettingStep(self):
        """ Create a settings file. """
        self.info(">>> Starting frame series settings creation...")
        tsMovies = self.getInputTSMovies()
        firstTSMovie = tsMovies.getFirstItem()
        fileName, extension = os.path.splitext(firstTSMovie.getFirstItem().getFileName())
        folderData = os.path.abspath(os.path.dirname(fileName))
        processingFolder = os.path.abspath(self._getExtraPath(FRAMESERIES_FOLDER))
        exposure = -1 * tsMovies.getAcquisition().getDosePerFrame()
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

        if extension == '.eer':
            if self.eer_ngroups.get() is not None:
                argsDict['--eer_ngroups'] = self.eer_ngroups.get()
            if self.eer_groupexposure.get():
                argsDict['--eer_groupexposure'] = self.eer_groupexposure.get()

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

    def createTiltSeriesSettingStep(self, tsId):
        self.info(">>> Starting tilt-series settings creation (%s)..." % tsId)
        setOfTSMovies = self.inputTSMovies.get()
        sr = setOfTSMovies.getSamplingRate()
        exposure = setOfTSMovies.getAcquisition().getDosePerFrame()
        firstTSMovie = setOfTSMovies.getFirstItem()
        fileName, extension = os.path.splitext(firstTSMovie.getFirstItem().getFileName())
        settingsFolder = os.path.abspath(self._getExtraPath(SETTINGS_FOLDER))
        pwutils.makePath(settingsFolder)
        processingFolder = os.path.abspath(self._getExtraPath(TILTSERIES_FOLDER))
        pwutils.makePath(processingFolder)
        tsSettingFile = tsId + '_' + TILTSERIE_SETTINGS
        tsSettingFilePath = os.path.abspath(os.path.join(self._getExtraPath(settingsFolder), tsSettingFile))
        argsDict = {
            "--folder_data": os.path.abspath(self._getExtraPath(TOMOSTAR_FOLDER)),
            "--extension": "%s.tomostar" % tsId,
            "--folder_processing": processingFolder,
            "--bin": self.getBinFactor(),
            '--angpix': sr,
            "--output": tsSettingFilePath
        }

        if exposure is not None:
            argsDict['--exposure'] = -1 * exposure

        if hasattr(self, 'tomo_thickness'):
            z = self.tomo_thickness.get()
            x = self.x_dimension.get() or setOfTSMovies.getDimensions()[0]
            y = self.y_dimension.get() or setOfTSMovies.getDimensions()[1]

            argsDict['--tomo_dimensions'] = f'{x}x{y}x{z}'

        if extension == '.eer':
            argsDict['--eer_ngroups'] = self.eer_ngroups.get()
            if self.eer_groupexposure.get():
                argsDict['--eer_groupexposure'] = self.eer_groupexposure.get()

        cmd = ' '.join(['%s %s' % (k, v) for k, v in argsDict.items()])

        self.runJob(Plugin.getProgram(CREATE_SETTINGS), cmd, executable='/bin/bash')

    def dataPrepare(self, tsMovie):
        """Creates the setting file that will be used by the different programs.
           It also extracts the tiltimages from the tiltseries and generates the *.tomostar files based on
           the tiltimages."""
        starFolder = self._getExtraPath(TOMOSTAR_FOLDER)
        pwutils.makePath(starFolder)
        imagesFolder = self._getExtraPath(FRAMES_FOLDER)
        invertTiltAngle = 1
        pwutils.makePath(imagesFolder)

        if tsMovie.isEnabled():
            tsId = tsMovie.getTsId()
            tiValues = {}
            for ti in tsMovie.iterItems():
                if ti.isEnabled():  # Excluding views
                    dose = 0
                    maskedFraction = 0
                    shiftX = 0
                    shiftY = 0
                    axisAngle = 0
                    amplitudeContrast = 0

                    if tsMovie.hasAcquisition():
                        axisAngle = tsMovie.getAcquisition().getTiltAxisAngle()
                    if ti.getAcquisition():
                        amplitudeContrast = ti.getAcquisition().getAmplitudeContrast()
                        dose = ti.getAcquisition().getAccumDose()
                    fileName = ti.getFileName()
                    newBinaryName = os.path.basename(fileName)
                    os.symlink(os.path.abspath(fileName), os.path.join(imagesFolder, os.path.basename(fileName)))

                    tiValues[ti.getTiltAngle() * invertTiltAngle] = [newBinaryName, ti.getTiltAngle() * invertTiltAngle,
                                                                     axisAngle, shiftX, shiftY, dose,
                                                                     amplitudeContrast, maskedFraction]

            tomoStarGenerate(tsId, tiValues, starFolder, 0)

    def tsDefocusHandStep(self):
        """Defocus handedness"""
        self.info(">>> Starting defocus handedness...")
        objSet = self.inputTSMovies.get()
        settingsFolder = os.path.abspath(self._getExtraPath(SETTINGS_FOLDER))
        tsId = objSet.getFirstItem().getTsId()
        tsSettingFile = tsId + '_' + TILTSERIE_SETTINGS
        tsSettingFilePath = os.path.abspath(os.path.join(self._getExtraPath(settingsFolder), tsSettingFile))
        argsDict = {
            "--settings": tsSettingFilePath,
        }
        cmd = ' '.join(['%s %s' % (k, v) for k, v in argsDict.items()])
        cmd += ' --check'
        self.runJob(self.getPlugin().getProgram(TS_DEFOCUS_HAND), cmd, executable='/bin/bash')
        self.createOutputDefocusHand()

    def proccessTSMoviesStep(self, tsId) -> None:
        """Estimate motion in frame series, produce aligned averages and register the output"""
        tsMovie = self.getInputTSMovies().getItem(TiltSeries.TS_ID_FIELD, tsId)
        warpMoviesNamesList = [os.path.abspath(tiName.getFileName()) for tiName in tsMovie.iterItems()]
        warpMoviesNamesList = " ".join(warpMoviesNamesList)
        self.info(">>> Starting estimate motion for %s..." % tsId)
        self.fsMotionAndCTF(tsMovie, warpMoviesNamesList)
        if self.estimateCTF.get():
            self.createTiltSeriesSettingStep(tsId)
            self.dataPrepare(tsMovie)
            self.tsCtfEstimationStep(tsId)
        with self._lock:
            self.createOutputTS(tsMovie)
            if self.estimateCTF.get():
                self.createOutputCTF(tsId)

    def insertFinalSteps(self, proccessTSMoviesSteps) -> list:
        """The final steps inserted into the protocol"""
        finalSteps = []
        if self.handedness.get():
            finalStep = self._insertFunctionStep(self.tsDefocusHandStep, prerequisites=proccessTSMoviesSteps,
                                                 needsGPU=True)
            finalSteps.append(finalStep)
        return finalSteps

    def fsMotionAndCTF(self, tsMovie, warpMoviesNamesList):
        # Prepare a list of absolute paths for the movies to process
        # Each movie name in micNamesList is converted to an absolute path and join them into a
        # single string separated by spaces (warp specification)
        self.info(">>> Starting align motion process (%s) ..." % tsMovie.getTsId())
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

    def tsCtfEstimationStep(self, tsId):
        """CTF estimation"""
        self.info(">>> Starting ctf estimation to %s" % tsId)
        inputTSAdquisition = self.inputTSMovies.get().getFirstItem().getAcquisition()
        settingFile = self._getExtraPath(SETTINGS_FOLDER, tsId + '_' + TILTSERIE_SETTINGS)
        argsDict = {
            "--settings": os.path.abspath(settingFile),
            "--window": self.window.get(),
            "--range_low": self.range_min.get(),
            "--range_high": self.range_max.get(),
            "--defocus_min": self.defocus_min.get(),
            "--defocus_max": self.defocus_max.get(),
            "--voltage": int(inputTSAdquisition.getVoltage()),
            "--cs": inputTSAdquisition.getSphericalAberration(),
            "--amplitude": inputTSAdquisition.getAmplitudeContrast(),
        }

        gpuList = self.getGpuList()
        if gpuList:
            argsDict['--device_list'] = ' '.join(map(str, gpuList))

        cmd = ' '.join(['%s %s' % (k, v) for k, v in argsDict.items()])
        self.runJob(self.getPlugin().getProgram(TS_CTF), cmd, executable='/bin/bash')

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

    def createOutputCTF(self, tsId):
        self.info(">>> Generating outputs to %s" % tsId)
        processingFolder = os.path.abspath(self._getExtraPath(TILTSERIES_FOLDER))
        tsSet = self.TiltSeries
        if tsSet:
            psdStack = os.path.join(processingFolder, POWERSPECTRUM_FOLDER, tsId + '.mrc')
            ts = self.TiltSeries.getItem(TiltSeries.TS_ID_FIELD, tsId)
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
                        newCTFTomo.setPsdFile(f"{index}@" + psdStack)
                        newCTFTomoSeries.append(newCTFTomo)
                        index += 1

                outputSetOfCTFTomoSeries.update(newCTFTomoSeries)
                outputSetOfCTFTomoSeries.write()
                self._store(outputSetOfCTFTomoSeries)

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
        else:
            self.averageCorrelation = Float()
        summary.append(f"Aligned tiltseries: {tilseriesSize} of {self.inputTSMovies.get().getSize()}")

        if self.hasAttribute(OUTPUT_CTF_SERIE):
            ctfSize = self.CTFTomoSeries.getSize()
        summary.append(f"CTF estimated: {ctfSize} of {self.inputTSMovies.get().getSize()}")

        if self.handedness.get():
            if self.averageCorrelation.get():
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
            tsSet = self.TiltSeries
            outputSetOfCTFTomoSeries.setSetOfTiltSeries(tsSet)
            outputSetOfCTFTomoSeries.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(**{outputSetName: outputSetOfCTFTomoSeries})
            self._defineCtfRelation(outputSetOfCTFTomoSeries, tsSet)

        return outputSetOfCTFTomoSeries

    def getBinFactor(self):
        import math
        return math.floor(math.log2(self.binFactor.get()))
