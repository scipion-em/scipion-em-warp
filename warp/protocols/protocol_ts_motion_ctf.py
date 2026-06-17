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
from pyworkflow.object import Set, Float, Boolean, Integer
from tomo.objects import (SetOfTiltSeriesM, SetOfTiltSeries, TiltImage,
                          TiltSeries, SetOfCTFTomoSeries, CTFTomoSeries,
                          CTFTomo)
from tomo.protocols import ProtTomoBase

from warp import Plugin
from warp.protocols.protocol_base import ProtTSMovieAlignBase
from warp.constants import *
from warp.utils import parseCtfXMLFile, tomoStarGenerate


class ProtWarpTSMotionCorr(ProtTomoBase, ProtTSMovieAlignBase):
    """
    Performs tilt-series movie motion correction and CTF estimation using WarpTools.
    The protocol aligns frame-series movies, generates aligned tilt-series averages,
    optionally estimates per-tilt CTF parameters, and evaluates defocus handedness.

    AI Generated:

    Tilt-Series Motion and CTF Estimation (ProtWarpTSMotionCorr) — User Manual
        Overview

        This protocol wraps external WarpTools programs to process cryo-electron
        tomography tilt-series movies. Its main goal is to correct beam-induced
        motion at the movie level, generate aligned tilt images, and optionally
        estimate CTF parameters for each tilt image.

        In a typical cryo-ET workflow, this protocol is used after importing raw
        tilt-series movies and before tomogram reconstruction. The output aligned
        tilt-series can be directly used for downstream alignment, reconstruction,
        or subtomogram analysis.

        General Workflow

        The protocol processes the input tilt-series in several stages:

        1. Frame-series settings creation
           A Warp settings file is generated using acquisition parameters such as
           sampling rate, exposure, binning, gain reference, and EER-specific options.

        2. Frame motion correction
           For each tilt-series, Warp estimates movie motion and generates aligned
           average images.

        3. Tilt-series preparation
           Temporary metadata and tomostar files are generated to describe the
           aligned tilt images.

        4. Optional CTF estimation
           If enabled, Warp estimates defocus parameters for each tilt image.

        5. Output registration
           The protocol assembles aligned tilt images into a new tilt-series stack
           and optionally creates a CTF series object.

        Input Parameters

        Input Tilt-Series Movies
            The protocol requires a previously imported set of tilt-series movies.

        Binning Factor
            Controls Fourier-space binning during loading.
            Larger values reduce data size and computation time, but also lower
            the final sampling resolution.

        Motion Fit Resolution
            Defines the spatial frequency range used during motion fitting.
            A wide range improves robustness, while high-frequency fitting may
            increase sensitivity to noise.

        B-Factor
            Downweights high spatial frequencies during motion estimation.

        Motion Model Grid
            Defines the spatial and temporal complexity of the motion model.
            Higher values allow more flexible correction but increase runtime.

        Even/Odd Averages
            Optionally exports independent averages from odd and even frames.
            This can be useful for denoiser training or validation procedures.

        Gain and Detector Defects

        Gain Transpose / Flip
            Allows correction of gain-reference orientation mismatches.

        EER Options
            For EER movies, virtual frame fractionation can be defined either by
            grouping frames or by specifying exposure per group.

        CTF Estimation

        Estimate CTF
            Enables per-tilt CTF fitting after motion correction.

        Window Size
            Defines the patch size used during CTF estimation.

        Resolution Range
            Sets the frequency interval used for fitting the CTF model.

        Defocus Search Range
            Defines the explored underfocus interval in microns.

        Defocus Model Grid
            Allows spatial or temporal modeling of defocus variation.

        Fit Phase
            Enables phase-shift estimation for phase plate data.

        Use Movie Average
            Uses the average movie spectrum instead of averaging individual-frame
            spectra. This can improve stability in low-signal datasets.

        Handedness Check
            Optionally evaluates defocus handedness consistency across the dataset.

        Processing Logic

        Initial Step
            The protocol initializes execution, stores sampling information, and
            creates the global frame-series Warp settings file.

        Per Tilt-Series Processing
            For each tilt-series:

            - movie file paths are collected
            - motion correction is executed
            - tilt metadata are prepared
            - tomostar metadata are generated
            - optional CTF estimation is performed
            - aligned outputs are written

        Output Generation

        Aligned Tilt-Series
            The protocol creates a new aligned tilt-series where:

            - tilt images are sorted by tilt angle
            - aligned averages are stacked into a single MRC stack
            - sampling rate is updated according to binning

        Even/Odd Stacks
            If enabled, separate odd and even aligned stacks are also generated.

        CTF Output
            When CTF estimation is enabled:

            - Warp XML output is parsed
            - defocus values are extracted
            - one CTF object per tilt image is created
            - PSD references are assigned

        Defocus Handedness

        If handedness evaluation is enabled:

            - Warp checks the global defocus handedness
            - the average correlation is parsed from stdout
            - a boolean output indicates whether handedness is consistent

        Validation

        Before execution, the protocol verifies that the selected CTF fitting
        resolution is not beyond the Nyquist frequency imposed by the selected
        binning factor.

        If the requested fitting resolution exceeds Nyquist limits, a validation
        warning is returned.

        Summary Output

        The protocol reports:

            - number of aligned tilt-series generated
            - number of CTF series estimated
            - handedness evaluation result (if requested)

        Practical Recommendations

        For routine cryo-ET processing:

            - start with binning = 1 unless data size is limiting
            - keep default motion fitting parameters for most datasets
            - enable CTF estimation when downstream reconstruction requires it
            - use even/odd averages only when explicitly needed
            - check handedness only for full dataset validation

        Final Perspective

        This protocol acts as a bridge between Scipion and WarpTools for
        high-throughput cryo-electron tomography preprocessing.

        Its main strength lies in combining movie motion correction, tilt-series
        assembly, CTF estimation, and metadata generation into a single
        reproducible workflow suitable for tomographic pipelines.
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
        form.addParallelSection(threads=2, mpi=0)
        form.addHidden(params.GPU_LIST, params.StringParam, default='0',
                       expertLevel=params.LEVEL_ADVANCED,
                       label="Choose GPU IDs",
                       help="Space-separated list of GPU IDs to use for processing. Default: all GPUs in the system."
                            " Warp can use multiple GPUs - in that case"
                            " set to i.e. *0 1 2*.")

    def _defineAlignmentParams(self, form):
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

    # --------------------------- STEPS functions -----------------------------

    def insertInitialSteps(self):
        self.numberOfThreads = Integer(2)
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

        self.runJob(Plugin.getProgram(WARP_TOOLS, CREATE_SETTINGS), cmd, executable='/bin/bash')

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
            argsDict['--exposure'] = exposure

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

        self.runJob(Plugin.getProgram(WARP_TOOLS, CREATE_SETTINGS), cmd, executable='/bin/bash')

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
        self.runJob(self.getPlugin().getProgram(WARP_TOOLS, TS_DEFOCUS_HAND), cmd, executable='/bin/bash')
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

        self.runJob(self.getPlugin().getProgram(WARP_TOOLS, FS_MOTION_AND_CTF), cmd, executable='/bin/bash')

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
        self.runJob(self.getPlugin().getProgram(WARP_TOOLS, TS_CTF), cmd, executable='/bin/bash')

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
    
    def _validate(self):
        errorMsg = []
        movieSampling = self.inputTSMovies.get().getSamplingRate()
        outputSampling = movieSampling*self.binFactor.get()
        nyquistFreq = 2 * outputSampling
        resToFit = self.range_max.get()

        if resToFit < nyquistFreq:
            warnMs = (f'The resolution to fit should be greater than nyquist. Currently, resolution to fit is {resToFit} '
                      f'and at this binning factor Nyquist is {nyquistFreq}\n')
            errorMsg.append(warnMs)
        return errorMsg

    def allowsDelete(self, obj):
        return True

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
