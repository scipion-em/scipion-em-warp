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
import logging
import math
import traceback
from enum import Enum
from os.path import splitext, abspath, dirname, join, basename
from typing import Union, Tuple, List
from pwem.emlib.image.image_readers import ImageStack, ImageReadersRegistry
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.object import Set, Float, Pointer
from pyworkflow.protocol import GPU_LIST, PointerParam, StringParam, LEVEL_ADVANCED, EnumParam, \
    FloatParam, IntParam, BooleanParam, LabelParam
from pyworkflow.utils import cyanStr, createLink, Message, makePath, redStr, replaceBaseExt, yellowStr
from tomo.objects import (SetOfTiltSeriesM, SetOfTiltSeries, TiltImage,
                          TiltSeries, SetOfCTFTomoSeries,
                          TiltSeriesM, TiltImageM, CTFTomoSeries, CTFTomo)
from warp import Plugin
from warp.constants import *
from warp.convert import writeTsStar
from warp.utils import tomoStarGenerate, parseCtfXMLFile

logger = logging.getLogger(__name__)

# Binning modes
BINNING_FACTOR = 0
TARGET_SAMPLING_RATE = 1

# EER grouping modes
EER_NGROUPS = 0
EER_GROUP_EXPOSURE = 1

# Gain flip options
NO_FLIP = 0
FLIP_X = 1
FLIP_Y = 2

# Gain transpose options
NO_SWAP = 0
TRANSPOSE_XY = 1


class WarpTsMcorrOutputs(Enum):
    tiltSeries = SetOfTiltSeries
    ctfs = SetOfCTFTomoSeries


class ProtWarpTSMotionCorr(EMProtocol):  # , ProtTSMovieAlignBase):
    """ This protocol wraps WarpTools programs.
        Estimate motion in frame series, produce aligned averages, estimate CTF
    """

    _label = 'tilt-series motion correction and ctf estimation'
    _devStatus = BETA
    _possibleOutputs = WarpTsMcorrOutputs
    evenOddCapable = True

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.samplingRate = None
        self.outSamplingRate = None
        self.exposure = None
        self.tsMDict = None
        self.tsmDir = None
        self.tsDir = None
        self.tsStarDir = None
        self.gainFile = None
        self.darkFile = None
        self.averageCorrelation = Float()
        self.failedTsIds = []

    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputTSMovies', PointerParam, pointerClass=SetOfTiltSeriesM,
                      important=True,
                      label=Message.LABEL_INPUT_MOVS,
                      help='Select a set of previously imported tilt series movies.')
        form.addSection('Alignment')
        self._defineAlignmentParams(form)
        # ProtTSMovieAlignBase._defineStreamingParams(self, form)
        form.addParallelSection(threads=2, mpi=0)
        form.addHidden(GPU_LIST, StringParam, default='0',
                       expertLevel=LEVEL_ADVANCED,
                       label="Choose GPU IDs",
                       help="Space-separated list of GPU IDs to use for processing. Default: all GPUs in the system."
                            " Warp can use multiple GPUs - in that case"
                            " set to i.e. *0 1 2*.")

    @staticmethod
    def _defineAlignmentParams(form):
        form.addParam('binFactorMode', EnumParam,
                      display=EnumParam.DISPLAY_HLIST,
                      choices=['Binning factor', 'Target sampling rate (Angst/pix)'],
                      default=BINNING_FACTOR,
                      label='Binning mode')
        form.addParam('binFactor', FloatParam,
                      condition=f'binFactorMode == {BINNING_FACTOR}',
                      allowsNull=True,
                      default=1,
                      label="Binning factor",
                      help="It supports non-integer values.")
        form.addParam('binTarget', FloatParam,
                      condition=f'binFactorMode == {TARGET_SAMPLING_RATE}',
                      allowsNull=True,
                      label="Sampling rate target (Angst/pix) for binning",
                      help="Choose the binning exponent automatically to match "
                           "this target pixel size in angstroms.")

        line = form.addLine('Resolution to fit',
                            help='Resolution in Angstrom to consider in fit.')
        line.addParam('m_range_min', FloatParam, default=500,
                      label='Min', help='Minimum resolution in Angstrom to consider in fit')
        line.addParam('m_range_max', FloatParam, default=10,
                      label='Max', help='Maximun resolution in Angstrom to consider in fit')

        form.addParam('bfactor', FloatParam, default=-500,
                      label="B-factor",
                      help="Downweight higher spatial frequencies using a "
                           "B-factor, in Angstrom^2")

        line = form.addLine('Motion model grid',
                            help="Resolution of the motion model grid in X, Y, and temporal dimensions, "
                                 "separated by 'x': e.g. 5x5x40; empty = auto")
        line.addParam('x', IntParam, default=None,
                      allowsNull=True,
                      label='X')
        line.addParam('y', IntParam,
                      default=None,
                      allowsNull=True,
                      label='Y')
        line.addParam('z', IntParam, default=None,
                      allowsNull=True,
                      label='Temporal')

        form.addParam('average_halves', BooleanParam,
                      default=False,
                      label='Do even and odd ?',
                      help='Export aligned averages of odd and even frames separately, e.g. for denoiser training')

        form.addSection(label="Gain and defects")
        form.addParam('gainSwap', EnumParam,
                      choices=['no swap', 'transpose X/Y'],
                      label="Transpose gain reference:",
                      default=NO_SWAP,
                      display=EnumParam.DISPLAY_COMBO)

        form.addParam('gainFlip', EnumParam,
                      choices=['no flip', 'flip X', 'flip Y'],
                      label="Flip gain reference:",
                      default=NO_FLIP,
                      display=EnumParam.DISPLAY_COMBO)

        form.addSection("EER")
        form.addParam('EERtext', LabelParam,
                      label="These options are ignored for non-EER files.")
        form.addParam('eerGroupMode', EnumParam,
                      display=EnumParam.DISPLAY_HLIST,
                      choices=['EER fractionation', 'EER group exposure'],
                      default=EER_NGROUPS,
                      label='ERR grouping mode')
        form.addParam('eer_ngroups', IntParam,
                      default=40,
                      allowsNull=True,
                      condition=f'eerGroupMode == {EER_NGROUPS}',
                      label='EER fractionation',
                      help="Number of groups to combine raw EER frames into, i.e. number of 'virtual' "
                           "frames in resulting stack; use negative value to specify the number of "
                           "frames per virtual frame instead")
        form.addParam('eer_groupexposure', FloatParam,
                      default=0.5,
                      allowsNull=True,
                      condition=f'eerGroupMode == {EER_GROUP_EXPOSURE}',
                      label='EER group exposure',
                      help="As an alternative to EER fractionation, fractionate the frames so that a group will "
                           "have this exposure in e-/A^2; this overrides EER fractionation"
                           "\nFractionate such that each fraction "
                           "has about 0.5 to 1.25 e/A2.")
        form.addSection(label="CTF")

        form.addParam('estimateCTF', BooleanParam, default=True,
                      label='Estimate the CTF ?',
                      help='Estimate the CTF')

        form.addParam('window', IntParam, default=512,
                      condition='estimateCTF',
                      label='Windows', help='Patch size for CTF estimation in binned pixels')

        line = form.addLine('Resolution (Å)',
                            condition='estimateCTF',
                            help='Resolution in Angstrom to consider in fit.')

        line.addParam('range_min', FloatParam, default=30,
                      condition='estimateCTF',
                      label='Min', help='Lowest (worst) resolution in Angstrom to consider in fit')

        line.addParam('range_max', FloatParam, default=4,
                      condition='estimateCTF',
                      label="Max",
                      help="Highest (best) resolution in Angstrom to consider in fit")

        line = form.addLine('Defocus search range (um)',
                            condition='estimateCTF',
                            help='Defocus values in um to explore during fitting (positive = underfocus). '
                                 'The units are microns!!')
        line.addParam('defocus_min', FloatParam, default=0.5,
                      condition='estimateCTF',
                      label='Min', help='Minimum defocus value in um to explore during fitting (positive = underfocus)')
        line.addParam('defocus_max', FloatParam, default=5,
                      condition='estimateCTF',
                      label='Max', help='Maximum defocus value in um to explore during fitting (positive = underfocus)')

        line = form.addLine('Defocus model grid',
                            condition='estimateCTF',
                            help="Resolution of the defocus model grid in X, Y, and temporal dimensions, "
                                 "separated by x: e.g. 5x5x40; empty = auto; Z > 1 is purely experimental")

        line.addParam('c_x', IntParam, default=None,
                      condition='estimateCTF',
                      allowsNull=True, label='X')
        line.addParam('c_y', IntParam, default=None,
                      condition='estimateCTF',
                      allowsNull=True, label='Y')
        line.addParam('c_z', IntParam, default=None, allowsNull=True,
                      condition='estimateCTF',
                      label='Temporal')

        form.addParam('fit_phase', BooleanParam, default=False,
                      condition='estimateCTF',
                      label='Fit phase', help='Fit the phase shift of a phase plate')

        form.addParam('use_sum', BooleanParam, default=False,
                      condition='estimateCTF',
                      label='Use the movie average',
                      help='Use the movie average spectrum instead of the average of individual '
                           'frames spectra. Can help in the absence of an energy filter, or when signal is low')

        form.addParam('handedness', BooleanParam, default=False,
                      condition='estimateCTF',
                      expertLevel=LEVEL_ADVANCED,
                      label='Check the handedness ?',
                      help='Checking defocus handedness across a dataset ')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._initialize()
        beginCtfReqDeps = []
        closeSetStepDeps = []
        pId = self._insertFunctionStep(self.createTsMovSettingsStep,
                                       prerequisites=[],
                                       needsGPU=False)

        pId = self._insertFunctionStep(self.createTsSettingsStep,
                                       prerequisites=pId,
                                       needsGPU=False)
        for tsId, tsM in self.tsMDict.items():
            pId = self._insertFunctionStep(self.processTsMStep, tsId,
                                           prerequisites=pId,
                                           needsGPU=True)

            pId = self._insertFunctionStep(self.createOutputTsStep, tsId,
                                           prerequisites=pId,
                                           needsGPU=False)
            beginCtfReqDeps.append(pId)
            closeSetStepDeps.append(pId)

        if self.estimateCTF.get():
            for tsId, tsM in self.tsMDict.items():
                pId = self._insertFunctionStep(self.createTsMStarStep, tsId,
                                               prerequisites=beginCtfReqDeps,
                                               needsGPU=False)

                pId = self._insertFunctionStep(self.estimateCtfStep, tsId,
                                               prerequisites=pId,
                                               needsGPU=True)

                pId = self._insertFunctionStep(self.createOutputCtfStep, tsId,
                                               prerequisites=pId,
                                               needsGPU=False)

            closeSetStepDeps.append(pId)

        self._insertFunctionStep(self.closeOutputStep,
                                 prerequisites=closeSetStepDeps,
                                 needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        tsMSet = self.getInputTSMovies()
        self.samplingRate = tsMSet.getSamplingRate()
        self.exposure = tsMSet.getAcquisition().getDosePerFrame()
        self.outSamplingRate = tsMSet.getSamplingRate() * self.binFactor.get()
        self.tsMDict = {tsM.getTsId(): tsM.clone() for tsM in tsMSet.iterItems()}
        # Gain and defects
        gainFile = tsMSet.getGain()
        darkFile = tsMSet.getDark()
        self.gainFile = abspath(gainFile) if gainFile else None
        self.darkFile = abspath(darkFile) if darkFile else None
        # Make paths
        self.tsmDir = self._getFrameSeriesDir()
        self.tsDir = self._getTiltSeriesDir()
        self.tsStarDir = self._getTsStarDir()
        makePath(*[self.tsmDir, self.tsDir, self.tsStarDir])

    def createTsMovSettingsStep(self):
        logger.info(cyanStr(">>> Creating frame-series settings..."))
        try:
            cmd = self._genCreateSettingsArgs()
            self.runJob(Plugin.getProgram(WARP_TOOLS, CREATE_SETTINGS), cmd, executable='/bin/bash')
        except Exception as e:
            raise Exception(f"{WARP_TOOLS} {CREATE_SETTINGS} failed with the exception --> {e}")

    def createTsSettingsStep(self):
        logger.info(cyanStr(">>> Creating tilt-series settings..."))
        try:
            cmd = self._genCreateSettingsTsArgs()
            self.runJob(Plugin.getProgram(WARP_TOOLS, CREATE_SETTINGS), cmd, executable='/bin/bash')
        except Exception as e:
            raise Exception(f"{WARP_TOOLS} {CREATE_SETTINGS} failed with the exception --> {e}")

    def processTsMStep(self, tsId: str):
        try:
            logger.info(cyanStr(f">>> tsId = {tsId} - performing the motion-correction..."))
            tsMovie = self.tsMDict[tsId]
            warpMoviesNamesList = [abspath(tiName.getFileName()) for tiName in tsMovie.iterItems()]
            warpMoviesNamesList = " ".join(warpMoviesNamesList)
            args = self._genFsMotionAndCTFArgs(tsMovie, warpMoviesNamesList)
            self.runJob(Plugin.getProgram(WARP_TOOLS, FS_MOTION_AND_CTF), args, executable='/bin/bash')
        except Exception as e:
            logger.error(redStr(f"{WARP_TOOLS} {FS_MOTION_AND_CTF} failed with the exception --> {e}"))
            traceback.print_exc()
            self.failedTsIds.append(tsId)

    def createOutputTsStep(self, tsId: str):
        if tsId in self.failedTsIds:
            return

        try:
            logger.info(cyanStr(f">>> tsId = {tsId} - creating the output tilt-series..."))
            newTs, tiList = self._prepareOutputTs(tsId)
            self._registerOutputTs(newTs, tiList)

        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to register the output TS with exception {e}. Skipping... '))
            logger.error(traceback.format_exc())
            self.failedTsIds.append(tsId)

    def createTsMStarStep(self, tsId: str):
        if tsId in self.failedTsIds:
            return

        try:
            logger.info(cyanStr(f">>> {tsId} - Creating the star file..."))
            writeTsStar(self.tsMDict[tsId], self._getAvgDir(), self._getOutTsStarFile(tsId))

        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to create the TS star file '
                                f'with exception {e}. Skipping... '))
            logger.error(traceback.format_exc())
            self.failedTsIds.append(tsId)

    def estimateCtfStep(self, tsId: str):
        if tsId in self.failedTsIds:
            return

        try:
            logger.info(cyanStr(f">>> {tsId} - estimating the CTF..."))
            tsMovie = self.tsMDict[tsId]
            args = self._getCtfEstimationArgs(tsMovie)
            self.runJob(Plugin.getProgram(WARP_TOOLS, TS_CTF), args, executable='/bin/bash')

        except Exception as e:
            logger.error(redStr(f"{WARP_TOOLS} {TS_CTF} failed with the exception --> {e}"))
            logger.error(traceback.format_exc())
            self.failedTsIds.append(tsId)

    def createOutputCtfStep(self, tsId: str):
        if tsId in self.failedTsIds:
            return

        try:
            logger.info(cyanStr(f">>> tsId = {tsId} - creating the output CTF..."))
            newCtfTomoSeries, ctfList = self._prepareOutputCtf(tsId)
            self._registerOutputCtf(newCtfTomoSeries, ctfList)

        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to register the output CTF with exception {e}. Skipping... '))
            logger.error(traceback.format_exc())
            self.failedTsIds.append(tsId)

    def _prepareOutputCtf(self, tsId: str) -> Tuple[CTFTomoSeries, List[CTFTomo]]:
        tsDir = self._getTiltSeriesDir()
        psdStack = join(tsDir, POWERSPECTRUM_FOLDER, tsId + '.mrc')
        outTsSet = getattr(self, self._possibleOutputs.tiltSeries.name)
        ts = outTsSet.getItem(TiltSeries.TS_ID_FIELD, tsId)
        newCTFTomoSeries = CTFTomoSeries(tsId=tsId)
        newCTFTomoSeries.copyInfo(ts)
        newCTFTomoSeries.setTiltSeries(ts)

        defocusFilePath = join(tsDir, tsId + '.xml')
        ctfData, gridCtfData = parseCtfXMLFile(defocusFilePath)
        defaultDefocusDelta = float(ctfData['DefocusDelta']) * 1e4
        defaultDefocusAngle = float(ctfData['DefocusAngle'])

        ctfTomoList = []
        acqOrderList = [ti.getAcquisitionOrder() for ti in ts.iterItems(orderBy=TiltImage.TILT_ANGLE_FIELD)]
        for index, acqOrder in enumerate(acqOrderList):
            if index not in gridCtfData["Nodes"]:
                logger.info(yellowStr(f'\t>>> tsId = {tsId} - No Warp CTF defocus found for tilt index {index}'))
                continue

            defocus = gridCtfData["Nodes"][index]
            defocusDelta = gridCtfData["DeltaNodes"].get(index, defaultDefocusDelta)
            defocusAngle = gridCtfData["AngleNodes"].get(index, defaultDefocusAngle)
            defocusU = defocus + defocusDelta
            defocusV = defocus - defocusDelta

            itemIndex = index + 1
            newCTFTomo = CTFTomo()
            newCTFTomo.setAcquisitionOrder(acqOrder)
            newCTFTomo.setIndex(itemIndex)
            newCTFTomo.setDefocusU(defocusU)
            newCTFTomo.setDefocusV(defocusV)
            newCTFTomo.setDefocusAngle(defocusAngle)
            newCTFTomo.setResolution(0)
            newCTFTomo.setFitQuality(0)
            newCTFTomo.standardize()
            newCTFTomo.setPsdFile(f"{itemIndex}@{psdStack}")

            ctfTomoList.append(newCTFTomo)

    def _registerOutputCtf(self, newCTFTomoSeries: CTFTomoSeries, ctfTomoList: List[CTFTomo]):
        # CTFTomoSeries set
        outputCtfSet = self.getOutputSetOfCTFTomoSeries(OUTPUT_CTF_SERIES)
        # CTFTomoSeries
        outputCtfSet.append(newCTFTomoSeries)
        # Tilt-images
        for ctfTomo in ctfTomoList:
            newCTFTomoSeries.append(ctfTomo)
        # Data persistence
        newCTFTomoSeries.write()
        outputCtfSet.update(newCTFTomoSeries)
        outputCtfSet.write()
        self._store(outputCtfSet)

    def closeOutputStep(self):
        super()._closeOutputSet()
        outputTsSet = getattr(self, self._possibleOutputs.tiltSeries.name, ())
        if not outputTsSet or outputTsSet and len(outputTsSet) == 0:
            raise Exception('No outputs were generated. Please check the logs run.stdout and run.stderr.')

    # --------------------------- UTILS functions -----------------------------
    def getInputTSMovies(self, asPointer: bool = False) -> Union[SetOfTiltSeriesM, Pointer]:
        inTsMSetPointer = self.inputTSMovies
        return inTsMSetPointer if asPointer else inTsMSetPointer.get()

    def _getFrameSeriesSettingsFn(self) -> str:
        return abspath(self._getExtraPath(FRAMESERIES_SETTINGS))

    def _getFrameSeriesDir(self) -> str:
        return abspath(self._getExtraPath(FRAMESERIES_FOLDER))

    def _getAvgDir(self) -> str:
        return join(self._getExtraPath(FRAMESERIES_FOLDER), AVERAGE_FOLDER)

    def _getTsStarDir(self) -> str:
        return abspath(self._getExtraPath(TOMOSTAR_FOLDER))

    def _getOutTsStarFile(self, tsId: str) -> str:
        return join(self.tsStarDir, f'{tsId}{TOMOSTAR_EXT}')

    def _getTiltSeriesSettingsFn(self) -> str:
        return abspath(self._getExtraPath(TILTSERIES_SETTINGS))

    def _getTiltSeriesDir(self) -> str:
        return abspath(self._getExtraPath(TILTSERIES_FOLDER))

    def _genCreateSettingsArgs(self) -> str:
        ts_movies = self.getInputTSMovies()
        first_ts_movie = ts_movies.getFirstItem()
        file_name, extension = splitext(first_ts_movie.getFirstItem().getFileName())
        folder_data = abspath(dirname(file_name))

        # Initialize argument list with base parameters
        args = [
            f"--folder_data {folder_data}",
            f"--extension '*{extension}'",
            f"--folder_processing {self.tsmDir}",
            f"--bin {self._getBinFactor()}",
            f"--angpix {self.samplingRate}",
            f"--exposure {self.exposure}",
            f"--output {self._getFrameSeriesSettingsFn()}",

        ]
        if self.fit_phase.get():
            args.append("--fit_phase")

        # Handle EER specific arguments
        if extension == ".eer":
            eer_grouping_mode = self.eerGroupMode.get()
            eer_groups = self.eer_ngroups.get()
            eer_group_exp = self.eer_groupexposure.get()

            if eer_grouping_mode == EER_NGROUPS and eer_groups is not None:
                args.append(f"--eer_ngroups {eer_groups}")
            elif eer_grouping_mode == EER_GROUP_EXPOSURE and eer_group_exp is not None:
                args.append(f"--eer_groupexposure {eer_group_exp}")

        # Handle gain file and transformations
        self._genGainAndDarkCmd(args)

        return " ".join(args)

    def _genCreateSettingsTsArgs(self) -> str:
        x, y, z = self._getTomoDims()
        args = [
            f'--output {self._getTiltSeriesSettingsFn()}',
            f'--folder_processing {self.tsDir}',
            f'--folder_data {self.tsStarDir}',
            f'--extension *{TOMOSTAR_EXT}',
            f'--angpix {self.samplingRate}',
            f'--exposure {self.exposure}',
            f'--tomo_dimensions {x}x{y}x{z}'
        ]
        self._genGainAndDarkCmd(args)
        return ' '.join(args)

    def _genFsMotionAndCTFArgs(self, tsMovie: TiltSeriesM, warpMoviesNames: str) -> str:
        """Execute frame series motion correction and CTF estimation using Warp tools.

        Args:
            tsMovie (TiltSeriesM): The tilt series movie object containing acquisition parameters.
            warpMoviesNames (str): Space-separated list or pattern of target movie absolute paths.
        """
        acquisition = tsMovie.getAcquisition()

        # Construct base argument list directly
        args = [
            f"--settings {self._getFrameSeriesSettingsFn()}",
            f"--m_range_min {self.m_range_min.get()}",
            f"--m_range_max {self.m_range_max.get()}",
            f"--m_bfac {self.bfactor.get()}",
            f"--c_window {self.window.get()}",
            f"--c_range_min {self.range_min.get()}",
            f"--c_range_max {self.range_max.get()}",
            f"--c_defocus_min {self.defocus_min.get()}",
            f"--c_defocus_max {self.defocus_max.get()}",
            f"--c_voltage {int(acquisition.getVoltage())}",
            f"--c_cs {acquisition.getSphericalAberration()}",
            f"--c_amplitude {acquisition.getAmplitudeContrast()}",
            f"--input_data {warpMoviesNames}",
            f"--output_processing {self.tsmDir}",
        ]

        # GPU and threads configuration
        self._addThreadsAndGpusArgs(args)

        # Default output parameters
        args.append("--out_averages")

        if self.average_halves.get():
            args.append("--out_average_halves")

        # Motion grid dimensions
        grid_x, grid_y, grid_z = self.x.get(), self.y.get(), self.z.get()
        if grid_x and grid_y and grid_z:
            args.append(f"--m_grid {grid_x}x{grid_y}x{grid_z}")

        # CTF grid dimensions
        c_grid_x, c_grid_y, c_grid_z = self.c_x.get(), self.c_y.get(), self.c_z.get()
        if c_grid_x and c_grid_y and c_grid_z:
            args.append(f"--c_grid {c_grid_x}x{c_grid_y}x{c_grid_z}")

        # Optional processing flags
        if self.fit_phase.get():
            args.append("--c_fit_phase")

        if self.use_sum.get():
            args.append("--c_use_sum")

        return " ".join(args)

    def _getCtfEstimationArgs(self, tsMovie: TiltSeriesM) -> str:
        acq = tsMovie.getAcquisition()
        args = [
            f'--settings {self._getTiltSeriesSettingsFn()}',
            f'--window {self.window.get()}',
            f'--range_low {self.range_min.get()}',
            f'--range_high {self.range_max.get()}',
            f'--defocus_min {self.defocus_min.get()}',
            f'--defocus_max {self.defocus_max.get()}',
            f'--voltage {acq.getVoltage():.0f}',
            f'--cs {acq.getSphericalAberration()}',
            f'--amplitude {acq.getAmplitudeContrast()}',
            f'--auto_hand {self._getAutoHandNTs()}'
        ]

        # GPU and threads configuration
        self._addThreadsAndGpusArgs(args)

        return " ".join(args)

    def _prepareOutputTs(self, tsId: str) -> Tuple[TiltSeries, List[TiltImage]]:
        tsMovie = self.tsMDict[tsId]
        averageFolder = self._getAvgDir()
        properties = {"sr": self.outSamplingRate}
        newStack = ImageStack(properties=properties)
        oddFileNames = ImageStack(properties=properties)
        evenFileNames = ImageStack(properties=properties)
        newBinaryName = join(averageFolder, f'{tsId}{MRCS_EXT}')
        hasAverageHalves = self.average_halves.get()
        newOddBinaryName = join(averageFolder, ODD, f'{tsId}_{ODD}{MRCS_EXT}')
        newEvenBinaryName = join(averageFolder, EVEN, f'{tsId}_{EVEN}{MRCS_EXT}')

        tiList = []
        for i, tiM in enumerate(tsMovie.iterItems(orderBy=TiltImageM.TILT_ANGLE_FIELD)):
            newTi = TiltImage()
            newTi.copyInfo(tiM)
            newTi.setFileName(newBinaryName)
            newTi.setIndex(i + 1)
            newTi.setSamplingRate(self.outSamplingRate)

            # Mount the stacks
            averageFn = join(averageFolder, replaceBaseExt(tiM.getFileName(), 'mrc'))
            newStack.append(ImageReadersRegistry.open(averageFn))
            if hasAverageHalves:
                newTi.setOddEven([newOddBinaryName, newEvenBinaryName])
                averageBaseName = basename(averageFn)
                oddFileNames.append(ImageReadersRegistry.open(join(averageFolder, ODD, averageBaseName)))
                evenFileNames.append(ImageReadersRegistry.open(join(averageFolder, EVEN, averageBaseName)))

            tiList.append(newTi)

        # Write the stacks
        ImageReadersRegistry.write(newStack, newBinaryName,
                                   isStack=True,
                                   samplingRate=self.outSamplingRate)
        if hasAverageHalves:
            ImageReadersRegistry.write(oddFileNames, newOddBinaryName,
                                       isStack=True,
                                       samplingRate=self.outSamplingRate)
            ImageReadersRegistry.write(evenFileNames, newEvenBinaryName,
                                       isStack=True,
                                       samplingRate=self.outSamplingRate)

        newTs = TiltSeries(tsId=tsId)
        newTs.copyInfo(tsMovie)

        return newTs, tiList

    def _registerOutputTs(self, newTs: TiltSeries, tiList: List[TiltImage]):
        # TS set
        outTsSet = self.getOutputSetOfTS()
        # TS
        outTsSet.append(newTs)
        # Tilt-images
        for newTi in tiList:
            newTs.append(newTi)
        # Data persistence
        newTs.write()
        outTsSet.update(newTs)
        outTsSet.write()
        self._store(outTsSet)

    # def createTiltSeriesSettingStep(self, tsId):
    #     self.info(">>> Starting tilt-series settings creation (%s)..." % tsId)
    #     setOfTSMovies = self.inputTSMovies.get()
    #     sr = setOfTSMovies.getSamplingRate()
    #     exposure = setOfTSMovies.getAcquisition().getDosePerFrame()
    #     firstTSMovie = setOfTSMovies.getFirstItem()
    #     fileName, extension = splitext(firstTSMovie.getFirstItem().getFileName())
    #     settingsFolder = abspath(self._getExtraPath(SETTINGS_FOLDER))
    #     makePath(settingsFolder)
    #     processingFolder = abspath(self._getExtraPath(TILTSERIES_FOLDER))
    #     makePath(processingFolder)
    #     tsSettingFile = tsId + '_' + TILTSERIES_SETTINGS
    #     tsSettingFilePath = abspath(join(self._getExtraPath(settingsFolder), tsSettingFile))
    #     argsDict = {
    #         "--folder_data": abspath(self._getExtraPath(TOMOSTAR_FOLDER)),
    #         "--extension": "%s.tomostar" % tsId,
    #         "--folder_processing": processingFolder,
    #         '--angpix': sr,
    #         "--output": tsSettingFilePath
    #     }
    #
    #     if self.binfactorMode.get() == BINNING_FACTOR:
    #         argsDict["--bin"] = self._getBinFactor(),
    #     else:
    #         argsDict["--bin_angpix"] = self.binTarget.get()
    #
    #     if exposure is not None:
    #         argsDict['--exposure'] = exposure
    #
    #     if hasattr(self, 'tomo_thickness'):
    #         z = self.tomo_thickness.get()
    #         x = self.x_dimension.get() or setOfTSMovies.getDimensions()[0]
    #         y = self.y_dimension.get() or setOfTSMovies.getDimensions()[1]
    #
    #         argsDict['--tomo_dimensions'] = f'{x}x{y}x{z}'
    #
    #     if extension == '.eer':
    #         argsDict['--eer_ngroups'] = self.eer_ngroups.get()
    #         if self.eer_groupexposure.get():
    #             argsDict['--eer_groupexposure'] = self.eer_groupexposure.get()
    #
    #     cmd = ' '.join(['%s %s' % (k, v) for k, v in argsDict.items()])
    #
    #     self.runJob(Plugin.getProgram(WARP_TOOLS, CREATE_SETTINGS), cmd, executable='/bin/bash')
    #
    # # def dataPrepare(self, tsMovie):
    #     """Creates the setting file that will be used by the different programs.
    #        It also extracts the tiltimages from the tiltseries and generates the *.tomostar files based on
    #        the tiltimages."""
    #     starFolder = self._getExtraPath(TOMOSTAR_FOLDER)
    #     makePath(starFolder)
    #     imagesFolder = self._getExtraPath(FRAMES_FOLDER)
    #     makePath(imagesFolder)
    #
    #     if tsMovie.isEnabled():
    #         tsId = tsMovie.getTsId()
    #         tiValues = {}
    #         for ti in tsMovie.iterItems():
    #             if ti.isEnabled():  # Excluding views
    #                 dose = 0
    #                 maskedFraction = 0
    #                 shiftX = 0
    #                 shiftY = 0
    #                 axisAngle = 0
    #                 amplitudeContrast = 0
    #
    #                 if tsMovie.hasAcquisition():
    #                     axisAngle = tsMovie.getAcquisition().getTiltAxisAngle()
    #                 if ti.getAcquisition():
    #                     amplitudeContrast = ti.getAcquisition().getAmplitudeContrast()
    #                     dose = ti.getAcquisition().getAccumDose()
    #                 fileName = ti.getFileName()
    #                 newBinaryName = basename(fileName)
    #                 createLink(abspath(fileName), join(imagesFolder, basename(fileName)))
    #
    #                 tiltAngle = ti.getTiltAngle()
    #
    #                 tiValues[tiltAngle] = [
    #                     newBinaryName,
    #                     -tiltAngle,
    #                     axisAngle,
    #                     shiftX,
    #                     shiftY,
    #                     dose,
    #                     amplitudeContrast,
    #                     maskedFraction
    #                 ]
    #
    #         tomoStarGenerate(tsId, tiValues, starFolder, 0)
    #
    # def createHandednessSetting(self):
    #     """Create a Warp settings file containing all tilt-series."""
    #     tsMovies = self.inputTSMovies.get()
    #     sr = tsMovies.getSamplingRate()
    #     exposure = tsMovies.getAcquisition().getDosePerFrame()
    #
    #     settingsPath = abspath(
    #         self._getExtraPath('defocus_hand.settings')
    #     )
    #     processingFolder = abspath(
    #         self._getExtraPath(TILTSERIES_FOLDER)
    #     )
    #
    #     argsDict = {
    #         "--folder_data": abspath(
    #             self._getExtraPath(TOMOSTAR_FOLDER)
    #         ),
    #         "--extension": "'*.tomostar'",
    #         "--folder_processing": processingFolder,
    #         "--bin": self.getBinFactor(),
    #         "--angpix": sr,
    #         "--output": settingsPath
    #     }
    #
    #     if exposure is not None:
    #         argsDict['--exposure'] = exposure
    #
    #     if hasattr(self, 'tomo_thickness'):
    #         z = self.tomo_thickness.get()
    #         x = self.x_dimension.get() or tsMovies.getDimensions()[0]
    #         y = self.y_dimension.get() or tsMovies.getDimensions()[1]
    #         argsDict['--tomo_dimensions'] = f'{x}x{y}x{z}'
    #
    #     cmd = ' '.join(
    #         ['%s %s' % (k, v) for k, v in argsDict.items()]
    #     )
    #
    #     self.runJob(
    #         Plugin.getProgram(WARP_TOOLS, CREATE_SETTINGS),
    #         cmd,
    #         executable='/bin/bash'
    #     )
    #
    #     return settingsPath
    #
    # def tsDefocusHandStep(self):
    #     """Check defocus handedness across the complete dataset."""
    #     self.info(">>> Starting defocus handedness...")
    #
    #     settingsPath = self.createHandednessSetting()
    #
    #     argsDict = {
    #         "--settings": settingsPath,
    #     }
    #
    #     cmd = ' '.join(
    #         ['%s %s' % (k, v) for k, v in argsDict.items()]
    #     )
    #     cmd += ' --check'
    #
    #     self.runJob(
    #         self.getPlugin().getProgram(WARP_TOOLS, TS_DEFOCUS_HAND),
    #         cmd,
    #         executable='/bin/bash'
    #     )
    #
    #     self.createOutputDefocusHand()
    #
    # def proccessTSMoviesStep(self, tsId) -> None:
    #     """Estimate motion in frame series, produce aligned averages and register the output"""
    #     tsMovie = self.getInputTSMovies().getItem(TiltSeries.TS_ID_FIELD, tsId)
    #     warpMoviesNamesList = [abspath(tiName.getFileName()) for tiName in tsMovie.iterItems()]
    #     warpMoviesNamesList = " ".join(warpMoviesNamesList)
    #     self.info(">>> Starting estimate motion for %s..." % tsId)
    #     self.fsMotionAndCTF(tsMovie, warpMoviesNamesList)
    #     if self.estimateCTF.get():
    #         self.createTiltSeriesSettingStep(tsId)
    #         self.dataPrepare(tsMovie)
    #         self.tsCtfEstimationStep(tsId)
    #     with self._lock:
    #         self.createOutputTS(tsMovie)
    #         if self.estimateCTF.get():
    #             self.createOutputCTF(tsId)
    #
    # def insertFinalSteps(self, proccessTSMoviesSteps) -> list:
    #     """The final steps inserted into the protocol"""
    #     finalSteps = []
    #     if self.handedness.get():
    #         finalStep = self._insertFunctionStep(self.tsDefocusHandStep, prerequisites=proccessTSMoviesSteps,
    #                                              needsGPU=True)
    #         finalSteps.append(finalStep)
    #     return finalSteps
    #
    # def tsCtfEstimationStep(self, tsId):
    #     """CTF estimation"""
    #     self.info(">>> Starting ctf estimation to %s" % tsId)
    #     inputTSAdquisition = self.inputTSMovies.get().getFirstItem().getAcquisition()
    #     settingFile = self._getExtraPath(SETTINGS_FOLDER, tsId + '_' + TILTSERIE_SETTINGS)
    #     argsDict = {
    #         "--settings": abspath(settingFile),
    #         "--window": self.window.get(),
    #         "--range_low": self.range_min.get(),
    #         "--range_high": self.range_max.get(),
    #         "--defocus_min": self.defocus_min.get(),
    #         "--defocus_max": self.defocus_max.get(),
    #         "--voltage": int(inputTSAdquisition.getVoltage()),
    #         "--cs": inputTSAdquisition.getSphericalAberration(),
    #         "--amplitude": inputTSAdquisition.getAmplitudeContrast(),
    #     }
    #
    #     gpuList = self.getGpuList()
    #     if gpuList:
    #         argsDict['--device_list'] = ' '.join(map(str, gpuList))
    #
    #     cmd = ' '.join(['%s %s' % (k, v) for k, v in argsDict.items()])
    #     self.runJob(self.getPlugin().getProgram(WARP_TOOLS, TS_CTF), cmd, executable='/bin/bash')
    #
    # def createOutputTS(self, tsMovie):
    #     self.info(">>> Generating output for %s..." % tsMovie.getTsId())
    #     outputTS = self.getOutputSetOfTS(OUTPUT_TILTSERIES)
    #     averageFolder = join(self._getExtraPath(FRAMESERIES_FOLDER), AVERAGE_FOLDER)
    #
    #     tsId = tsMovie.getTsId()
    #     newTs = TiltSeries(tsId=tsId)
    #     outputTS.append(newTs)
    #
    #     tiOrderDict = {}
    #     properties = {"sr": tsMovie.getSamplingRate()}
    #     newStack = ImageStack(properties=properties)
    #     oddFileNames = ImageStack(properties=properties)
    #     evenFileNames = ImageStack(properties=properties)
    #     newBinaryName = join(averageFolder, tsId + '.mrcs')
    #     hasAverageHalves = self.average_halves.get()
    #     newOddBinaryName = join(averageFolder, 'odd', tsId + '_odd.mrcs')
    #     newEvenBinaryName = join(averageFolder, 'even', tsId + '_even.mrcs')
    #
    #     for index, tiM in enumerate(tsMovie):
    #         fileName = splitext(basename(tiM.getFileName()))[0] + '.mrc'
    #         newTi = TiltImage(location=(index + 1, newBinaryName))
    #         newTi.copyInfo(tiM)
    #         newTi.setAcquisition(tiM.getAcquisition().clone())
    #         newTi.setSamplingRate(tiM.getSamplingRate() * self.binFactor.get())
    #         tiOrderDict[newTi.getTiltAngle()] = (newTi, fileName)
    #         if hasAverageHalves:
    #             newTi.setOddEven([newOddBinaryName, newEvenBinaryName])
    #
    #     sortedTiltAngle = sorted(tiOrderDict.keys())
    #
    #     for angle in sortedTiltAngle:
    #         fileName = tiOrderDict[angle][1]
    #         newStack.append(ImageReadersRegistry.open(join(averageFolder, fileName)))
    #         if hasAverageHalves:
    #             oddFileNames.append(ImageReadersRegistry.open(join(averageFolder, 'odd', fileName)))
    #             evenFileNames.append(ImageReadersRegistry.open(join(averageFolder, 'even', fileName)))
    #
    #     ImageReadersRegistry.write(newStack, newBinaryName, isStack=True)
    #     if hasAverageHalves:
    #         ImageReadersRegistry.write(oddFileNames, newOddBinaryName, isStack=True)
    #         ImageReadersRegistry.write(evenFileNames, newEvenBinaryName, isStack=True)
    #
    #     for index, angle in enumerate(sortedTiltAngle):
    #         ti = tiOrderDict[angle][0]
    #         ti.setIndex(index + 1)
    #         ti.setObjId(index + 1)
    #         newTs.append(ti)
    #
    #     outputTS.update(newTs)
    #     outputTS.write()
    #     self._store(outputTS)
    #
    # def deleteIntermediateOutputsStep(self):
    #     try:
    #         averageFolder = join(self._getExtraPath(FRAMESERIES_FOLDER), AVERAGE_FOLDER)
    #         if not exists(averageFolder):
    #             logger.info(f"The directory {averageFolder} does not exist.")
    #             return
    #         for filename in os.listdir(averageFolder):
    #             if filename.endswith(".mrc") or filename.endswith(".json"):
    #                 file_path = join(averageFolder, filename)
    #                 os.remove(file_path)
    #
    #         logger.info("All .mrc files have been deleted.")
    #
    #     except Exception as e:
    #         logger.error(f"An error occurred: {e}")
    #
    # def createOutputCTF(self, tsId):
    #     self.info(">>> Generating outputs to %s" % tsId)
    #     processingFolder = abspath(self._getExtraPath(TILTSERIES_FOLDER))
    #     tsSet = self.TiltSeries
    #
    #     if tsSet:
    #         psdStack = join(processingFolder, POWERSPECTRUM_FOLDER, tsId + '.mrc')
    #         ts = tsSet.getItem(TiltSeries.TS_ID_FIELD, tsId)
    #
    #         if ts.isEnabled():
    #             tsId = ts.getTsId()
    #             outputSetOfCTFTomoSeries = self.getOutputSetOfCTFTomoSeries(OUTPUT_CTF_SERIE)
    #
    #             newCTFTomoSeries = CTFTomoSeries(tsId=tsId)
    #             newCTFTomoSeries.copyInfo(ts)
    #             newCTFTomoSeries.setTiltSeries(ts)
    #             outputSetOfCTFTomoSeries.append(newCTFTomoSeries)
    #
    #             defocusFilePath = join(processingFolder, tsId + '.xml')
    #             ctfData, gridCtfData = parseCtfXMLFile(defocusFilePath)
    #
    #             defaultDefocusDelta = float(ctfData['DefocusDelta']) * 1e4
    #             defaultDefocusAngle = float(ctfData['DefocusAngle'])
    #
    #             tiltImages = sorted(
    #                 (ti for ti in ts.iterItems(iterate=False) if ti.isEnabled()),
    #                 key=lambda ti: ti.getTiltAngle()
    #             )
    #
    #             for index, ti in enumerate(tiltImages):
    #                 if index not in gridCtfData["Nodes"]:
    #                     raise ValueError(
    #                         f'No Warp CTF defocus found for tilt index {index} '
    #                         f'in tilt-series {tsId}'
    #                     )
    #
    #                 defocus = gridCtfData["Nodes"][index]
    #                 defocusDelta = gridCtfData["DeltaNodes"].get(
    #                     index, defaultDefocusDelta
    #                 )
    #                 defocusAngle = gridCtfData["AngleNodes"].get(
    #                     index, defaultDefocusAngle
    #                 )
    #
    #                 defocusU = defocus + defocusDelta
    #                 defocusV = defocus - defocusDelta
    #
    #                 itemIndex = index + 1
    #
    #                 newCTFTomo = CTFTomo()
    #                 newCTFTomo.setAcquisitionOrder(ti.getAcquisitionOrder())
    #                 newCTFTomo.setIndex(itemIndex)
    #                 newCTFTomo.setDefocusU(defocusU)
    #                 newCTFTomo.setDefocusV(defocusV)
    #                 newCTFTomo.setDefocusAngle(defocusAngle)
    #                 newCTFTomo.setResolution(0)
    #                 newCTFTomo.setFitQuality(0)
    #                 newCTFTomo.standardize()
    #                 newCTFTomo.setPsdFile(f"{itemIndex}@{psdStack}")
    #
    #                 newCTFTomoSeries.append(newCTFTomo)
    #
    #             outputSetOfCTFTomoSeries.update(newCTFTomoSeries)
    #             outputSetOfCTFTomoSeries.write()
    #             self._store(outputSetOfCTFTomoSeries)
    #
    # def createOutputDefocusHand(self):
    #     stdoutFile = abspath(
    #         join(self.getPath(), 'logs', 'run.stdout')
    #     )
    #
    #     correlation = None
    #
    #     with open(stdoutFile, 'r', encoding='utf-8') as file:
    #         for line in reversed(file.readlines()):
    #             if 'Average correlation:' in line:
    #                 correlation = float(line.split()[-1])
    #                 break
    #
    #     if correlation is None:
    #         raise RuntimeError('Warp did not report an average defocus-hand correlation.')
    #
    #     self.averageCorrelation.set(correlation)
    #
    #     self._defineOutputs(
    #         **{OUTPUT_HANDEDNESS: Boolean(correlation > 0)}
    #     )
    #
    #     self._store(self.averageCorrelation)

    def _summary(self):
        summary = []
        tilseriesSize = 0
        ctfSize = 0
        if self.hasAttribute(OUTPUT_TILTSERIES):
            tilseriesSize = self.TiltSeries.getSize()
        else:
            self.averageCorrelation = Float()
        summary.append(f"Aligned tiltseries: {tilseriesSize} of {self.inputTSMovies.get().getSize()}")

        if self.hasAttribute(OUTPUT_CTF_SERIES):
            ctfSize = self.CTFTomoSeries.getSize()
        summary.append(f"CTF estimated: {ctfSize} of {self.inputTSMovies.get().getSize()}")

        if self.handedness.get():
            correlation = self.averageCorrelation.get()

            if correlation > 0:
                summary.append(f'Defocus handedness: {correlation:.3f} '
                               '(no flip required)')
            elif correlation < 0:
                summary.append(f'Defocus handedness: {correlation:.3f} '
                               '(flip required)'
                               )
            else:
                summary.append('Defocus handedness: Not ready')

        return summary

    def _validate(self):
        errorMsg = []
        movieSampling = self.inputTSMovies.get().getSamplingRate()
        outputSampling = movieSampling * self.binFactor.get()
        nyquistFreq = 2 * outputSampling
        resToFit = self.range_max.get()

        if resToFit < nyquistFreq:
            warnMs = (
                f'The resolution to fit should be greater than nyquist. Currently, resolution to fit is {resToFit} '
                f'and at this binning factor Nyquist is {nyquistFreq}\n')
            errorMsg.append(warnMs)

        binMode = self.binFactorMode.get()
        binFactor = self.binFactor.get()
        binTarget = self.binTarget.get()
        if binMode == BINNING_FACTOR:
            if binFactor is None:
                errorMsg.append('Binning factor cannot be empty.')
            elif binFactor < 1:
                errorMsg.append('Binning factor must be greater or equal than 1.')
        if binMode == TARGET_SAMPLING_RATE:
            if binTarget is None:
                errorMsg.append('Target sampling rate for binning cannot be empty.')
            elif binTarget <= 0:
                errorMsg.append('Target sampling rate for binning must be greater than 0.')

        return errorMsg

    # def allowsDelete(self, obj):
    #     return True
    #
    def getOutputSetOfTS(self) -> SetOfTiltSeries:
        outputName = self._possibleOutputs.tiltSeries.name
        outputSetOfTiltSeries = getattr(self, outputName, None)

        if outputSetOfTiltSeries:
            outputSetOfTiltSeries.enableAppend()
        else:
            outputSetOfTiltSeries = SetOfTiltSeries.create(self._getPath(), template='tiltseries')
            tsMSetPointer = self.getInputTSMovies(asPointer=True)
            tsMovieSet = tsMSetPointer.get()
            outputSetOfTiltSeries.setSamplingRate(self.outSamplingRate)
            outputSetOfTiltSeries.setAcquisition(tsMovieSet.getAcquisition())
            outputSetOfTiltSeries.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(**{outputName: outputSetOfTiltSeries})
            self._defineSourceRelation(tsMSetPointer, outputSetOfTiltSeries)

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

    def _getBinFactor(self) -> int:
        """2^x pre-binning factor, applied in Fourier space when loading raw data.
        0 = no binning, 1 = 2x2 binning, 2 = 4x4 binning"""
        return math.floor(math.log2(self.binFactor.get()))

    def _getTomoDims(self) -> Tuple[int, int, int]:
        tsMSet = self.getInputTSMovies()
        doSwap = 45 < abs(tsMSet.getAcquisition().getTiltAxisAngle()) < 135
        x = int(tsMSet.getDimensions()[0])
        y = int(tsMSet.getDimensions()[1])
        z = int(max(x, y) / 3)
        if doSwap:
            x, y = y, x
        return x, y, z

    def _genGainAndDarkCmd(self, args: List[str]) -> None:
        if self.gainFile:
            args.append(f'--gain_path {self.gainFile}')
            if self.gainFlip.get() == FLIP_X:
                args.append('--gain_flip_x')
            elif self.gainFlip.get() == FLIP_Y:
                args.append('--gain_flip_y')
            if self.gainSwap.get() == TRANSPOSE_XY:
                args.append('--gain_transpose')
        if self.darkFile:
            args.append(f'--defects_path {self.darkFile}')

    def _addThreadsAndGpusArgs(self, args: List[str]) -> None:
        args.append(f"--perdevice {self.numberOfThreads.get()}")
        gpu_list = self.getGpuList()
        if gpu_list:
            device_ids = " ".join(map(str, gpu_list))
            args.append(f"--device_list {device_ids}")

    def _getAutoHandNTs(self) -> int:
        tsMSet = self.getInputTSMovies()
        nTsM = len(tsMSet)
        return nTsM if nTsM < 10 else 10
