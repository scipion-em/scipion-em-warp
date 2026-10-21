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
from collections import Counter
from enum import Enum
from os.path import splitext, abspath, dirname, join, basename, exists
from typing import Union, Tuple, List, Optional
from pwem.emlib.image.image_readers import ImageStack, ImageReadersRegistry
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.object import Set, Pointer
from pyworkflow.protocol import GPU_LIST, PointerParam, StringParam, LEVEL_ADVANCED, EnumParam, \
    FloatParam, IntParam, BooleanParam, LabelParam
from pyworkflow.utils import cyanStr, Message, makePath, redStr, replaceBaseExt, yellowStr
from tomo.objects import (SetOfTiltSeriesM, SetOfTiltSeries, TiltImage,
                          TiltSeries, SetOfCTFTomoSeries,
                          TiltSeriesM, TiltImageM, CTFTomoSeries, CTFTomo)
from warp import Plugin, WARP_TOOLS, CREATE_SETTINGS, TS_CTF, TS_DEFOCUS_HAND, FRAMESERIES_FOLDER, \
    AVERAGE_FOLDER, FRAMESERIES_SETTINGS, TOMOSTAR_FOLDER, TILTSERIES_SETTINGS, TOMOSTAR_EXT, TILTSERIES_FOLDER, \
    MRCS_EXT, EVEN, ODD, POWERSPECTRUM_FOLDER, FS_MOTION, FS_CTF
from warp.convert import writeTsStar
from warp.utils import parseCtfXMLFile, _getDefocusHand

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

# Fixed motion-correction patch size used by Warp fs_motion (WarpLib Movie.ProcessShift.cs).
# If binning shrinks the frame below this on either axis, the patch grid degenerates and the
# native worker crashes with a cuFFT CUFFT_INVALID_SIZE error.
MOTION_WINDOW = 768


class WarpTsMcorrOutputs(Enum):
    tiltSeries = SetOfTiltSeries
    ctfs = SetOfCTFTomoSeries


class ProtWarpTSMotionCorr(EMProtocol):
    """ This protocol wraps WarpTools programs.
        Estimate motion in frame series, produce aligned averages, estimate CTF
    """

    _label = 'tilt-series motion correction and ctf estimation'
    _devStatus = BETA
    _possibleOutputs = WarpTsMcorrOutputs

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
        self.failedTsIds = []
        # Per-tilt-series defocus handedness (+1 / -1)
        self.tsHandednessDict = {}

    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        form.addSection(Message.LABEL_INPUT)
        form.addParam('inputTSMovies', PointerParam, pointerClass=SetOfTiltSeriesM,
                      important=True,
                      label=Message.LABEL_INPUT_MOVS,
                      help='Select a set of previously imported tilt series movies.')
        self._defineAlignmentParams(form)
        self._defineCtfParams(form)
        self._defineEerParams(form)
        form.addParallelSection(threads=2, mpi=0)
        form.addHidden(GPU_LIST, StringParam, default='0',
                       expertLevel=LEVEL_ADVANCED,
                       label="Choose GPU IDs",
                       help="Space-separated list of GPU IDs to use for processing. Default: all GPUs in the system."
                            " Warp can use multiple GPUs - in that case"
                            " set to i.e. *0 1 2*.")

    @staticmethod
    def _defineAlignmentParams(form):
        form.addSection('Alignment')
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
                      label='Do even and odd?',
                      help='Export aligned averages of odd and even frames separately, e.g. for '
                           'noise2noise-based denoisers.')

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

    @staticmethod
    def _defineEerParams(form):
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

    @staticmethod
    def _defineCtfParams(form):
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

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._initialize()
        closeSetStepDeps = []
        pId = self._insertFunctionStep(self.createTsMovSettingsStep,
                                       prerequisites=[],
                                       needsGPU=False)

        pId = self._insertFunctionStep(self.createTsSettingsStep,
                                       prerequisites=pId,
                                       needsGPU=False)

        for tsId, tsM in self.tsMDict.items():
            pId = self._insertFunctionStep(self.motionCorrStep, tsId,
                                           prerequisites=pId,
                                           needsGPU=True)

            pId = self._insertFunctionStep(self.createTsMStarStep, tsId,
                                           prerequisites=pId,
                                           needsGPU=False)

            pId = self._insertFunctionStep(self.estimateTsCtfStep, tsId,
                                           prerequisites=pId,
                                           needsGPU=True)

            pId = self._insertFunctionStep(self.estimateHandednessStep, tsId,
                                           prerequisites=pId,
                                           needsGPU=True)

            pId = self._insertFunctionStep(self.createOutputTsStep, tsId,
                                           prerequisites=pId,
                                           needsGPU=False)

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
        self.outSamplingRate = self._getOutputSamplingRate()
        self.exposure = tsMSet.getAcquisition().getDosePerFrame()
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
            # Second, unbinned settings used only for the frame-series CTF, so that the CTF
            # (and the defocus handedness it feeds) is estimated at full resolution regardless
            # of the motion-correction binning.
            if self.estimateCTF.get():
                logger.info(cyanStr(">>> Creating unbinned frame-series CTF settings..."))
                cmd = self._genCreateSettingsArgs(settingsFn=self._getFrameSeriesCtfSettingsFn(),
                                                  binningArg="--bin 0")
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

    def motionCorrStep(self, tsId: str):
        """Motion correction (+ aligned averages) at the requested binning."""
        try:
            logger.info(cyanStr(f">>> tsId = {tsId} - performing the motion-correction..."))
            tsMovie = self.tsMDict[tsId]
            warpMoviesNamesList = " ".join([abspath(ti.getFileName()) for ti in tsMovie.iterItems()])
            args = self._genMcArgs(tsMovie, warpMoviesNamesList)
            self.runJob(Plugin.getProgram(WARP_TOOLS, FS_MOTION), args, executable='/bin/bash')

        except Exception as e:
            logger.error(redStr(f"{WARP_TOOLS} {FS_MOTION} failed with the exception --> {e}"))
            traceback.print_exc()
            self.failedTsIds.append(tsId)
            return

    def estimateTsMCtfParams(self, tsId: str):
        """Frame-series CTF, estimated UNBINNED (decoupled from the motion binning) so that the
        defocus gradient and handedness stay reliable. Writes physical CTF into the same movie
        XMLs as the motion step."""
        if not self.estimateCTF.get() or tsId in self.failedTsIds:
            return

        try:
            logger.info(cyanStr(f">>> tsId = {tsId} - estimating the frame-series CTF (unbinned)..."))
            tsMovie = self.tsMDict[tsId]
            warpMoviesNamesList = " ".join([abspath(ti.getFileName()) for ti in tsMovie.iterItems()])
            args = self._genFsCtfArgs(tsMovie, warpMoviesNamesList)
            self.runJob(Plugin.getProgram(WARP_TOOLS, FS_CTF), args, executable='/bin/bash')

        except Exception as e:
            logger.error(redStr(f"{WARP_TOOLS} {FS_CTF} failed with the exception --> {e}"))
            traceback.print_exc()
            self.failedTsIds.append(tsId)

    def createTsMStarStep(self, tsId: str):
        if not self.estimateCTF.get() or tsId in self.failedTsIds:
            return
        try:
            logger.info(cyanStr(f">>> {tsId} - Creating the star file..."))
            writeTsStar(self.tsMDict[tsId], self._getOutTsStarFile(tsId))

        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to create the TS star file '
                                f'with exception {e}. Skipping... '))
            logger.error(traceback.format_exc())
            self.failedTsIds.append(tsId)

    def estimateTsCtfStep(self, tsId: str):
        if not self.estimateCTF.get() or tsId in self.failedTsIds:
            return

        try:
            # CTF estimation for this single tilt-series (.tomostar)
            logger.info(cyanStr(f">>> tsId = {tsId} - estimating the tilt-series CTF..."))
            args = self._getCtfEstimationArgs(tsId)
            self.runJob(Plugin.getProgram(WARP_TOOLS, TS_CTF), args, executable='/bin/bash')

        except Exception as e:
            logger.error(redStr(f"tsId = {tsId} -> {WARP_TOOLS} {TS_CTF} failed with the exception --> {e}"))
            logger.error(traceback.format_exc())
            self.failedTsIds.append(tsId)

    def estimateHandednessStep(self, tsId: str):
        if not self.estimateCTF.get() or tsId in self.failedTsIds:
            return

        try:
            # Defocus handedness for this single tilt-series
            logger.info(cyanStr(f">>> tsId = {tsId} - estimating the defocus handedness..."))
            args = [
                f'--settings {self._getTiltSeriesSettingsFn()}',
                f'--input_data {self._getOutTsStarFile(tsId)}',
                '--set_auto',
            ]
            self.runJob(Plugin.getProgram(WARP_TOOLS, TS_DEFOCUS_HAND), ' '.join(args), executable='/bin/bash')
            xmlPath = join(self._getTiltSeriesDir(), tsId + '.xml')
            self.tsHandednessDict[tsId] = _getDefocusHand(xmlPath)

        except Exception as e:
            logger.warning(yellowStr(f"tsId = {tsId} -> {WARP_TOOLS} {TS_DEFOCUS_HAND} failed with the exception "
                                     f"--> {e}. This tilt-series will not contribute to the handedness vote."))
            logger.warning(traceback.format_exc())

    def createOutputTsStep(self, tsId: str):
        if tsId in self.failedTsIds:
            return

        try:
            logger.info(cyanStr(f">>> tsId = {tsId} - creating the output tilt-series..."))
            newTs, tiList = self._prepareOutputTs(tsId)
            perTsHand = None
            majorityHand = None
            if self.estimateCTF.get():
                # Per-tilt-series handedness (this series' own XML) for the TS, and the
                # dataset-wide majority vote for the output set metadata object.
                xmlPath = join(self._getTiltSeriesDir(), tsId + '.xml')
                if exists(xmlPath):
                    perTsHand = _getDefocusHand(xmlPath)
                majorityHand = self._getMajorityHandedness()
            # Fall back to the majority if this single series' value is unavailable
            newTs.setHandedness(perTsHand if perTsHand is not None else majorityHand)
            self._registerOutputTs(newTs, tiList, majorityHand)

        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to register the output TS with exception {e}. Skipping... '))
            logger.error(traceback.format_exc())
            self.failedTsIds.append(tsId)

    def createOutputCtfStep(self, tsId: str):
        if not self.estimateCTF.get() or tsId in self.failedTsIds:
            return

        try:
            logger.info(cyanStr(f">>> tsId = {tsId} - creating the output CTF..."))
            newCtfTomoSeries, ctfList = self._prepareOutputCtf(tsId)
            self._registerOutputCtf(newCtfTomoSeries, ctfList)

        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to register the output CTF with exception {e}. Skipping... '))
            logger.error(traceback.format_exc())
            self.failedTsIds.append(tsId)

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

    def _getFrameSeriesCtfSettingsFn(self) -> str:
        """Unbinned frame-series settings used only for the frame CTF, so that the
        CTF (and hence the defocus handedness) is estimated at full resolution
        regardless of the motion-correction binning."""
        return abspath(self._getExtraPath('warp_frameseries_ctf.settings'))

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

    def _genBinningArgs(self):
        binMode = self.binFactorMode.get()
        binTarget = self.binTarget.get()
        return f"--bin {self._getBinFactor()}" if binMode == BINNING_FACTOR else f"--bin_angpix {binTarget}"

    def _genCreateSettingsArgs(self, settingsFn: Optional[str] = None,
                               binningArg: Optional[str] = None) -> str:
        ts_movies = self.getInputTSMovies()
        first_ts_movie = ts_movies.getFirstItem()
        file_name, extension = splitext(first_ts_movie.getFirstItem().getFileName())
        folder_data = abspath(dirname(file_name))

        # Defaults: the main (possibly binned) frame-series settings.
        if settingsFn is None:
            settingsFn = self._getFrameSeriesSettingsFn()
        if binningArg is None:
            binningArg = self._genBinningArgs()

        # Initialize argument list with base parameters
        args = [
            f"--folder_data {folder_data}",
            f"--extension '*{extension}'",
            f"--folder_processing {self.tsmDir}",
            binningArg,
            f"--angpix {self.samplingRate}",
            f"--exposure {self.exposure}",
            f"--output {settingsFn}",
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
            # self._genBinningArgs(),
            f'--exposure {self.exposure}',
            f'--tomo_dimensions {x}x{y}x{z}'
        ]
        self._genGainAndDarkCmd(args)
        return ' '.join(args)

    def _genMcArgs(self, tsMovie: TiltSeriesM, warpMoviesNames: str) -> str:
        """fs_motion arguments: motion correction + aligned averages at the requested binning.
        The frame-series CTF is handled separately (unbinned) by _genFsCtfArgs."""
        commonArgs = [
            f"--settings {self._getFrameSeriesSettingsFn()}",
            f"--input_data {warpMoviesNames}",
            f"--output_processing {self.tsmDir}",
        ]

        # GPU and threads configuration
        self._addThreadsAndGpusArgs(commonArgs)

        args = [f"--range_min {self.m_range_min.get()}",
                f"--range_max {self.m_range_max.get()}",
                f"--bfac {self.bfactor.get()}",
                '--averages']

        # Motion grid dimensions
        grid_x, grid_y, grid_z = self.x.get(), self.y.get(), self.z.get()
        if grid_x and grid_y and grid_z:
            commonArgs.append(f"--grid {grid_x}x{grid_y}x{grid_z}")
        if self.average_halves.get():
            commonArgs.append("--average_halves")
        return " ".join(commonArgs + args)

    def _genFsCtfArgs(self, tsMovie: TiltSeriesM, warpMoviesNames: str) -> str:
        """Frame-series CTF (fs_ctf) on the UNBINNED settings. Decoupling the CTF from the
        motion-correction binning keeps the local defocus gradient (and the defocus
        handedness read by ts_defocus_hand) reliable. Writes physical CTF values into the
        same movie XMLs as the motion step, which the tilt-series steps then consume."""
        acq = tsMovie.getAcquisition()
        args = [
            f"--settings {self._getFrameSeriesCtfSettingsFn()}",
            f"--input_data {warpMoviesNames}",
            f"--output_processing {self.tsmDir}",
            f"--window {self.window.get()}",
            f"--range_min {self.range_min.get()}",
            f"--range_max {self.range_max.get()}",
            f"--defocus_min {self.defocus_min.get()}",
            f"--defocus_max {self.defocus_max.get()}",
            f"--voltage {int(acq.getVoltage())}",
            f"--cs {acq.getSphericalAberration()}",
            f"--amplitude {acq.getAmplitudeContrast()}",
        ]
        self._addThreadsAndGpusArgs(args)

        # Defocus model grid. ts_defocus_hand needs a local defocus gradient, i.e. a grid of
        # at least 2x2x1, so fall back to that minimum when no explicit grid is given.
        c_grid_x, c_grid_y, c_grid_z = self.c_x.get(), self.c_y.get(), self.c_z.get()
        if c_grid_x and c_grid_y and c_grid_z:
            args.append(f"--grid {c_grid_x}x{c_grid_y}x{c_grid_z}")
        else:
            args.append("--grid 2x2x1")

        if self.fit_phase.get():
            args.append("--fit_phase")
        if self.use_sum.get():
            args.append("--use_sum")

        return " ".join(args)

    def _getCtfEstimationArgs(self, tsId: str) -> str:
        acq = self.getInputTSMovies().getAcquisition()
        args = [
            f'--settings {self._getTiltSeriesSettingsFn()}',
            f'--input_data {self._getOutTsStarFile(tsId)}',
            f'--window {self.window.get()}',
            f'--range_low {self.range_min.get()}',
            f'--range_high {self.range_max.get()}',
            f'--defocus_min {self.defocus_min.get()}',
            f'--defocus_max {self.defocus_max.get()}',
            f'--voltage {acq.getVoltage():.0f}',
            f'--cs {acq.getSphericalAberration()}',
            f'--amplitude {acq.getAmplitudeContrast()}',
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

        inTiList = [ti.clone() for ti in tsMovie.iterItems(orderBy=TiltImageM.TILT_ANGLE_FIELD)]
        tiList = []
        for i, tiM in enumerate(inTiList):
            newTi = TiltImage()
            newTi.copyInfo(tiM)
            newTi.setFileName(newBinaryName)
            newTi.setIndex(i + 1)
            newTi.setSamplingRate(self.outSamplingRate)
            newTi.setAcquisition(tiM.getAcquisition())

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

    def _registerOutputTs(self,
                          newTs: TiltSeries,
                          tiList: List[TiltImage],
                          defocusHand: Optional[int]):
        # TS set
        outTsSet = self.getOutputSetOfTS(defocusHand)
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

    def _prepareOutputCtf(self, tsId: str) -> Tuple[CTFTomoSeries, List[CTFTomo]]:
        tsDir = self._getTiltSeriesDir()
        psdStack = join(tsDir, POWERSPECTRUM_FOLDER, tsId + '.mrc')
        outTsSet = self._getOutputTsSet()
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

        return newCTFTomoSeries, ctfTomoList

    def _registerOutputCtf(self, newCTFTomoSeries: CTFTomoSeries, ctfTomoList: List[CTFTomo]):
        # CTFTomoSeries set
        outputCtfSet = self.getOutputSetOfCTFTomoSeries()
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

    def getOutputSetOfTS(self, defocusHand: Optional[int]) -> SetOfTiltSeries:
        outputName = self._possibleOutputs.tiltSeries.name
        outputSetOfTiltSeries = self._getOutputTsSet()

        if outputSetOfTiltSeries:
            outputSetOfTiltSeries.enableAppend()
        else:
            outputSetOfTiltSeries = SetOfTiltSeries.create(self._getPath(), template='tiltseries')
            tsMSetPointer = self.getInputTSMovies(asPointer=True)
            tsMovieSet = tsMSetPointer.get()
            outputSetOfTiltSeries.setSamplingRate(self.outSamplingRate)
            outputSetOfTiltSeries.setAcquisition(tsMovieSet.getAcquisition())
            outputSetOfTiltSeries.setHandedness(defocusHand)
            self._defineOutputs(**{outputName: outputSetOfTiltSeries})
            self._defineSourceRelation(tsMSetPointer, outputSetOfTiltSeries)

        return outputSetOfTiltSeries

    def getOutputSetOfCTFTomoSeries(self):
        outputSetName = self._possibleOutputs.ctfs.name
        outputSetOfCTFTomoSeries = self._getOutputCtfSet()

        if outputSetOfCTFTomoSeries:
            outputSetOfCTFTomoSeries.enableAppend()
        else:
            outputSetOfCTFTomoSeries = SetOfCTFTomoSeries.create(self._getPath(),
                                                                 template='CTFmodels%s.sqlite')

            tsSet = self._getOutputTsSet()
            outputSetOfCTFTomoSeries.setSetOfTiltSeries(tsSet)
            self._defineOutputs(**{outputSetName: outputSetOfCTFTomoSeries})
            self._defineCtfRelation(outputSetOfCTFTomoSeries, tsSet)

        return outputSetOfCTFTomoSeries

    def _getBinFactor(self) -> int:
        """2^x pre-binning factor, applied in Fourier space when loading raw data.
        0 = no binning, 1 = 2x2 binning, 2 = 4x4 binning"""
        return math.floor(math.log2(self.binFactor.get()))

    def _getTomoDims(self) -> Tuple[int, int, int]:
        tsMSet = self.getInputTSMovies()
        x = int(tsMSet.getDimensions()[0])
        y = int(tsMSet.getDimensions()[1])
        # The binary swap is exact only near φ∈{0,90}. Near 45°, the detector rectangle rotated into
        # the tomogram frame needs a box larger than either W×H or H×W. The exact axis-aligned bounding
        # box of the rotated rectangle covers every imaged pixel at any angle and self-pads near 45°I
        tiltAxisAngle = tsMSet.getAcquisition().getTiltAxisAngle()
        c, s = abs(math.cos(tiltAxisAngle)), abs(math.sin(tiltAxisAngle))
        tomoX = math.ceil(x * c + y * s)
        tomoY = math.ceil(x * s + y * c)
        tomoZ = int(max(tomoX, tomoY) / 3)
        return tomoX, tomoY, tomoZ

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

    def _getMajorityHandedness(self) -> Optional[int]:
        """Majority vote of the per-tilt-series defocus handedness (+1 / -1).

        Uses the values recorded during CTF estimation; on a continued execution
        (where those steps are skipped) it falls back to reading each ``<tsId>.xml``.
        Returns None if no handedness could be determined. A tie falls back to -1
        (AreAnglesInverted = False, Warp's 'no flip' default)."""
        votes = list(self.tsHandednessDict.values())
        if not votes:  # continued run: recover from the persisted XMLs
            tsDir = self._getTiltSeriesDir()
            for tsId in self.tsMDict:
                if tsId in self.failedTsIds:
                    continue
                xmlPath = join(tsDir, tsId + '.xml')
                if not exists(xmlPath):
                    continue
                try:
                    votes.append(_getDefocusHand(xmlPath))
                except Exception:
                    continue

        if not votes:
            return None

        counts = Counter(votes)
        plus, minus = counts.get(1, 0), counts.get(-1, 0)
        logger.info(cyanStr(f">>> Defocus handedness vote across {len(votes)} tilt-series: "
                            f"flip(+1)={plus}, no-flip(-1)={minus}"))
        if plus == minus:
            logger.warning(yellowStr(">>> Handedness vote tied; defaulting to -1 (no flip)."))
            return -1
        return 1 if plus > minus else -1

    def _getOutputTsSet(self) -> Optional[SetOfTiltSeries]:
        return getattr(self, self._possibleOutputs.tiltSeries.name, None)

    def _getOutputCtfSet(self) -> Optional[SetOfCTFTomoSeries]:
        return getattr(self, self._possibleOutputs.ctfs.name, None)

    def _getOutputSamplingRate(self) -> float:
        binMode = self.binFactorMode.get()
        binFactor = self.binFactor.get()
        binTarget = self.binTarget.get()
        movieSampling = self.inputTSMovies.get().getSamplingRate()
        return movieSampling * binFactor if binMode == BINNING_FACTOR else binTarget

    # -------------------------- INFO functions -------------------------------
    def _summary(self):
        summary = []
        if self.isFinished():
            tilseriesSize = 0
            ctfSize = 0
            inputMoviesSize = len(self.inputTSMovies.get())
            outTsSet = self._getOutputTsSet()
            outCtfSet = self._getOutputCtfSet()
            if outTsSet:
                tilseriesSize = len(outTsSet)
            summary.append(f"Aligned tiltseries: {tilseriesSize} of {inputMoviesSize}")

            if outCtfSet:
                ctfSize = len(outCtfSet)
            summary.append(f"CTF estimated: {ctfSize} of {inputMoviesSize}")

            if self.estimateCTF.get():
                if not outCtfSet:
                    summary.append("*Error: No CTF estimated!*\n"
                                   "Check the Output Log > run.stderr to check the error.")
                elif outCtfSet and tilseriesSize > ctfSize:
                    summary.append("*Warning: some CTF were not estimated!*\n"
                                   "Check the Output Log > run.stderr to check the details.")

        return summary

    def _validate(self):
        errorMsg = []
        binMode = self.binFactorMode.get()
        binFactor = self.binFactor.get()
        binTarget = self.binTarget.get()
        outputSampling = self._getOutputSamplingRate()
        nyquistFreq = 2 * outputSampling
        resToFit = self.range_max.get()

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

        if self.estimateCTF.get() and resToFit < nyquistFreq:
            warnMs = (
                f'The resolution to fit should be greater than nyquist. Currently, resolution to fit is {resToFit} '
                f'and at this binning factor Nyquist is {nyquistFreq}\n')
            errorMsg.append(warnMs)

        # Warp tiles the frame with fixed windows: fs_motion uses a MOTION_WINDOW (768) patch, and
        # when the CTF is also estimated (fs_motion_and_ctf) it additionally tiles with a 2 x window
        # region. If the binned frame falls below the required window on either axis, the patch grid
        # degenerates and the native worker crashes with a cuFFT CUFFT_INVALID_SIZE error.
        binParamsOk = ((binMode == BINNING_FACTOR and binFactor is not None and binFactor >= 1) or
                       (binMode == TARGET_SAMPLING_RATE and binTarget is not None and binTarget > 0))
        if binParamsOk:
            tsMSet = self.inputTSMovies.get()
            movieSampling = tsMSet.getSamplingRate()
            dims = tsMSet.getDimensions()
            if dims is not None and movieSampling:
                rawX, rawY = int(dims[0]), int(dims[1])
                # Downsample factor as Warp actually applies it: 2^floor(log2(binFactor)) in
                # binning-factor mode, or target/raw pixel size in target-sampling mode.
                if binMode == BINNING_FACTOR:
                    downsample = 2 ** self._getBinFactor()
                else:
                    downsample = binTarget / movieSampling
                if downsample > 0:
                    binnedX = round(rawX / downsample)
                    binnedY = round(rawY / downsample)
                    # The binned frame is only used for motion correction (fs_motion); the CTF is
                    # estimated separately on the unbinned frame, so only the motion window applies.
                    if min(binnedX, binnedY) < MOTION_WINDOW:
                        maxBinAngpix = movieSampling * min(rawX, rawY) / MOTION_WINDOW
                        errorMsg.append(
                            f"Binning is too strong: the binned frame would be {binnedX}x{binnedY} px, "
                            f"smaller than the {MOTION_WINDOW}-px window required by Warp {FS_MOTION}. "
                            f"This makes it crash with a cuFFT CUFFT_INVALID_SIZE error.\nReduce the "
                            f"binning so the binned frame stays >= {MOTION_WINDOW} px on both axes, i.e. "
                            f"use an output pixel size <= {maxBinAngpix:.2f} A/px (raw frame "
                            f"{rawX}x{rawY} px at {movieSampling:.4f} A/px).\n")

        return errorMsg


