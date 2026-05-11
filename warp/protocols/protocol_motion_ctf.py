# ******************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@gmail.com) [1]
# *              Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [2]
# *
# * [1] St.Jude Children's Research Hospital, Memphis, TN
# * [2] MRC Laboratory of Molecular Biology (MRC-LMB)
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

from pyworkflow import BETA
import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params

from warp.protocols.protocol_base import ProtMovieAlignBase
from warp import Plugin, WARP_TOOLS
from warp.constants import (CREATE_SETTINGS, FS_MOTION, FRAMESERIES_FOLDER,
                            FRAMESERIES_SETTINGS, AVERAGE_FOLDER)


class ProtWarpMotionCorr(ProtMovieAlignBase):
    """
    Performs movie motion correction using WarpTools by estimating beam-induced
    motion in frame series and generating aligned average micrographs.

    AI Generated:

    Motion Correction (ProtWarpMotionCorr) — User Manual
        Overview

        This protocol performs beam-induced motion correction on cryo-EM movie
        frame series using WarpTools.

        Its main purpose is to estimate frame-to-frame motion produced during
        image acquisition and generate aligned average micrographs suitable for
        downstream processing such as CTF estimation, particle picking, and
        high-resolution reconstruction.

        In biological cryo-EM workflows, motion correction is one of the
        earliest and most critical preprocessing steps because uncorrected
        specimen drift and beam-induced motion strongly reduce high-resolution
        information.

        Inputs and General Workflow

        The protocol takes as input a set of raw movies.

        From the input metadata, the protocol extracts:

        - Sampling rate.
        - Per-frame exposure.
        - Gain reference (if available).
        - File format and movie location.

        The workflow follows these stages:

        1. Creation of Warp processing settings.
        2. Optional EER-specific fractionation setup.
        3. Motion estimation on selected movies.
        4. Generation of aligned average micrographs.
        5. Registration of output micrographs.

        Warp Settings Preparation

        Before motion estimation begins, the protocol creates a Warp settings
        file describing the experimental dataset.

        The settings include:

        - Raw movie folder.
        - File extension.
        - Processing directory.
        - Binning factor.
        - Pixel size.
        - Exposure per frame.

        If a gain reference is available, the protocol also applies optional:

        - Gain transpose.
        - Gain flip along X.
        - Gain flip along Y.

        From a practical cryo-EM perspective, correct gain handling is
        essential because a wrong gain orientation can introduce systematic
        artifacts that propagate through the entire workflow.

        Motion Estimation

        During processing, the protocol launches Warp motion estimation on the
        selected input movies.

        The estimation uses several user-controlled parameters:

        Resolution range

        The user defines the minimum and maximum spatial frequencies used to
        estimate motion.

        In practical terms:

        - Lower resolution emphasizes global motion behavior.
        - Higher resolution allows finer local motion fitting.

        B-factor weighting

        High spatial frequencies can be downweighted using a B-factor.

        This often stabilizes motion estimation when movies are noisy.

        Motion model grid

        Motion can be modeled spatially and temporally through a grid in X, Y,
        and temporal dimensions.

        Biologically, this becomes important when local motion differs across
        the field of view, which is common in thin ice, uneven support films,
        or large particles distributed across the micrograph.

        GPU Acceleration

        The protocol supports GPU execution.

        Multiple GPUs can be assigned, which is especially useful in
        high-throughput cryo-EM facility environments where large movie
        datasets must be processed efficiently.

        EER Movie Support

        For EER movies, the protocol provides dedicated fractionation options.

        Two strategies are supported:

        - EER fractionation by grouping frames into virtual frames.
        - Exposure-based grouping.

        Exposure-based grouping overrides fixed frame grouping.

        In biological practice, EER fractionation directly affects temporal
        sampling and therefore the balance between motion correction quality
        and computational cost.

        A practical guideline is to choose fractions corresponding roughly to
        about 0.5–1.25 e-/Å² per grouped frame.

        Output Generation

        After motion correction, the protocol generates aligned average
        micrographs.

        For each input movie:

        - Motion is estimated.
        - Frames are aligned.
        - A corrected average micrograph is written as output.

        These micrographs are then registered as protocol outputs and become
        available for downstream Scipion workflows.

        Biological Interpretation

        The output micrographs should show:

        - Sharper high-resolution Thon rings.
        - Better particle contrast.
        - Reduced blurring caused by drift.

        From a biological perspective, motion correction does not create new
        information, but it preserves high-resolution signal that would
        otherwise be lost during acquisition.

        Practical Recommendations

        In most cryo-EM workflows:

        - Start with default motion grid values.
        - Use moderate B-factor downweighting when data are noisy.
        - Verify gain orientation carefully.
        - Use EER grouping conservatively to preserve temporal information.

        For very stable datasets, coarse grids are often sufficient.

        For difficult datasets with strong local motion, increasing the motion
        model grid may substantially improve alignment quality.

        Validation and Streaming

        The protocol supports streaming execution.

        During streaming, movies are processed in batches.

        A validation check ensures that the batch size is larger than zero.

        This makes the protocol appropriate for on-the-fly cryo-EM data
        acquisition pipelines.

        Output Summary

        After execution, the protocol reports how many micrographs have been
        generated relative to the number of input movies.

        This provides a quick operational overview of processing completion.

        Final Perspective

        Motion correction is one of the most important early preprocessing
        steps in cryo-EM.

        Reliable motion estimation directly impacts all downstream analyses,
        including CTF fitting, particle alignment, and final resolution.

        For most biological users, careful handling of gain reference,
        fractionation strategy, and motion grid selection will usually have
        the strongest practical effect on final data quality.
    """

    _label = 'motion correction'
    _devStatus = BETA
    evenOddCapable = True

    def __init__(self, **kwargs):
        ProtMovieAlignBase.__init__(self, **kwargs)

    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        form.addSection('Input')
        super()._defineInputMoviesParam(form)
        form.addSection('Alignment')
        self._defineAlignmentParams(form)
        form.addSection("EER")
        form.addParam('EERtext', params.LabelParam,
                      label="These options are ignored for non-EER movies.")
        form.addParam('eer_ngroups', params.IntParam, default=16,
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
        self._defineStreamingParams(form)
        form.addParallelSection(threads=3, mpi=0)

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
        line.addParam('range_min', params.FloatParam, default=500,
                      label='Min')
        line.addParam('range_max', params.FloatParam, default=10,
                      label='Max')

        form.addParam('bfactor', params.FloatParam, default=-500,
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

        # form.addParam('average_halves', params.BooleanParam,
        #               default=False,
        #               label='Do even and odd ?',
        #               help='Export aligned averages of odd and even frames separately, e.g. for denoiser training')

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
    def insertInitialSteps(self):
        self.samplingRate = self.getInputMovies().getSamplingRate()
        createSettingStep = self._insertFunctionStep(self.createSettingStep,
                                                     prerequisites=[], needsGPU=False)
        return [createSettingStep]

    def createSettingStep(self):
        """ Create a settings file. """
        movies = self.getInputMovies()
        firstMovie = movies.getFirstItem()
        fileName, extension = os.path.splitext(firstMovie.getFileName())
        folderData = os.path.abspath(os.path.dirname(fileName))
        processingFolder = os.path.abspath(self._getExtraPath(FRAMESERIES_FOLDER))
        sr = firstMovie.getSamplingRate()
        exposure = -1 * movies.getAcquisition().getDosePerFrame()
        gainPath = os.path.abspath(movies.getGain()) if movies.getGain() else None
        pwutils.makePath(processingFolder)
        argsDict = {
            "--folder_data": folderData,
            "--extension": "*%s" % extension,
            "--folder_processing": processingFolder,
            "--bin": self.getBinFactor(),
            "--angpix": sr,
            "--output": os.path.abspath(self._getExtraPath(FRAMESERIES_SETTINGS)),
        }

        if exposure is not None:
            argsDict['--exposure'] = exposure

        if extension == '.eer':
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

    def proccessMoviesStep(self, micNamesList) -> None:
        """Estimate motion in frame series, produce aligned averages and register the output"""
        self.info(">>> Starting estimate motion...")
        # Prepare a list of absolute paths for the movies to process
        # Each movie name in micNamesList is converted to an absolute path and join them into a
        # single string separated by spaces (warp specification)
        warpMoviesNamesList = [os.path.abspath(self.getFileName(micName)) for micName in micNamesList]
        warpMoviesNamesList = " ".join(warpMoviesNamesList)
        argsDict = {
            "--settings": os.path.abspath(self._getExtraPath(FRAMESERIES_SETTINGS)),
            "--range_min": self.range_min.get(),
            "--range_max": self.range_max.get(),
            "--bfac": self.bfactor.get(),
            "--input_data": warpMoviesNamesList
        }
        gpuList = self._stepsExecutor.getGpuList()
        if gpuList:
            argsDict['--device_list'] = ' '.join(map(str, gpuList))

        cmd = ' '.join(['%s %s' % (k, v) for k, v in argsDict.items()])
        cmd += ' --averages'

        if self.x.get() and self.y.get() and self.z.get():
            cmd += ' --grid %sx%sx%s' % (self.x.get(), self.y.get(), self.z.get())

        self.runJob(self.getPlugin().getProgram(WARP_TOOLS, FS_MOTION), cmd, executable='/bin/bash')

        processingFolder = self._getExtraPath(FRAMESERIES_FOLDER, AVERAGE_FOLDER)
        # Generate a list of output micrograph locations based on the original micrograph names
        micLocations = [os.path.join(processingFolder, os.path.splitext(micNames)[0] + '.mrc')
                        for micNames in micNamesList]

        # Register the output micrographs along with their corresponding output locations
        self.addMicrographs(micNamesList, micLocations)

    def _summary(self):
        summary = []
        if self.hasAttribute(self.OUT_MICS):
            summary.append(f"Micrographs: {self.Micrographs.getSize()} of {self.inputMovies.get().getSize()}\n")
        else:
            summary.append("Outputs are not ready yet.")
        return summary

    def _validate(self):
        errors = []
        if self.streamingBatchSize.get() < 1:
            errors.append('The batch size value must be greater than 1')
        return errors

    def getBinFactor(self):
        import math
        return math.floor(math.log2(self.binFactor.get()))