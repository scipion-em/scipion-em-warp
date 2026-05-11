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
    Performs motion correction on cryo-EM movie frame series in order
    to compensate for beam-induced specimen drift and produce aligned
    micrograph averages suitable for downstream image processing.

    AI Generated:

    Motion Correction (ProtWarpMotionCorr) - User Manual
        Overview

        The Motion Correction protocol estimates and compensates for
        sample motion present during movie acquisition in cryo-electron
        microscopy. Its main purpose is to align the successive frames
        of each recorded movie so that the final averaged micrograph
        preserves as much high-resolution structural information as
        possible.

        In modern cryo-EM workflows, movies are collected instead of
        single exposure images because beam interaction with the frozen
        specimen induces subtle but significant movement over time.
        Without motion correction, this drift causes image blurring,
        weakens high-frequency signal, and limits the quality of all
        subsequent analysis steps such as contrast transfer estimation,
        particle picking, classification, and high-resolution
        reconstruction.

        Biological Importance

        For biological users, motion correction is one of the earliest
        and most critical preprocessing stages. A well-corrected
        micrograph retains structural detail that may otherwise be lost
        permanently at the beginning of the workflow. This directly
        affects the interpretability of macromolecular complexes,
        membrane proteins, viral particles, and other fragile
        biological assemblies.

        Even small residual motion can reduce the visibility of fine
        structural features. For this reason, careful motion correction
        often has a disproportionate impact on the final attainable
        resolution of a cryo-EM project.

        Inputs and General Workflow

        The protocol operates on a set of input movies acquired as
        frame series. Each movie is treated as a temporal record of the
        same exposure, and the objective is to estimate how the image
        shifts during acquisition and to compensate for those
        displacements before generating the final aligned average.

        The input dataset should ideally have consistent acquisition
        conditions, including stable pixel size, dose information, and
        detector geometry. Accurate acquisition metadata helps the
        refinement process preserve physically meaningful motion
        estimates.

        Motion Modeling

        Motion is modeled across both space and time. This means that
        the protocol can account not only for global specimen drift but
        also for local differential motion across the field of view.
        Such local behavior is common in cryo-EM, particularly in thin
        ice, large fields of view, or specimens with uneven support
        properties.

        For many biological datasets, local motion correction is more
        effective than global alignment because different regions of
        the image may move in slightly different ways during exposure.
        Correctly accounting for this behavior helps preserve local
        structural detail that would otherwise be smeared out.

        Resolution Range and Frequency Weighting

        The refinement can be guided by a chosen spatial frequency
        range. This allows the alignment to focus on the signal most
        informative for motion estimation while reducing sensitivity to
        noise or irrelevant low-frequency intensity variation.

        Frequency weighting is also important. High-resolution
        information is often weaker and more noise sensitive, so
        downweighting unstable frequencies can improve alignment
        robustness. For biological specimens with weak contrast, this
        often leads to more stable and reliable corrected averages.

        Binning and Practical Tradeoffs

        The protocol allows Fourier-space binning of the raw movie
        frames before alignment. Binning can reduce computational cost
        and may improve robustness in noisy datasets. This is often
        useful during exploratory processing, screening sessions, or
        large-scale facility pipelines.

        For high-resolution projects, however, users usually prefer
        minimal binning so that fine structural information remains
        available. The appropriate choice depends on data quality,
        particle size, and the biological resolution goals of the
        experiment.

        EER Movie Support

        The protocol also supports electron event representation movie
        formats. In this context, temporal fractionation becomes
        biologically relevant because it determines how finely the
        exposure is divided during motion estimation.

        Finer temporal sampling can better capture rapid motion early
        in the exposure, where beam-induced drift is often strongest.
        However, excessive temporal subdivision may reduce the signal
        available in each fraction. In practical biological work, the
        best choice is usually a balance between temporal precision and
        sufficient per-frame signal.

        Gain Reference Considerations

        Proper detector gain handling is essential because detector
        normalization artifacts can propagate directly into aligned
        micrographs. When gain orientation differs from the movie
        orientation, correcting that mismatch ensures that intensity
        normalization remains physically meaningful.

        Although this is often considered a technical preprocessing
        detail, its biological consequence is important because
        inaccurate normalization can subtly degrade particle contrast
        and affect all downstream interpretation.

        Outputs and Their Interpretation

        The main output of the protocol is a set of aligned micrograph
        averages. These corrected images serve as the standard starting
        point for subsequent cryo-EM analysis.

        From a practical biological perspective, successful motion
        correction is often recognized by sharper particle boundaries,
        improved visibility of high-frequency features, and more stable
        downstream contrast transfer estimation. These corrected
        micrographs are generally the images that will be inspected for
        data quality and used throughout the remainder of the workflow.

        Practical Recommendations

        In routine cryo-EM processing, it is often advisable to begin
        with conservative alignment settings and visually inspect the
        corrected micrographs. If the dataset shows strong local drift,
        increasing the flexibility of the motion model may improve
        results.

        For noisy datasets or rapid screening, moderate binning can
        provide a useful balance between speed and reliability. For
        high-resolution biological studies, users generally favor
        finer temporal sampling and reduced binning whenever signal
        quality permits.

        Final Perspective

        For most cryo-EM users, motion correction is not simply a
        technical preprocessing task. It is the stage where raw movie
        data are converted into structurally meaningful images.
        Careful handling of motion, detector normalization, temporal
        sampling, and spatial modeling strongly influences the quality
        of every downstream biological interpretation.
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