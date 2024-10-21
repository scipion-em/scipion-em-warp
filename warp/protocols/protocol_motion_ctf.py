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
from .protocol_base import ProtMovieAlignBase, ProtWarpBase

from .. import Plugin, CREATE_SETTINGS, FS_MOTION, FRAMESERIES_FOLDER, FRAMESERIES_SETTINGS, AVERAGE_FOLDER


class ProtWarpMotionCorr(ProtMovieAlignBase):
    """ This protocol wraps WarpTools programs.
        Estimate motion in frame series, produce aligned averages
    """

    _label = 'motioncorr'
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
        self._defineStreamingParams(form)
        form.addParallelSection(threads=3, mpi=0)

    def _defineAlignmentParams(self, form):
        form.addHidden(params.GPU_LIST, params.StringParam, default='0',
                       expertLevel=params.LEVEL_ADVANCED,
                       label="Choose GPU IDs",
                       help="Space-separated list of GPU IDs to use for processing. Default: all GPUs in the system."
                            " Warp can use multiple GPUs - in that case"
                            " set to i.e. *0 1 2*.")

        form.addParam('binFactor', params.FloatParam, default=0.,
                      label="Binning factor",
                      help="2^x pre-binning factor, applied in Fourier "
                           "space when loading raw data. 0 = no binning, "
                           "1 = 2x2 binning, 2 = 4x4 binning, supports "
                           "non-integer values")

        form.addParam('eerGroup', params.IntParam, default=40,
                      label='EER fractionation',
                      help="The number of hardware frames to group into one "
                           "fraction. This option is relevant only for Falcon 4 "
                           "movies in the EER format. Fractionate such that each fraction "
                           "has about 0.5 to 1.25 e/A2.")

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

        line.addParam('sumFrame0', params.IntParam, default=1,
                      label='from')
        line.addParam('sumFrameN', params.IntParam, default=0,
                      label='to')

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
        createSettingStep = self._insertFunctionStep(self.createSettingStep, prerequisites=[], needsGPU=False)
        return [createSettingStep]

    def createSettingStep(self):
        """ Create a settings file. """
        movies = self.getInputMovies()
        firstMovie = movies.getFirstItem()
        fileName, extension = os.path.splitext(firstMovie.getFileName())
        folderData = os.path.abspath(os.path.dirname(fileName))
        processingFolder = os.path.abspath(self._getExtraPath(FRAMESERIES_FOLDER))
        sr = firstMovie.getSamplingRate()
        gainPath = os.path.abspath(movies.getGain()) if movies.getGain() else None
        pwutils.makePath(processingFolder)
        argsDict = {
            "--folder_data": folderData,
            "--extension": "*%s" % extension,
            "--folder_processing": processingFolder,
            "--bin": self.binFactor.get(),
            "--angpix": sr,
            # "--exposure": -30.648,  # exposure/frame
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
        gpuList = self.getGpuList()
        if gpuList:
            argsDict['--device_list'] = ' '.join(map(str, gpuList))

        cmd = ' '.join(['%s %s' % (k, v) for k, v in argsDict.items()])
        cmd += ' --averages'

        if self.x.get() and self.y.get() and self.z.get():
            cmd += ' --grid %sx%sx%s' % (self.x.get(), self.y.get(), self.z.get())

        self.runJob(self.getPlugin().getProgram(FS_MOTION), cmd, executable='/bin/bash')

        processingFolder = self._getExtraPath(FRAMESERIES_FOLDER, AVERAGE_FOLDER)
        # Generate a list of output micrograph locations based on the original micrograph names
        micLocations = [os.path.join(processingFolder, os.path.splitext(micNames)[0] + '.mrc')
                        for micNames in micNamesList]

        # Register the output micrographs along with their corresponding output locations
        self.addMicrographs(micNamesList, micLocations)

















