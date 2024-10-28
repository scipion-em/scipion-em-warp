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
from pyworkflow import BETA
import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from tomo.objects import SetOfTiltSeriesM
from .protocol_motion_ctf import ProtWarpMotionCorr
from .protocol_base import ProtMovieAlignBase

from .. import Plugin, CREATE_SETTINGS, FS_MOTION, FRAMESERIES_FOLDER, FRAMESERIES_SETTINGS, AVERAGE_FOLDER


class ProtWarpTSMotionCorr(ProtWarpMotionCorr):
    """ This protocol wraps WarpTools programs.
        Align tilt-series movies
    """

    _label = 'align tilt-series movies'
    _devStatus = BETA
    evenOddCapable = True

    def __init__(self, **kwargs):
        ProtWarpMotionCorr.__init__(self, **kwargs)

    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputMovies', params.PointerParam, pointerClass=SetOfTiltSeriesM,
                      important=True,
                      label=pwutils.Message.LABEL_INPUT_MOVS,
                      help='Select a set of previously imported tilt series movies.')
        form.addSection('Alignment')
        self._defineAlignmentParams(form)
        self._defineStreamingParams(form)
        form.addParallelSection(threads=3, mpi=0)

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

        # processingFolder = self._getExtraPath(FRAMESERIES_FOLDER, AVERAGE_FOLDER)
        # # Generate a list of output micrograph locations based on the original micrograph names
        # micLocations = [os.path.join(processingFolder, os.path.splitext(micNames)[0] + '.mrc')
        #                 for micNames in micNamesList]
        #
        # # Register the output micrographs along with their corresponding output locations
        # self.addMicrographs(micNamesList, micLocations)


















