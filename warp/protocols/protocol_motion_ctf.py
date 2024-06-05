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
import threading
import time
from uuid import uuid4
from collections import OrderedDict
from emtools.utils import Timer, Pipeline

from pyworkflow import SCIPION_DEBUG_NOCLEAN, BETA
import pyworkflow.object as pwobj
import pyworkflow.utils as pwutils
from pyworkflow.protocol.constants import STEPS_SERIAL
import pyworkflow.protocol.params as params
from pwem.objects import Float, SetOfMovies
from pwem.protocols import ProtAlignMovies

from .. import Plugin


class ProtMotionCorrTasks(ProtAlignMovies):
    """ This protocol wraps WarpTools programs.
    """

    _label = 'tasks'
    _devStatus = BETA
    evenOddCapable = True

    def __init__(self, **kwargs):
        ProtAlignMovies.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_SERIAL
        # Disable parallelization options just take into account GPUs
        self.numberOfMpi.set(0)
        self.numberOfThreads.set(0)
        self.allowMpi = False
        self.allowThreads = False

    # We are not using the steps mechanism for parallelism from Scipion
    def _stepsCheck(self):
        pass

    # -------------------------- DEFINE param functions -----------------------
    def _defineAlignmentParams(self, form):
        form.addHidden(params.GPU_LIST, params.StringParam, default='0',
                       expertLevel=params.LEVEL_ADVANCED,
                       label="Choose GPU IDs",
                       help="GPU may have several cores. Set it to zero"
                            " if you do not know what we are talking about."
                            " First core index is 0, second 1 and so on."
                            " Warp can use multiple GPUs - in that case"
                            " set to i.e. *0 1 2*.")

        form.addParam('binFactor', params.FloatParam, default=0.,
                      label="'Binning factor",
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

        form.addParam('bfactor', params.IntParam, default=500,
                      label="B-factor",
                      help="Downweight higher spatial frequencies using a "
                           "B-factor, in Angstrom^2")

        line = form.addLine('Motion model grid',
                            help="Resolution of the motion model grid in "
                                 "X, Y, and temporal dimensions, e.g. 5x5x40; "
                                 "0 = auto")
        line.addParam('x', params.IntParam, default=5, label='X')
        line.addParam('y', params.IntParam, default=5, label='Y')
        line.addParam('z', params.IntParam, default=0, label='Temporal')

        line.addParam('sumFrame0', params.IntParam, default=1,
                      label='from')
        line.addParam('sumFrameN', params.IntParam, default=0,
                      label='to')

        form.addParam('splitEvenOdd', params.BooleanParam,
                      default=False,
                      label='Split & sum odd/even frames?',
                      help='Generate odd and even sums using odd and even frames '
                           'respectively when this option is enabled.')

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

        form.addParam('defectFile', params.FileParam, allowsNull=True,
                      label='Camera defects file',
                      help='')

        self._defineStreamingParams(form)

        # Make default 1 minute for sleeping when no new input movies
        form.getParam('streamingSleepOnWait').setDefault(60)

    # --------------------------- STEPS functions -----------------------------
    def _insertAllSteps(self):
        self.samplingRate = self.getInputMovies().getSamplingRate()
        self._insertFunctionStep(self._convertInputStep)
        self._insertFunctionStep(self._processAllMoviesStep)

    def _convertInputStep(self):
        """ Create settings file. """
        argsDict = {
            "--folder_data": "tmp",
            "--extension": "'*.tiff'",
            "--folder_processing": ".",
            "--angpix": 0.885,
            "--gain_path": "",
            "--defects_path": "",
            "--exposure": -30.648,  # exposure/frame
            "--output": "warp_frameseries.settings",
        }

        cmd = ' '.join(['%s %s' % (k, v) for k, v in argsDict.items()])
        self.runJob(Plugin.getProgram("create_settings"), cmd)

    def _loadInputMovies(self):
        inputMovies = OrderedDict()
        # If there are already some output movies, added them to avoid
        # re-processing them
        outputMovies = getattr(self, 'outputMovies', None)
        if outputMovies is not None:
            self._firstTimeOutput = False
            self.error(f"Existing output: {outputMovies.getSize()} movies")
            for m in outputMovies:
                inputMovies[m.getObjId()] = m.clone()
        else:
            self.error("No output movies ")

        while True:
            moviesFile = self.getInputMovies().getFileName()
            movieSet = SetOfMovies(filename=moviesFile)
            movieSet.loadAllProperties()
            for m in movieSet:
                mid = m.getObjId()
                if mid not in inputMovies:
                    inputMovies[mid] = newMovie = m.clone()
                    yield newMovie
            movieSet.close()
            if movieSet.isStreamClosed():
                break

            time.sleep(self.streamingSleepOnWait.get())
        self.info(f"No more movies, stream closed. Total: {len(inputMovies)}")

    def _processAllMoviesStep(self):
        self.lock = threading.Lock()
        self.processed = 0
        self.registered = 0

        self.error(">>> Starting processing movies...")
        self.program = Plugin.getProgram("fs_motion")
        self.command = self._getCmd()
        self._firstTimeOutput = True

        mc = Pipeline()
        g = mc.addGenerator(self._generateBatches)
        gpus = self.getGpuList()
        mc1 = mc.addProcessor(g.outputQueue, self._getMcProcessor(gpus[0]))
        for gpu in gpus[1:]:
            mc.addProcessor(g.outputQueue, self._getMcProcessor(gpu),
                            outputQueue=mc1.outputQueue)

        o1 = mc.addProcessor(mc1.outputQueue, self._moveBatchOutput)
        mc.run()
        # Mark the output as closed
        self._updateOutputSets([], pwobj.Set.STREAM_CLOSED)

    def _generateBatches(self):
        """ Check for new input movies and generate new tasks. """
        # FIXME: Take streaming into account
        self._batchCount = 0

        def _createBatch(movies):
            batch_id = str(uuid4())
            batch_path = self._getTmpPath(batch_id)
            pwutils.cleanPath(batch_path)
            os.mkdir(batch_path)

            for movie in movies:
                fn = movie.getFileName()
                baseName = os.path.basename(fn)
                os.symlink(os.path.abspath(fn),
                           os.path.join(batch_path, baseName))
            self._batchCount += 1
            return {
                'movies': movies,
                'id': batch_id,
                'path': batch_path,
                'index': self._batchCount
            }

        batchSize = self.streamingBatchSize.get()
        movies = []
        for movie in self._loadInputMovies():
            movies.append(movie)

            if len(movies) == batchSize:
                yield _createBatch(movies)
                movies = []

        if movies:
            yield _createBatch(movies)

    def _getMcProcessor(self, gpu):
        def _processBatch(batch):
            try:
                n = len(batch['movies'])
                t = Timer()

                cmd = self.command.replace('-Gpu #', f'-Gpu {gpu}')
                self.runJob(self.program, cmd, cwd=batch['path'])

                elapsed = t.getToc()
                t.toc(f'Ran fs_motion batch of {n} movies')

                with self.lock:
                    self.processed += n
                    batch_ids = [m.getObjId() for m in batch['movies']]
                    self.error(f"Batch {batch['index']}:{batch_ids}, "
                               f"{elapsed}, Processed: {self.processed}")

                return batch

            except Exception as e:
                self.error("ERROR: fs_motion has failed for batch %s. --> %s\n"
                           % (batch['id'], str(e)))
                import traceback
                traceback.print_exc()

        return _processBatch

    def _moveBatchOutput(self, batch):
        t = Timer()
        srcDir = batch['path']
        doClean = not pwutils.envVarOn(SCIPION_DEBUG_NOCLEAN)
        applyDose = self.doApplyDoseFilter
        saveUnweighted = self._doSaveUnweightedMic()
        usePatches = self.patchX != 0 or self.patchY != 0
        logSuffix = '%s-Full.log' % ('-Patch' if usePatches else '')
        newDone = []

        def _moveToExtra(src, dst):
            srcFn = os.path.join(srcDir, src)
            dstFn = self._getExtraPath(dst)
            print(f"Moving {srcFn} -> {dstFn}\n\texists source: {os.path.exists(srcFn)}")
            if os.path.exists(srcFn):
                pwutils.moveFile(srcFn, dstFn)
                return True
            return False

        def _moveMovieFiles(movie):
            movieRoot = 'output_' + ProtMotionCorr._getMovieRoot(self, movie)
            if applyDose:
                _moveToExtra(movieRoot + '_DW.mrc', self._getOutputMicWtName(movie))

            if not applyDose or saveUnweighted:
                _moveToExtra(movieRoot + '.mrc', self._getOutputMicName(movie))

            if self.splitEvenOdd:
                _moveToExtra(movieRoot + '_EVN.mrc',
                             self._getOutputMicEvenName(movie))
                _moveToExtra(movieRoot + '_ODD.mrc',
                             self._getOutputMicOddName(movie))

            done = _moveToExtra(ProtMotionCorr._getMovieRoot(self, movie) + logSuffix,
                                self._getMovieLogFile(movie))
            if done:
                newDone.append(movie)

        for movie in batch['movies']:
            _moveMovieFiles(movie)

        if newDone:
            self._updateOutputSets(newDone, pwobj.Set.STREAM_OPEN)

        elapsed = t.getToc()
        t.toc('Registered outputs')

        with self.lock:
            self.registered += len(newDone)
            batch_ids = [m.getObjId() for m in batch['movies']]
            self.error(f"OUTPUT: Batch {batch['index']}:{batch_ids}, "
                       f"{elapsed}, "
                       f"New done {len(newDone)}, "
                       f"Registered {self.registered}, "
                       f"Processed {self.processed}")

        self._firstTimeOutput = False

        # Clean batch folder if not in debug mode
        if doClean:
            os.system('rm -rf %s' % batch['path'])

        t.toc(f"Moved output for batch {batch['id']}")

        return batch

    # --------------------------- INFO functions ------------------------------

    # --------------------------- UTILS functions -----------------------------
    def _getCmd(self):
        """ Set return a command string that will be used for each batch. """
        inputMovies = self.getInputMovies()
        argsDict = self._getMcArgs()
        argsDict['-Gpu'] = '#'
        argsDict['-Serial'] = 1


        if inputMovies.getGain():
            argsDict.update({'-Gain': f'"{inputMovies.getGain()}"',
                             '-RotGain': self.gainRot.get(),
                             '-FlipGain': self.gainFlip.get()})

        if inputMovies.getDark():
            argsDict['-Dark'] = inputMovies.getDark()

        # Get input format, but for the batch
        firstMovie = inputMovies.getFirstItem()
        ext = pwutils.getExt(firstMovie.getFileName()).lower()
        if ext in ['.mrc', '.mrcs']:
            inprefix = '-InMrc'
        elif ext in ['.tif', '.tiff']:
            inprefix = '-InTiff'
        else:
            raise Exception(f"Unsupported format '{ext}' for batch processing "
                            f"in Motioncor protocol. ")

        argsDict[inprefix] = './'
        argsDict['-OutMrc'] = 'output_'

        cmd = ' '.join(['%s %s' % (k, v) for k, v in argsDict.items()])
        cmd += self.extraParams2.get()

        return cmd

    def _setPlotInfo(self, movie, mic):
        # FIXME: For now not support PSD or Thumbnail
        if self.doApplyDoseFilter:
            total, early, late = self.calcFrameMotion(movie)
            mic._rlnAccumMotionTotal = Float(total)
            mic._rlnAccumMotionEarly = Float(early)
            mic._rlnAccumMotionLate = Float(late)

    def _getMovieRoot(self, movie):
        return "mic_%06d" % movie.getObjId()

    def _getOutputMovieName(self, movie):
        """ Returns the name of the output movie.
        (relative to micFolder)
        """
        return "movie_%06d" % movie.getObjId()

    def _getOutputMicName(self, movie):
        """ Returns the name of the output micrograph
        (relative to micFolder)
        """
        return self._getMovieRoot(movie) + '.mrc'

    def _getOutputMicWtName(self, movie):
        """ Returns the name of the output dose-weighted micrograph
        (relative to micFolder)
        """
        return self._getMovieRoot(movie) + '_DW.mrc'

    def _getOutputMicEvenName(self, movie):
        """ Returns the name of the output EVEN micrograph
        (relative to micFolder)
        """
        return self._getMovieRoot(movie) + '_EVN.mrc'

    def _getOutputMicOddName(self, movie):
        """ Returns the name of the output EVEN micrograph
        (relative to micFolder)
        """
        return self._getMovieRoot(movie) + '_ODD.mrc'

    def _getOutputMicThumbnail(self, movie):
        return self._getExtraPath(self._getMovieRoot(movie) + '_thumbnail.png')

    def _getMovieLogFile(self, movie):
        usePatches = self.patchX != 0 or self.patchY != 0
        return '%s%s-Full.log' % (self._getMovieRoot(movie),
                                  '-Patch' if usePatches else '')

    def getInputMovies(self):
        return self.inputMovies.get()
