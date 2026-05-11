# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *              Yunior C. Fonseca Reyna (cfonseca@cnb.csic.es)
# *
# * MRC Laboratory of Molecular Biology (MRC-LMB)
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
# **************************************************************************

import os
import mrcfile
import numpy as np

from pyworkflow.constants import ID_ATTRIBUTE
from pyworkflow.protocol import ProtStreamingBase
from pyworkflow.protocol.params import (FloatParam, Positive, LEVEL_ADVANCED, PointerParam, IntParam)
import pyworkflow.utils as pwutils
from pwem.objects import SetOfMicrographs, SetOfMovies, Micrograph
from pwem.protocols import EMProtocol
from pwem.emlib.image import ImageHandler
from pwem.emlib.image.image_readers import ImageStack, ImageReadersRegistry
from tomo.objects import SetOfTiltSeries, SetOfCTFTomoSeries, SetOfTiltSeriesM
from warp import Plugin, WARP_TOOLS, MCORE
from warp.constants import (CREATE_SETTINGS, TOMOSTAR_FOLDER, TILTIMAGES_FOLDER,
                            AVERAGE_FOLDER, TILTSERIES_FOLDER, TILTSERIE_SETTINGS,
                            SETTINGS_FOLDER, WARP_TOOLS_GPU_ALGORITHMS, RELION_FOLDER)
from warp.utils import tom_deconv, tomoStarGenerate


class ProtWarpBase(EMProtocol):
    """
    Base infrastructure for Warp-related preprocessing, deconvolution,
    and streaming alignment workflows for both movies and tilt-series.

    AI Generated:

    Warp Processing and Streaming Alignment (ProtWarpBase) — User Manual
        Overview

        This protocol family provides a common computational framework for
        preprocessing cryo-EM data before downstream reconstruction,
        refinement, or tomographic analysis.

        The provided code combines three major functional areas:

        1. Deconvolution of micrographs and tilt-series.
        2. Streaming processing of movie datasets.
        3. Streaming processing of tilt-series movie datasets.

        Although the original source defines several classes, their behavior
        can be understood as one unified processing framework whose goal is
        to transform incoming experimental data into corrected, analysis-ready
        outputs while preserving compatibility with Scipion streaming logic.

        General Purpose

        In cryo-EM workflows, raw images often contain strong contrast-transfer
        effects, low-frequency background, and acquisition artifacts that
        complicate interpretation.

        This protocol addresses those limitations by providing:

        - Wiener-like deconvolution using CTF information.
        - GPU-aware execution of Warp tools.
        - Preparation of tomography-specific metadata.
        - Streaming ingestion of incoming movies or tilt-series.
        - Incremental registration of processed outputs.

        This design is particularly useful in facility pipelines, automated
        acquisition environments, or high-throughput cryo-EM workflows.

        --------------------------------------------------------------------
        DECONVOLUTION WORKFLOW
        --------------------------------------------------------------------

        Input Requirements

        The deconvolution workflow expects:

        - Input images or stacks.
        - Acquisition metadata.
        - CTF information indexed by image identity.

        Each input image is matched against a CTF dictionary. Only inputs
        with valid defocus information are processed.

        Main Deconvolution Parameters

        deconvstrength
            Controls the strength of the deconvolution filter.

        snrfalloff
            Controls how rapidly signal-to-noise decreases with frequency.

        highpassnyquist
            Defines the low-frequency cutoff expressed as a fraction of
            Nyquist frequency.

        Practical Interpretation

        These parameters regulate how strongly low-frequency background
        is suppressed and how much high-frequency information is restored.

        Conservative values are recommended initially. Excessively strong
        deconvolution can amplify noise and introduce artifacts.

        Processing Logic

        For every input item:

        - Its metadata key is matched against the CTF dictionary.
        - Defocus is extracted and converted.
        - Output filename is generated.
        - The image is processed either as:

            * a single micrograph
            * a tilt-series stack

        Image Processing Modes

        _processImage
            Reads one image, applies tomographic deconvolution, and writes
            a corrected MRC file.

        _processStack
            Reads all slices from an image stack, processes them one by one,
            and writes a corrected MRC stack.

        Output Naming

        Output files are stored in the protocol extra directory and receive
        the suffix:

            _deconv.mrc

        --------------------------------------------------------------------
        GPU EXECUTION
        --------------------------------------------------------------------

        GPU-aware command execution is handled automatically.

        If GPUs are available:

        - Warp GPU algorithms receive a device list.
        - MCore-compatible programs receive Warp-specific GPU flags.

        The framework also warns when thread counts are incompatible with
        GPU execution.

        Practical Note

        In GPU mode, CPU thread count may be ignored depending on the
        underlying Warp implementation.

        --------------------------------------------------------------------
        TOMOGRAPHY PREPARATION
        --------------------------------------------------------------------

        Tilt-Series Data Preparation

        For tomographic workflows, the protocol can prepare per-tilt metadata.

        For each tilt image:

        - Individual images are extracted from the tilt-series.
        - Acquisition metadata are collected.
        - Alignment shifts are converted to physical units.
        - Tilt angles, dose, shifts, and contrast metadata are stored.

        A tomostar file is then generated.

        This preparation allows Warp to understand both geometry and
        acquisition conditions of the tilt-series.

        Auxiliary Folder Structure

        The preparation stage creates:

        - tilt-image folders
        - tomostar folders
        - symbolic links required by Warp

        This organization is necessary because Warp expects both original
        image headers and average-accessible image paths.

        Tilt-Series Settings Generation

        The protocol can also create Warp tilt-series settings files.

        These settings may include:

        - sampling rate
        - exposure per frame
        - tomogram dimensions
        - EER grouping information

        These settings are required before downstream Warp tomography
        programs can be executed.

        --------------------------------------------------------------------
        STREAMING MOVIE PROCESSING
        --------------------------------------------------------------------

        Streaming Purpose

        The streaming movie framework is designed for acquisition-time
        processing of incoming movie data.

        Instead of waiting for a full dataset to finish acquisition,
        the protocol continuously monitors the input set.

        Streaming Workflow

        During execution:

        - new movies are detected
        - only unseen inputs are selected
        - items are grouped into batches
        - processing steps are inserted dynamically
        - outputs are appended incrementally

        Batch Size Behavior

        streamingBatchSize = 1
            Process one movie at a time.

        streamingBatchSize > 1
            Group several movies into one processing step.

        streamingBatchSize = 0
            Process all currently available items together.

        Output Registration

        Processed outputs are converted into micrographs.

        For each generated micrograph:

        - input acquisition metadata are copied
        - accumulated dose is updated
        - sampling rate is propagated
        - the item is appended to the streaming output set

        Streaming State Management

        Two internal dictionaries keep track of state:

        _moviesToProcess
            Newly detected inputs waiting for processing.

        _moviesInProcess
            Inputs already dispatched for execution.

        This prevents duplicate processing.

        --------------------------------------------------------------------
        STREAMING TILT-SERIES PROCESSING
        --------------------------------------------------------------------

        The tilt-series streaming framework follows the same philosophy,
        but the processing unit becomes an entire tilt-series rather than
        a single movie.

        Workflow

        - new tilt-series are detected
        - already processed ones are ignored
        - each tilt-series launches its own processing step
        - optional finalization steps may be inserted
        - output sets are closed once streaming ends

        Output Types

        The framework supports tilt-series outputs and optionally
        CTF-tomography outputs.

        Compared with movie streaming, the logic is simpler because
        each streaming task corresponds to one tilt-series identifier.

        --------------------------------------------------------------------
        PATH NORMALIZATION UTILITIES
        --------------------------------------------------------------------

        Several helper methods normalize paths for downstream tools.

        They are used to locate:

        - particle files
        - tomograms
        - optimization files

        If relative paths are found, they are resolved against protocol
        output folders.

        This is particularly useful in iterative Warp/Relion workflows.

        --------------------------------------------------------------------
        EXTENSION POINTS FOR DEVELOPERS
        --------------------------------------------------------------------

        The framework intentionally leaves some methods unimplemented.

        These must be specialized by derived protocols.

        Main abstract extension points include:

        deconvolveStep()
            Defines how deconvolution is triggered.

        createOutputStep()
            Defines how outputs are created.

        proccessMoviesStep()
            Defines batch processing logic for streaming movies.

        proccessTSMoviesStep()
            Defines processing logic for streaming tilt-series.

        insertInitialSteps()
            Allows protocols to insert pre-processing steps.

        insertFinalSteps()
            Allows protocols to insert finalization steps.

        --------------------------------------------------------------------
        PRACTICAL USAGE RECOMMENDATIONS
        --------------------------------------------------------------------

        For standard micrograph deconvolution:

        - begin with conservative deconvolution strength
        - verify output visually before large-scale processing

        For streaming workflows:

        - use batch size 1 for immediate feedback
        - use larger batch sizes when minimizing process-launch overhead

        For tomography:

        - verify tilt metadata carefully
        - ensure acquisition geometry is correct before reconstruction

        For GPU execution:

        - prefer GPU when available
        - avoid relying on thread count for performance scaling

        --------------------------------------------------------------------
        FINAL PERSPECTIVE
        --------------------------------------------------------------------

        This protocol base is not a single algorithm but a shared execution
        architecture for Warp-oriented preprocessing.

        Its main role is to make cryo-EM preprocessing scalable,
        automation-friendly, and compatible with modern streaming acquisition.

        In practice, it provides the infrastructure needed to move from
        raw experimental movies toward reproducible, structured,
        downstream-ready cryo-EM data products.
    """
    _label = None

    # -------------------------- DEFINE param functions -----------------------
    @classmethod
    def defineProcessParams(cls, form):
        form.addParam('deconvstrength', FloatParam, validators=[Positive],
                      default=1.0,
                      expertLevel=LEVEL_ADVANCED,
                      label='Deconvolution strength',
                      help='Strength parameter for the deconvolution filter.')
        form.addParam('snrfalloff', FloatParam, validators=[Positive],
                      default=1.1,
                      expertLevel=LEVEL_ADVANCED,
                      label='SNR falloff',
                      help='SNR falloff parameter for the deconvolution filter.')
        form.addParam('highpassnyquist', FloatParam, validators=[Positive],
                      default=0.02,
                      expertLevel=LEVEL_ADVANCED,
                      label='High-pass fraction',
                      help='Fraction of Nyquist frequency to be cut off on '
                           'the lower end (since it will be boosted the most).')
        
    # --------------------------- STEPS functions -----------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.deconvolveStep)
        self._insertFunctionStep(self.createOutputStep)

    def deconvolveStep(self):
        raise NotImplementedError

    def createOutputStep(self):
        raise NotImplementedError

    def _deconvolve(self, pixSize, acquisition, inputList,
                    ctfDict, keyName="_micName", isTS=False):
        """ Main function to iterate over inputList items.
        Checks the match with ctfDict by keyName for each item.
        :param pixSize: input pixel size
        :param acquisition: input acquisition object
        :param inputList: list of dicts
        :param ctfDict: defocus dict
        :param keyName: ket for inputList
        :param isTS: flag to process TS mrcs stack
        .
        """
        pwutils.cleanPath(self._getExtraPath())
        pwutils.makePath(self._getExtraPath())

        voltage = acquisition.getVoltage()
        cs = acquisition.getSphericalAberration()
        ncpu = self.numberOfThreads.get()
        gpu = self.usesGpu()
        gpuid = self.getGpuList()[0]

        snrfalloff = self.snrfalloff.get()
        deconvstrength = self.deconvstrength.get()
        highpassnyquist = self.highpassnyquist.get()

        for item in inputList:
            key = item[keyName]
            if key in ctfDict:
                fileName = item["_filename"]
                defocus = ctfDict[key] / 10000
                outputFn = self._getOutputFn(fileName)
                self.info(f"Deconvolving {fileName}")

                func = self._processStack if isTS else self._processImage
                func(fileName, outputFn, angpix=pixSize, voltage=voltage,
                     cs=cs, defocus=defocus, ncpu=ncpu, gpu=gpu, gpuid=gpuid,
                     snrfalloff=snrfalloff, deconvstrength=deconvstrength,
                     highpassnyquist=highpassnyquist)
            else:
                self.warning(f"No CTF found for: {key}")

    # --------------------------- INFO functions ------------------------------
    def _warnings(self):
        warnings = []

        if self.usesGpu() and self.numberOfThreads > 1:
            warnings.append("Number of threads is ignored for GPU execution")

        return warnings

    # -------------------------- UTILS functions ------------------------------
    def _getOutputFn(self, micName):
        return self._getExtraPath(pwutils.removeBaseExt(micName) + "_deconv.mrc")

    def _updateItem(self, item, row):
        outputFn = self._getOutputFn(item.getFileName())
        if os.path.exists(outputFn):
            item.setFileName(outputFn)
        else:
            item._appendItem = False

    @staticmethod
    def _processImage(inputFn, outputFn, **kwargs):
        ih = ImageHandler()
        inputData = ih.read(inputFn).getData()
        result = tom_deconv(inputData, **kwargs)
        with mrcfile.new(outputFn) as mrcOut:
            mrcOut.set_data(result)
            mrcOut.voxel_size = kwargs["angpix"]
            mrcOut.update_header_from_data()

    @staticmethod
    def _processStack(inputFn, outputFn, **kwargs):
        ih = ImageHandler()
        x, y, z, n = ih.getDimensions(inputFn)
        nImages = max(z, n)  # Just handle ambiguity with mrc format
        stack_shape = (nImages, y, x)
        mrc = mrcfile.new_mmap(outputFn, shape=stack_shape,
                               mrc_mode=2, overwrite=True)
        for i in range(nImages):
            inputData = ih.read((i + 1, inputFn)).getData()
            mrc.data[i] = tom_deconv(inputData, **kwargs)
        mrc.reset_header_stats()
        mrc.update_header_from_data()
        mrc.set_image_stack()
        mrc.voxel_size = kwargs["angpix"]
        mrc.close()

    @staticmethod
    def _createDataImportSettings():
        """ Create data import settings"""
        pass

    def runProgram(self, argsDict, program, algorithm, othersCmds=None):
        gpuList = self.getGpuList()
        if gpuList:
            deviceList = ' '.join(map(str, gpuList))
            if algorithm in WARP_TOOLS_GPU_ALGORITHMS:
                argsDict['--device_list'] = deviceList
            elif program == MCORE:
                argsDict['--devicelist'] = deviceList

        cmd = ' '.join(['%s %s' % (k, v) for k, v in argsDict.items()])
        if othersCmds:
            cmd += ' %s' % othersCmds
        self.runJob(self.getPlugin().getProgram(program, algorithm), cmd, executable='/bin/bash')

    def tsDataPrepare(self, ts, perTs=False):
        """Creates the setting file that will be used by the different programs.
           It also extracts the tiltimages from the tiltseries and generates the *.tomostar files based on
           the tiltimages."""
        tsId = ts.getTsId()
        starFolder = self._getExtraPath(TOMOSTAR_FOLDER)
        pwutils.makePath(starFolder)
        if perTs:
            starFolder = self._getExtraPath(TOMOSTAR_FOLDER, tsId)
            pwutils.makePath(starFolder)
        objSet = self.inputSet.get()
        imagesFolder = self._getExtraPath(TILTIMAGES_FOLDER)
        invertTiltAngle = 1
        pwutils.makePath(imagesFolder)
        hasAlignment = objSet.hasAlignment()
        sr = objSet.getSamplingRate()
        fileDict = []

        if ts.isEnabled():
            tsId = ts.getTsId()
            properties = {"sr": sr}
            tiValues = {}
            for ti in ts.iterItems():
                if ti.isEnabled():  # Excluding views
                    dose = 0
                    maskedFraction = 0
                    shiftX = 0
                    shiftY = 0
                    axisAngle = 0
                    amplitudeContrast = 0

                    if ts.hasAcquisition():
                        axisAngle = ts.getAcquisition().getTiltAxisAngle()
                    if ti.getAcquisition():
                        amplitudeContrast = ti.getAcquisition().getAmplitudeContrast()
                        dose = ti.getAcquisition().getAccumDose()
                    newBinaryName = tsId + f'_TO_%02d.mrc' % ti.getAcquisitionOrder()
                    fileDict.append(newBinaryName)
                    newFrame = ImageStack(properties=properties)
                    newFrame.append(ImageReadersRegistry.open(str(ti.getIndex()) + '@' + ti.getFileName()))
                    ImageReadersRegistry.write(newFrame, os.path.join(self._getExtraPath(TILTIMAGES_FOLDER), newBinaryName))
                    # Taking angleTilt axisAngle Dose AverageIntensity MaskedFraction
                    if hasAlignment:
                        transform = ti.getTransform()
                        matrix = transform.getMatrix()
                        newMatrix = np.zeros((3, 3), dtype=float)
                        newMatrix[0, 0:2] = matrix[0, 0:2]
                        newMatrix[1, 0:2] = matrix[1, 0:2]
                        newMatrix[2, 2] = 1
                        tiShift = [-1 * matrix[0, 2], -1 * matrix[1, 2], 0]
                        transpose = np.transpose(newMatrix)
                        multShift = np.dot(transpose, tiShift)
                        shiftX = multShift[0] * sr
                        shiftY = multShift[1] * sr

                    tiValues[ti.getTiltAngle() * invertTiltAngle] = [newBinaryName,
                                                                     ti.getTiltAngle() * invertTiltAngle,
                                                                     axisAngle, shiftX, shiftY, dose,
                                                                     amplitudeContrast, maskedFraction]

            tomoStarGenerate(tsId, tiValues, starFolder, True, perTs=perTs)

        # 2. Create symbolic links from tiltseries folder to average folder
        #    We need to do this because warp needs both folders: The tiltimages folder to get
        #    the header of the images (it only needs the first tiltimage of each tiltseries),
        #    and the averages folder to read them.
        averagesFolder = os.path.join(imagesFolder, AVERAGE_FOLDER)
        pwutils.makePath(averagesFolder)
        for file in fileDict:
            fileName = os.path.basename(file)
            destFolder = os.path.join(averagesFolder, fileName)
            if not os.path.exists(destFolder):
                os.symlink('../' + fileName, destFolder)

    def createTiltSeriesSetting(self, ts, perTs=False):
        tsId = ts.getTsId() if ts is not None else ''
        extension = "*.tomostar" if perTs or ts is None else "%s.tomostar" % tsId
        tomoStarFolder = os.path.abspath(self._getExtraPath(TOMOSTAR_FOLDER))
        folderPath = tomoStarFolder if not perTs else os.path.join(tomoStarFolder, tsId)
        self.info(">>> Starting tilt-series settings creation (%s)..." % tsId)
        setOfTS = self.inputSet.get()
        sr = setOfTS.getSamplingRate()
        exposure = setOfTS.getAcquisition().getDosePerFrame()
        settingsFolder = os.path.abspath(self._getExtraPath())
        if ts is not None:
            settingsFolder = os.path.abspath(self._getExtraPath(SETTINGS_FOLDER))
            pwutils.makePath(settingsFolder)
        processingFolder = os.path.abspath(self._getExtraPath(TILTSERIES_FOLDER))
        pwutils.makePath(processingFolder)
        tsSettingFile = tsId + '_' + TILTSERIE_SETTINGS if ts is not None else TILTSERIE_SETTINGS
        tsSettingFilePath = os.path.abspath(os.path.join(self._getExtraPath(settingsFolder), tsSettingFile))
        argsDict = {
            "--folder_data": folderPath,
            "--extension": extension,
            "--folder_processing": processingFolder,
            '--angpix': sr,
            "--output": tsSettingFilePath
        }

        if exposure is not None:
            argsDict['--exposure'] = exposure

        if hasattr(self, 'tomo_thickness'):
            z = self.tomo_thickness.get()
            x = self.x_dimension.get() or setOfTS.getDimensions()[0]
            y = self.y_dimension.get() or setOfTS.getDimensions()[1]

            argsDict['--tomo_dimensions'] = f'{x}x{y}x{z}'

        if extension == '.eer':
            argsDict['--eer_ngroups'] = self.eer_ngroups.get()
            if self.eer_groupexposure.get():
                argsDict['--eer_groupexposure'] = self.eer_groupexposure.get()

        cmd = ' '.join(['%s %s' % (k, v) for k, v in argsDict.items()])

        self.runJob(Plugin.getProgram(WARP_TOOLS, CREATE_SETTINGS), cmd, executable='/bin/bash')

    def normalizeParticlesPath(self, particleValue):
        particleFile = particleValue
        if not os.path.exists(particleFile):
            particleFile = os.path.normpath(os.path.join(self._getExtraPath(TILTSERIES_FOLDER), particleValue))
        return particleFile

    def normalizeTomogramsPath(self, tomogramValue):
        tomogramFile = tomogramValue
        if not os.path.exists(tomogramFile):
            tomogramFile = os.path.join(self._getExtraPath(RELION_FOLDER), tomogramValue)
        return tomogramFile

    def normalizeOptimizationPath(self, optimizationValue):
        optimizationFile = optimizationValue
        if not os.path.exists(optimizationFile):
            optimizationFile = os.path.join(self._getExtraPath(RELION_FOLDER), optimizationValue)
        return optimizationFile


class ProtMovieAlignBase(EMProtocol, ProtStreamingBase):
    """
    Protocol base to use in streaming motion correction protocols
    """
    OUT_MICS = 'Micrographs'
    _possibleOutputs = {OUT_MICS: SetOfMicrographs}
    _lastInputId = 0
    MIC_NAME_ATTR = '_micName'
    IMAGE_FILENAME_ATTR = '_filename'

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        ProtStreamingBase.__init__(self, **kwargs)
        self._moviesToProcess = {}
        self._moviesInProcess = {}

    # --------------------------- DEFINE param functions ----------------------
    def _defineInputMoviesParam(self, form):
        """ Defines the 'inputMovies' parameter for the provided protocol form.
            This method adds a parameter to the protocol, allowing the user to select
            a set of previously imported movies. """
        form.addParam('inputMovies', PointerParam, pointerClass=SetOfMovies,
                      important=True,
                      label=pwutils.Message.LABEL_INPUT_MOVS,
                      help='Select a set of previously imported movies.')

    def _defineStreamingParams(self, form):
        """ This function can be called during the _defineParams method
        of some protocols that support stream processing.
        It will add a Streaming section together with the following
        params:
            streamingBatchSize: For some programs it is more efficient to process
                many items at once and not one by one. So this parameter will
                allow to group a number of items to be processed in the same
                protocol step. This can also reduce some IO overhead and spawning
                new OS processes.
        """
        ProtStreamingBase._defineStreamingParams(self, form)
        form.getSection("Streaming")
        form.addParam("streamingBatchSize", IntParam, default=1,
                      label="Batch size",
                      help="This value allows to group several items to be "
                           "processed inside the same protocol step. You can "
                           "use the following values: \n"
                           "*1*    The default behavior, the items will be "
                           "processed one by one.\n"
                           "*0*    Put in the same step all the items "
                           "available. If the sleep time is short, it could be "
                           "practically the same of one by one. If not, you "
                           "could have steps with more items. If the steps will "
                           "be executed in parallel, it is better not to use "
                           "this option.\n"
                           "*>1*   The number of items that will be grouped into "
                           "a step.")

    def stepsGeneratorStep(self) -> None:
        """
        Generates and inserts processing steps for streaming input movies, processing them in batches.
        This step should be implemented by any streaming protocol.
        It should check its input and when ready conditions are met, call the self._insertFunctionStep method.
        The method repeatedly checks for new input movies, processes them in batches, and ensures
        that all steps are completed before closing the output.
        """
        self.getNewInputMicNames(removeDone=True)
        batches = self._getStreamingBatchSize()  # The number of movies processed in each batch
        alive = True  # indicate whether the streaming input is still active.

        initialSteps = self.insertInitialSteps()  # the initial steps inserted into the protocol(must be re-implemented in the plugin).
        proccessMoviesSteps = []  # list that stores the processing steps for each batch of movies
        while alive:
            if self.getInputMovies(loadProperties=True).isStreamClosed():
                alive = False
            # split the list of micNames to be processed and divide them according to the value of the batch size
            micNamesToProcess = list(self._moviesToProcess.keys())
            splitMicNameList = [micNamesToProcess[i:i + batches] for i in range(0, len(micNamesToProcess), batches)]
            # move into the list and launch steps with micNames sublist.
            # we assume that the last sublist may be incomplete(remnant).
            for micNameList in splitMicNameList:
                if len(micNameList) == batches or not alive:
                    self.info('List of movies to be precessed: %s' % micNameList)
                    self.removeDoneMics(micNameList)
                    proccessMoviesStep = self._insertFunctionStep(self.proccessMoviesStep, micNameList,
                                                                  prerequisites=initialSteps, needsGPU=self.usesGpu())
                    proccessMoviesSteps.append(proccessMoviesStep)

            if not alive:
                break
            self._streamingSleepOnWait()
            self.getNewInputMicNames()
        self._insertFunctionStep(self.closeOutputStep, prerequisites=proccessMoviesSteps, needsGPU=False)

    def insertInitialSteps(self) -> list:
        """The initial steps inserted into the protocol"""
        return []

    def proccessMoviesStep(self, micNamesList) -> None:
        pass

    def addMicrographs(self, micNameList, micLocations) -> None:
        """ Adds micrographs to the output set based on the provided list of names and locations.
            Parameters:
            - micNameList (list): List of micrograph names to get the associated movie
            - micLocations (list): List of locations corresponding to each micrograph.
            """
        with self._lock:
            output = self.getOutputSet()
            for index, micName in enumerate(micNameList):
                movie = self.getInputMovies().getItem(self.MIC_NAME_ATTR, micName)
                newMic = Micrograph(location=micLocations[index])
                newMic.copyInfo(movie)
                output.append(newMic)
                movieAcq = movie.getAcquisition()
                newMic.getAcquisition().setDosePerFrame(movieAcq.getDosePerFrame() * movie.getNumberOfFrames() +
                                                        movieAcq.getDoseInitial())
                newMic.getAcquisition().setDoseInitial(0)
                output.update(newMic)
                output.write()
                self._store(output)

    def removeDoneMics(self, micNameList) -> None:
        """Moves the micNames from `_moviesToProcess`to `_moviesInProcess`,
           indicating that it has been processed."""
        for micName in micNameList:
            self._moviesInProcess[micName] = self._moviesToProcess[micName]
            self._moviesToProcess.pop(micName)

    def getOutputMicNames(self):
        """ Retrieves the unique names of the micrographs from the output set.
            Returns:
            - list: A list of unique micrograph names, or an empty list if the output set is not available.
            """
        outputSet = self.getOutputSet(createSet=False)
        if outputSet is not None:
            return outputSet.getUniqueValues(self.MIC_NAME_ATTR)  #TODO When pwem is released, use Micrograph.MIC_NAME
        return []

    def getNewInputMicNames(self, removeDone=False):
        """ Retrieves new micrograph names from the input set, with an option to exclude already processed ones.
            This method fetches unique micrograph names, IDs, and file paths from the input movie set,
            where the IDs are greater than the last processed input ID (`self._lastInputId`).

            Parameters:
            - removeDone (bool): if True, excludes micrographs that have already been processed (present in the output set).
            """
        inputSet = self.getInputMovies()
        results = inputSet.getUniqueValues([self.MIC_NAME_ATTR, ID_ATTRIBUTE, self.IMAGE_FILENAME_ATTR],
                                           '%s > %s' % (ID_ATTRIBUTE, self._lastInputId))
        newIds = results[ID_ATTRIBUTE]
        if newIds:
            self._lastInputId = max(results[ID_ATTRIBUTE])
            outputMicNames = self.getOutputMicNames() if removeDone else []
            for index, micName in enumerate(results[self.MIC_NAME_ATTR]):
                if not removeDone or micName not in outputMicNames:
                    self._moviesToProcess[micName] = results[self.IMAGE_FILENAME_ATTR][index]

    def closeOutputStep(self):
        """Finalizes and closes the output sets. """
        self._closeOutputSet()

    def getInputMovies(self, loadProperties=False) -> SetOfMovies:
        """ Retrieves the set of input movies.
            Returns:
            - SetOfMovies: The set of input movies to be processed.  """
        if loadProperties:
            self.inputMovies.get().loadAllProperties()
        return self.inputMovies.get()

    def getOutputSet(self, createSet=True) -> SetOfMicrographs:
        """ Retrieves the output set of micrographs, creating it if necessary.
            Parameters:
            - createSet (bool): If True, a new output set will be created if it doesn't already exist.
            Returns:
            - SetOfMicrographs: The output set of micrographs, either existing or newly created.
            """
        outputSet = getattr(self, self.OUT_MICS, None)
        if createSet and outputSet is None:
            outputSet = SetOfMicrographs.create(self._getPath())
            outputSet.setStreamState(outputSet.STREAM_OPEN)
            inputMovies = self.getInputMovies()
            outputSet.copyInfo(inputMovies)
            outputSet.setSamplingRate(inputMovies.getSamplingRate() * self.binFactor.get())

            self._defineOutputs(**{self.OUT_MICS: outputSet})
            self._defineSourceRelation(self.getInputMovies(), outputSet)

        return outputSet

    def _getStreamingBatchSize(self):
        """
            Retrieves the configured batch size for streaming operations.
            Returns:
            - int: The batch size for streaming operations.
            """
        return self.getAttributeValue('streamingBatchSize', 1)

    def getFileName(self, micName):
        """ Retrieves the file name associated with a given micrograph name.
            Parameters:
            - micName (str): The name of the micrograph for which the file name is requested.
            Returns:
            - str or None: The file name associated with the specified micrograph name,
              or None if the `micName` is not found in `_moviesToProcess`.
            """
        if micName in self._moviesInProcess:
            return self._moviesInProcess[micName]
        return None


class ProtTSMovieAlignBase(EMProtocol, ProtStreamingBase):
    """
    Protocol base to use in streaming ts motion correction protocols
    """
    OUT_TS = 'TiltSeries'
    OUT_CTFTOMO = 'CTFTomoSeries'
    _possibleOutputs = {OUT_TS: SetOfTiltSeries,
                        OUT_CTFTOMO: SetOfCTFTomoSeries}
    _lastInputId = 0
    TS_ID_ATTR = '_tsId'

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        ProtStreamingBase.__init__(self, **kwargs)
        self._moviesToProcess = {}
        self._moviesInProcess = {}

    def _defineStreamingParams(self, form):
        """ This function can be called during the _defineParams method
        of some protocols that support stream processing.
        It will add a Streaming section together with the following
        params:
            streamingSleepOnWait: Some streaming protocols are quite fast,
                so, checking input/output updates creates an IO overhead.
                This params allows them to sleep (without consuming resources)
                to wait for new work to be done.
        """
        ProtStreamingBase._defineStreamingParams(self, form)

    def stepsGeneratorStep(self) -> None:
        """
        Generates and inserts processing steps for streaming input tiltseries movies, processing them in batches.
        This step should be implemented by any streaming protocol.
        It should check its input and when ready conditions are met, call the self._insertFunctionStep method.
        The method repeatedly checks for new input movies, processes them in batches, and ensures
        that all steps are completed before closing the output.
        """
        self.getNewInputTSMovies(removeDone=True)
        batches = self._getStreamingBatchSize()  # The number of tilt images movies processed(per tiltseies movies) in each batch
        alive = True  # indicate whether the streaming input is still active.

        initialSteps = self.insertInitialSteps()  # the initial steps inserted into the protocol(must be re-implemented in the plugin).
        proccessTSMoviesSteps = []  # list that stores the processing steps for each batch of movies
        while alive:
            if self.getInputTSMovies(loadProperties=True).isStreamClosed():
                alive = False
            tsToProcess = list(self._moviesToProcess.keys())
            # move into the list and launch steps with tiltseries movies.
            # we assume that the last sublist may be incomplete(remnant).
            for tsId in tsToProcess:
                self.removeDoneTS(tsId)
                proccessTSMoviesStep = self._insertFunctionStep(self.proccessTSMoviesStep, tsId,
                                                                prerequisites=initialSteps, needsGPU=self.usesGpu())
                proccessTSMoviesSteps.append(proccessTSMoviesStep)

            if not alive:
                break
            self._streamingSleepOnWait()
            self.getNewInputTSMovies()
        finalSteps = self.insertFinalSteps(proccessTSMoviesSteps)
        self._insertFunctionStep(self.closeOutputStep, prerequisites=proccessTSMoviesSteps + finalSteps, needsGPU=False)

    def insertInitialSteps(self) -> list:
        """The initial steps inserted into the protocol"""
        return []

    def insertFinalSteps(self, proccessTSMoviesSteps) -> list:
        """The final steps inserted into the protocol"""
        return []

    def proccessTSMoviesStep(self, tsId) -> None:
        pass

    def removeDoneTS(self, tsId) -> None:
        """Moves the micNames from `_moviesToProcess`to `_moviesInProcess`,
           indicating that it has been processed."""
        self._moviesInProcess[tsId] = self._moviesToProcess[tsId]
        self._moviesToProcess.pop(tsId)

    def getOutputTS(self):
        """ Retrieves the unique names of the TS from the output set.
            Returns:
            - list: A list of unique TS ids, or an empty list if the output set is not available.
            """
        outputSet = self.getOutputSet(createSet=False)
        if outputSet is not None:
            return outputSet.getUniqueValues(self.TS_ID_ATTR)  #TODO When pwem is released, use Micrograph.MIC_NAME
        return []

    def getNewInputTSMovies(self, removeDone=False):
        """ Retrieves new tiltseries movies from the input set, with an option to exclude already processed ones.
            This method fetches unique tiltseries names, IDs, and file paths from the input tilt series movie set,
            where the IDs are greater than the last processed input ID (`self._lastInputId`).

            Parameters:
            - removeDone (bool): if True, excludes tilt series movies that have already been processed (present in the output set).
            """
        inputSet = self.getInputTSMovies()
        results = inputSet.getUniqueValues([self.TS_ID_ATTR, ID_ATTRIBUTE],
                                           '%s > %s' % (ID_ATTRIBUTE, self._lastInputId))
        newIds = results[ID_ATTRIBUTE]
        if newIds:
            self._lastInputId = max(results[ID_ATTRIBUTE])
            outputTS = self.getOutputTS() if removeDone else []
            for index, tsId in enumerate(results[self.TS_ID_ATTR]):
                if not removeDone or tsId not in outputTS:
                    self._moviesToProcess[tsId] = results[self.TS_ID_ATTR][index]

    def closeOutputStep(self):
        """Finalizes and closes the output sets. """
        self._closeOutputSet()

    def getInputTSMovies(self, loadProperties=False) -> SetOfTiltSeriesM:
        """ Retrieves the set of input ts movies.
            Returns:
            - SetOfTSMovies: The set of input tilt series movies to be processed.  """
        if loadProperties:
            self.inputTSMovies.get().loadAllProperties()
        return self.inputTSMovies.get()

    def getOutputSet(self, createSet=True) -> SetOfTiltSeries:
        """ Retrieves the output set of TS, creating it if necessary.
            Parameters:
            - createSet (bool): If True, a new output set will be created if it doesn't already exist.
            Returns:
            - SetOfTiltSeries: The output set of TS, either existing or newly created.
            """
        outputSet = getattr(self, self.OUT_TS, None)
        return outputSet


