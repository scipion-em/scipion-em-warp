# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
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
import time
import mrcfile
import glob
import numpy

from pyworkflow.constants import ID_ATTRIBUTE
from pyworkflow.protocol import ProtStreamingBase
from pyworkflow.protocol.params import (FloatParam, Positive, LEVEL_ADVANCED, PointerParam, IntParam)
import pyworkflow.utils as pwutils
from pwem.objects import SetOfMicrographs, SetOfMovies, Micrograph
from pwem.protocols import EMProtocol
from pwem.emlib.image import ImageHandler
from pwem.emlib.image.image_readers import ImageStack, ImageReadersRegistry
from tomo.objects import TiltSeries, SetOfTiltSeries
from warp.constants import (CREATE_SETTINGS, TOMOSTAR_FOLDER, TILTIMAGES_FOLDER,
                            AVERAGE_FOLDER, TILTSERIES_FOLDER, TILTSERIE_SETTINGS,
                            FRAMES_FOLDER, FRAMESERIES_FOLDER)
from warp.utils import tom_deconv


class ProtWarpBase(EMProtocol):
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

    def runProgram(self, argsDict, program, othersCmds=None):
        gpuList = self.getGpuList()
        if gpuList:
            argsDict['--device_list'] = ' '.join(map(str, gpuList))

        cmd = ' '.join(['%s %s' % (k, v) for k, v in argsDict.items()])
        if othersCmds:
            cmd += ' %s' % othersCmds
        self.runJob(self.getPlugin().getProgram(program), cmd, executable='/bin/bash')

    def dataPrepare(self, objSet):
        """Creates the setting file that will be used by the different programs.
           It also extracts the tiltimages from the tiltseries and generates the *.tomostar files based on
           the tiltimages."""
        plugin = self.getPlugin()
        starFolder = self._getExtraPath(TOMOSTAR_FOLDER)
        pwutils.makePath(starFolder)
        isTiltSeries = isinstance(objSet, SetOfTiltSeries)
        imagesFolder = self._getExtraPath(TILTIMAGES_FOLDER) if isTiltSeries else self._getExtraPath(FRAMES_FOLDER)
        invertTiltAngle = -1 if isTiltSeries else 1
        pwutils.makePath(imagesFolder)
        hasAlignment = objSet.hasAlignment()
        sr = objSet.getSamplingRate()
        exposure = objSet.getAcquisition().getDosePerFrame()

        for ts in objSet.iterItems():
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
                        fileName = ti.getFileName()
                        newBinaryName = tsId + f'_TO_%02d.mrc' % ti.getAcquisitionOrder()
                        if isTiltSeries:
                            newFrame = ImageStack(properties=properties)
                            newFrame.append(ImageReadersRegistry.open(str(ti.getIndex()) + '@' + ti.getFileName()))
                            ImageReadersRegistry.write(newFrame, os.path.join(self._getExtraPath(TILTIMAGES_FOLDER), newBinaryName))
                            # Taking angleTilt axisAngle Dose AverageIntensity MaskedFraction
                            if hasAlignment:
                                transform = ti.getTransform()
                                matrix = transform.getMatrix()
                                newMatrix = numpy.zeros((3, 3), dtype=float)
                                newMatrix[0, 0:2] = matrix[0, 0:2]
                                newMatrix[1, 0:2] = matrix[1, 0:2]
                                newMatrix[2, 2] = 1
                                tiShift = [-1 * matrix[0, 2], -1 * matrix[1, 2], 0]
                                transpose = numpy.transpose(newMatrix)
                                multShift = numpy.dot(transpose, tiShift)
                                shiftX = multShift[0] * sr
                                shiftY = multShift[1] * sr
                        else:
                            newBinaryName = os.path.basename(fileName)
                            os.symlink(os.path.abspath(fileName), os.path.join(imagesFolder, os.path.basename(fileName)))

                        tiValues[ti.getTiltAngle() * invertTiltAngle] = [newBinaryName, ti.getTiltAngle() * invertTiltAngle,
                                                                         axisAngle, shiftX, shiftY, dose,
                                                                         amplitudeContrast, maskedFraction]

                self.tomoStarGenerate(tsId, tiValues, starFolder, isTiltSeries)

        # 2. Create symbolic links from tiltseries folder to average folder
        #    We need to do this because warp needs both folders: The tiltimages folder to get
        #    the header of the images (it only needs the first tiltimage of each tiltseries),
        #    and the averages folder to read them.
        if isTiltSeries:
            tiltSeriesFiles = glob.glob(os.path.join(imagesFolder, '*'))
            averagesFolder = os.path.join(imagesFolder, AVERAGE_FOLDER)
            pwutils.makePath(averagesFolder)
            for file in tiltSeriesFiles:
                fileName = os.path.basename(file)
                destFolder = os.path.join(averagesFolder, fileName)
                os.symlink('../' + fileName, destFolder)

            # 3. Create warp_tiltseries.settings file
            processingFolder = os.path.abspath(self._getExtraPath(TILTSERIES_FOLDER))
            if isinstance(objSet, TiltSeries):
                processingFolder = os.path.abspath(self._getExtraPath(TILTSERIES_FOLDER))
                pwutils.makePath(processingFolder)
            argsDict = {
                "--folder_data": os.path.abspath(self._getExtraPath(TOMOSTAR_FOLDER)),
                "--extension": "*.tomostar",
                "--folder_processing": processingFolder,
                "--bin": self.getBinFactor(),
                '--angpix': sr,
                "--output": os.path.abspath(self._getExtraPath(TILTSERIE_SETTINGS))
            }

            if exposure is not None:
                argsDict['--exposure'] = exposure

            if hasattr(self, 'tomo_thickness'):
                z = self.tomo_thickness.get()
                x = self.x_dimension.get() or objSet.getDimensions()[0]
                y = self.y_dimension.get() or objSet.getDimensions()[1]

                argsDict['--tomo_dimensions'] = f'{x}x{y}x{z}'

            cmd = ' '.join(['%s %s' % (k, v) for k, v in argsDict.items()])

            self.runJob(plugin.getProgram(CREATE_SETTINGS), cmd, executable='/bin/bash')

    @staticmethod
    def tomoStarGenerate(tsId, tiValues, otputFolder, isTiltSeries):
        """Generate the .tomostar files from TS"""
        _fileName = os.path.abspath(otputFolder) + '/%s.tomostar' % tsId
        _file = open(_fileName, 'a+')
        header = """
data_

loop_
_wrpMovieName #1
_wrpAngleTilt #2
_wrpAxisAngle #3
_wrpAxisOffsetX #4
_wrpAxisOffsetY #5
_wrpDose #6
_wrpAverageIntensity #7  
_wrpMaskedFraction #8
"""
        _file.write(header)

        sortedTiValues = sorted(tiValues)
        imagesFolder = TILTIMAGES_FOLDER if isTiltSeries else FRAMESERIES_FOLDER
        for acqOrder in sortedTiValues:
            value = tiValues[acqOrder]
            tiPath = '../%s/' % imagesFolder + value[0]
            angleTilt = value[1]
            axisAngle = value[2]
            shiftX = f"{value[3]:.6f}"
            shiftY = f"{value[4]:.6f}"
            dose = value[5]
            averageIntensity = 3.721
            maskedFraction = value[7]
            _file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                tiPath, angleTilt, axisAngle, shiftX, shiftY, dose, averageIntensity, maskedFraction))

        _file.close()


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
            streamingSleepOnWait: Some streaming protocols are quite fast,
                so, checking input/output updates creates an IO overhead.
                This params allows them to sleep (without consuming resources)
                to wait for new work to be done.
            streamingBatchSize: For some programs it is more efficient to process
                many items at once and not one by one. So this parameter will
                allow to group a number of items to be processed in the same
                protocol step. This can also reduce some IO overhead and spawning
                new OS processes.
        """
        form.addSection("Streaming")

        form.addParam("streamingSleepOnWait", IntParam, default=10,
                      label="Sleep when waiting (secs)",
                      help="If you specify a value greater than zero, "
                           "it will be the number of seconds that the "
                           "protocol will sleep when waiting for new "
                           "input data in streaming mode. ")
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

    def _getStreamingSleepOnWait(self):
        """ Retrieves the configured sleep duration for waiting during streaming.
            Returns:
            - int: The sleep duration in seconds during streaming wait.
            """
        return self.getAttributeValue('streamingSleepOnWait', 0)

    def _getStreamingBatchSize(self):
        """
            Retrieves the configured batch size for streaming operations.
            Returns:
            - int: The batch size for streaming operations.
            """
        return self.getAttributeValue('streamingBatchSize', 1)

    def _streamingSleepOnWait(self):
        """ This method should be used by protocols that want to sleep
        when there is not more work to do.
        """
        sleepOnWait = self._getStreamingSleepOnWait()
        if sleepOnWait > 0:
            self.info("Not much work to do now, sleeping %s seconds." % sleepOnWait)
            time.sleep(sleepOnWait)

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
