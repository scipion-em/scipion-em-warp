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
import mrcfile
import glob

from pyworkflow.protocol import FloatParam, Positive, LEVEL_ADVANCED
import pyworkflow.utils as pwutils
from pwem.protocols import EMProtocol
from pwem.emlib.image import ImageHandler
from pwem.emlib.image.image_readers import ImageStack, ImageReadersRegistry
from warp import CREATE_SETTINGS

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
        stack_shape = (n, y, x)
        mrc = mrcfile.new_mmap(outputFn, shape=stack_shape,
                               mrc_mode=2, overwrite=True)
        for i in range(n):
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

    def dataPrepare(self, setOfTiltSeries):
        """Creates the setting file that will be used by the different programs.
           It also extracts the tiltimages from the tiltseries and generates the *.tomostar files based on
           the tiltimages."""
        starFolder = self._getExtraPath('tomostar')
        pwutils.makePath(starFolder)
        tiltimagesFolder = self._getExtraPath('tiltimages')
        pwutils.makePath(tiltimagesFolder)

        # 1. Extract all tilt images from the tiltseries
        for ts in setOfTiltSeries.iterItems():
            tsId = ts.getTsId()
            properties = {"sr": ts.getSamplingRate()}
            index = 1
            tiValues = {}
            for ti in ts.iterItems():
                newBinaryName = tsId + '_%s_%s.mrc' % (index, ti.getTiltAngle())
                newFrame = ImageStack(properties=properties)
                newFrame.append(ImageReadersRegistry.open(str(index) + '@' + ti.getFileName()))
                ImageReadersRegistry.write(newFrame, os.path.join(self._getExtraPath('tiltimages'), newBinaryName))
                # Taking angleTilt axisAngle Dose AverageIntensity MaskedFraction
                dose = 0
                if ts.hasAcquisition():
                    axisAngle = ts._acquisition.getTiltAxisAngle()
                else:
                    axisAngle = 0
                amplitudeContrast = 0
                maskedFraction = 0
                if ti.getAcquisition():
                    amplitudeContrast = ti.getAcquisition().getAmplitudeContrast()
                    dose = ti.getAcquisition().getDosePerFrame()

                tiValues[newBinaryName] = [ti.getTiltAngle(), axisAngle, dose, amplitudeContrast, maskedFraction]

                index += 1

            self.tomoStarGenerate(tsId, tiValues, starFolder)

        # 2. Create symbolic links from tiltseries folder to average folder
        #    We need to do this because warp needs both folders: The tiltimages folder to get
        #    the header of the images (it only needs the first tiltimage of each tiltseries),
        #    and the averages folder to read them.
        tiltSeriesFiles = glob.glob(os.path.join(tiltimagesFolder, '*'))
        averagesFolder = os.path.join(tiltimagesFolder, 'average')
        pwutils.makePath(averagesFolder)
        for file in tiltSeriesFiles:
            fileName = os.path.basename(file)
            destFolder = os.path.join(averagesFolder, fileName)
            os.symlink('../' + fileName, destFolder)

        # 3. Create warp_tiltseries.settings file
        self.processingFolder = os.path.abspath(self._getExtraPath('warp_tiltseries'))
        pwutils.makePath(self.processingFolder)
        argsDict = {
            "--folder_data": os.path.abspath(self._getExtraPath('tomostar')),
            "--extension": "*.tomostar",
            "--folder_processing": self.processingFolder,
            '--angpix': 0.7894,
            '--exposure': 2.64,
            "--output": os.path.abspath(self._getExtraPath("warp_tiltseries.settings")),
        }

        cmd = ' '.join(['%s %s' % (k, v) for k, v in argsDict.items()])
        self.runJob(self.getPlugin().getProgram(CREATE_SETTINGS), cmd, executable='/bin/bash')

    @staticmethod
    def tomoStarGenerate(tsId, tiValues, otputFolder):
        """Generate the .tomostar files from TS"""
        _fileName = os.path.abspath(otputFolder) + '/%s.tomostar' % tsId
        _file = open(_fileName, 'a+')
        header = """
data_

loop_
_wrpMovieName #1
_wrpAngleTilt #2
_wrpAxisAngle #3
_wrpDose #4
_wrpAverageIntensity #5
_wrpMaskedFraction #6
"""
        _file.write(header)

        for key, value in tiValues.items():
            tiPath = '../tiltimages/' + key
            angleTilt = value[0]
            axisAngle = value[1]
            dose = value[2]
            averageIntensity = value[3]
            maskedFraction = value[4]
            _file.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (
                tiPath, angleTilt, axisAngle, dose, averageIntensity, maskedFraction))

        _file.close()



