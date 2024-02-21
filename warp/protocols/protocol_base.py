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

import pyworkflow.utils as pwutils
from pwem.protocols import EMProtocol
from pwem.emlib.image import ImageHandler

from warp.utils import tom_deconv


class ProtWarpBase(EMProtocol):
    _label = None

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

        for item in inputList:
            key = item[keyName]
            if key in ctfDict:
                fileName = item["_filename"]
                defocus = ctfDict[key] / 10000
                outputFn = self._getOutputFn(fileName)
                self.info(f"Deconvolving {fileName}")

                func = self._processStack if isTS else self._processImage
                func(fileName, outputFn, angpix=pixSize, voltage=voltage,
                     cs=cs, defocus=defocus, ncpu=ncpu, gpu=gpu, gpuid=gpuid)
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

    @classmethod
    def _processImage(cls, inputFn, outputFn, **kwargs):
        ih = ImageHandler()
        inputData = ih.read(inputFn).getData()
        with mrcfile.new(outputFn) as mrcOut:
            result = tom_deconv(inputData, **kwargs)
            mrcOut.set_data(result)
            mrcOut.voxel_size = kwargs["angpix"]
            mrcOut.update_header_from_data()

    @classmethod
    def _processStack(cls, inputFn, outputFn, **kwargs):
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
