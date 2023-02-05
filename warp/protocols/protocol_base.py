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

from ..utils import tom_deconv


class ProtWarpBase(EMProtocol):
    _label = None

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)

    # --------------------------- STEPS functions -----------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.deconvolveStep)
        self._insertFunctionStep(self.createOutputStep)

    def deconvolveStep(self):
        raise NotImplementedError

    def createOutputStep(self):
        raise NotImplementedError

    def _deconvolve(self, pixSize, acquisition, inputList,
                    ctfDict, keyName="_micName"):
        """ Main function to iterate over inputList items.
        Checks the match with ctfDict by keyName for each item.
        """
        voltage = acquisition.getVoltage()
        cs = acquisition.getSphericalAberration()
        ih = ImageHandler()

        for item in inputList:
            key = item[keyName]
            if key in ctfDict:
                fileName = item["_filename"]
                defocus = ctfDict[key] / 10000
                inputData = ih.read(fileName).getData()

                with mrcfile.new(self._getOutputFn(fileName), overwrite=True) as mrcOut:
                    self.info(f"Deconvolving {fileName}")
                    result = tom_deconv(inputData, pixSize, voltage, cs, defocus,
                                        ncpu=self.numberOfThreads.get())
                    mrcOut.set_data(result)
                    mrcOut.voxel_size = pixSize
                    mrcOut.update_header_from_data()
            else:
                self.warning(f"No CTF found for: {key}")

    # -------------------------- UTILS functions ------------------------------
    def getInputMicrographs(self, pointer=False):
        if pointer:
            return self.inputMicrographs
        else:
            return self.inputMicrographs.get()

    def _getOutputFn(self, micName):
        return self._getExtraPath(pwutils.removeBaseExt(micName) + "_deconv.mrc")

    def _updateItem(self, item, row):
        outputFn = self._getOutputFn(item.getFileName())
        if os.path.exists(outputFn):
            item.setFileName(outputFn)
        else:
            item._appendItem = False
