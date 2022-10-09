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

import pyworkflow.protocol.params as params
from pyworkflow.utils.properties import Message
from pwem.protocols import ProtMicrographs
from pwem.objects import SetOfMicrographs
from pwem.emlib.image import ImageHandler
from pwem.constants import RELATION_CTF

from ..utils import tom_deconv


class ProtWarpDeconv2D(ProtMicrographs):
    """ Protocol to deconvolve (Wiener-like filter) a set of micrographs.
    See https://github.com/dtegunov/tom_deconv
    """
    _label = 'deconvolve 2D'
    _possibleOutputs = {'outputMicrographs': SetOfMicrographs}

    def __init__(self, **kwargs):
        ProtMicrographs.__init__(self, **kwargs)
        #self.stepsExecutionMode = params.STEPS_PARALLEL

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputMicrographs',
                      params.PointerParam,
                      pointerClass='SetOfMicrographs',
                      label='Input micrographs',
                      important=True)
        form.addParam('ctfRelations', params.RelationParam,
                      important=True,
                      relationName=RELATION_CTF,
                      attributeName='getInputMicrographs',
                      label='CTF estimation',
                      help='Choose a CTF estimation '
                           'related to the input micrographs.')

        form.addParallelSection(threads=1, mpi=0)

    # --------------------------- STEPS functions -----------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.deconvolveStep)
        self._insertFunctionStep(self.createOutputStep)

    def deconvolveStep(self):
        # Load CTFs
        ctfDict = dict()
        for ctf in self.ctfRelations.get():
            micKey = ctf.getMicrograph().getMicName()
            ctfDict[micKey] = 0.5 * (ctf.getDefocusU() + ctf.getDefocusU())

        input_mics = self.getInputMicrographs()
        acq = input_mics.getAcquisition()
        pix = input_mics.getSamplingRate()
        voltage = acq.getVoltage()
        cs = acq.getSphericalAberration()

        micsList = input_mics.aggregate(["COUNT"], "_micName", ["_micName", "_filename"])
        ih = ImageHandler()

        # Iterate over mics
        for mic in micsList:
            micKey = mic["_micName"]
            if micKey in ctfDict:
                micName = mic["_filename"]
                defocus = ctfDict[micKey] / 10000
                inputData = ih.read(micName).getData()

                with mrcfile.new(self._getOutputMic(micName), overwrite=True) as mrcOut:
                    self.info(f"Deconvolving {micName}")
                    result = tom_deconv(inputData, pix, voltage, cs, defocus,
                                        ncpu=self.numberOfThreads.get())
                    mrcOut.set_data(result)
                    mrcOut.voxel_size = pix
                    mrcOut.update_header_from_data()
            else:
                self.info(f"No CTF found for mic: {micKey}")

    def createOutputStep(self):
        in_mics = self.getInputMicrographs()
        out_mics = self._createSetOfMicrographs()
        out_mics.copyInfo(in_mics)
        out_mics.copyItems(in_mics, updateItemCallback=self._updateItem)

        self._defineOutputs(outputMicrographs=out_mics)
        self._defineTransformRelation(self.getInputMicrographs(pointer=True),
                                      out_mics)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []

        return summary

    # -------------------------- UTILS functions ------------------------------
    def getInputMicrographs(self, pointer=False):
        if pointer:
            return self.inputMicrographs
        else:
            return self.inputMicrographs.get()

    def _getOutputMic(self, micName):
        return self._getExtraPath(os.path.basename(micName))

    def _updateItem(self, item, row):
        current_out_mic = self._getOutputMic(item.getFileName())
        if os.path.exists(current_out_mic):
            item.setFileName(current_out_mic)
        else:
            item._appendItem = False
