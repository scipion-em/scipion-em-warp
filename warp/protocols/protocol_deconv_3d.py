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
from pwem.emlib.image import ImageHandler
from pwem.protocols import EMProtocol

from tomo.protocols import ProtTomoBase
from tomo.objects import SetOfTomograms

from ..utils import tom_deconv


class ProtWarpDeconv3D(EMProtocol, ProtTomoBase):
    """ Protocol to deconvolve (Wiener-like filter) a set of tomograms.
    See https://github.com/dtegunov/tom_deconv
    """
    _label = 'deconvolve 3D'
    _possibleOutputs = {'outputTomograms': SetOfTomograms}

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = params.STEPS_PARALLEL

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputTomograms',
                      params.PointerParam,
                      pointerClass='SetOfTomograms',
                      label='Input tomograms',
                      important=True)
        form.addParam('inputCTFs',
                      params.PointerParam,
                      pointerClass='SetOfCTFTomoSeries',
                      important=True,
                      label='CTF tomo series',
                      help='Set of CTFs that correspond to the '
                           'input above. The matching is done using tsId.')

        form.addParallelSection(threads=1, mpi=1)

    # --------------------------- STEPS functions -----------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.deconvolveStep)
        self._insertFunctionStep(self.createOutputStep)

    def deconvolveStep(self):
        tomoSet = self.getInputTomos()
        ctfSet = self.inputCTFs.get()
        tsIds_from_tomos = set(item.getTsId() for item in tomoSet)
        tsIds_from_ctfs = set(item.getTsId() for item in ctfSet)
        if not tsIds_from_ctfs.issubset(tsIds_from_tomos):
            self.warning("Found CTFs with tsId that did not match "
                         "provided tomograms: "
                         f"{set.difference(tsIds_from_ctfs, tsIds_from_tomos)}")

        # Load CTFs
        ctfDict = dict()
        for ctfSeries in ctfSet.iterItems():
            tsKey = ctfSeries.getTsId()
            ctfValues = [0.5 * (ctf.getDefocusU() + ctf.getDefocusU()) for ctf in ctfSeries]
            ctfDict[tsKey] = sum(ctfValues) / len(ctfValues)

        acq = tomoSet.getAcquisition()
        pix = tomoSet.getSamplingRate()
        voltage = acq.getVoltage()
        cs = acq.getSphericalAberration()

        tomoList = tomoSet.aggregate(["COUNT"], "_tsId", ["_tsId", "_filename"])
        ih = ImageHandler()

        # Iterate over tomos
        for tomo in tomoList:
            tsKey = tomo["_tsId"]
            if tsKey in ctfDict:
                micName = tomo["_filename"]
                defocus = ctfDict[tsKey] / 10000
                inputData = ih.read(micName).getData()

                with mrcfile.new(self._getOutputTomo(micName), overwrite=True) as mrcOut:
                    self.info(f"Deconvolving {micName}")
                    result = tom_deconv(inputData, pix, voltage, cs, defocus)
                    mrcOut.set_data(result)
                    mrcOut.update_header_from_data()
            else:
                self.info(f"No CTF found for tomo: {tsKey}")

    def createOutputStep(self):
        in_tomos = self.getInputTomos()
        out_tomos = self._createSetOfTomograms()
        out_tomos.copyInfo(in_tomos)
        out_tomos.copyItems(in_tomos, updateItemCallback=self._updateItem)

        self._defineOutputs(outputTomograms=out_tomos)
        self._defineTransformRelation(self.getInputTomos(pointer=True),
                                      out_tomos)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []

        return summary

    # -------------------------- UTILS functions ------------------------------
    def getInputTomos(self, pointer=False):
        if pointer:
            return self.inputTomograms
        else:
            return self.inputTomograms.get()

    def _getOutputTomo(self, micName):
        return self._getExtraPath(os.path.basename(micName))

    def _updateItem(self, item, row):
        current_out_fn = self._getOutputTomo(item.getFileName())
        if os.path.exists(current_out_fn):
            item.setLocation(current_out_fn)
        else:
            item._appendItem = False
