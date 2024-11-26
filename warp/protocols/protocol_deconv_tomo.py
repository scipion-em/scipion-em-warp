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
from enum import Enum

from pyworkflow.constants import PROD
import pyworkflow.protocol.params as params
from pyworkflow.utils.properties import Message

from tomo.protocols import ProtTomoBase
from tomo.objects import SetOfTomograms

from warp.protocols.protocol_base import ProtWarpBase


class outputs(Enum):
    Tomograms = SetOfTomograms


class ProtWarpDeconvTomo(ProtWarpBase, ProtTomoBase):
    """ Protocol to deconvolve (Wiener-like filter) a set of tomograms.
    """
    _label = 'deconvolve tomograms'
    _possibleOutputs = outputs
    _devStatus = PROD

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addHidden(params.USE_GPU, params.BooleanParam,
                       default=False,
                       label="Use GPU for execution",
                       help="This protocol has both CPU and GPU implementation. "
                            "Select the one you want to use.")
        form.addHidden(params.GPU_LIST, params.StringParam, default='0',
                       label="Choose GPU ID",
                       help="GPU may have several cores. Set it to zero"
                            " if you do not know what we are talking about."
                            " First core index is 0, second 1 and so on."
                            " You can use only single GPU.")
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

        self.defineProcessParams(form)
        form.addParallelSection(threads=8, mpi=0)
        
    # --------------------------- STEPS functions -----------------------------
    def deconvolveStep(self):
        """ Load CTF and tomograms sets, match by tsId before processing. """
        tomoSet = self.getInputTomos()
        ctfSet = self.inputCTFs.get()
        acq = tomoSet.getAcquisition()
        pix = tomoSet.getSamplingRate()

        tsDict, ctfDict = {}, {}

        for tomo in tomoSet.iterItems():
            tsDict[tomo.getTsId()] = {
                "_tsId": tomo.getTsId(),
                "_filename": tomo.getFileName()
            }

        for ctfSeries in ctfSet.iterItems():
            ctfValues = [0.5 * (ctf.getDefocusU() + ctf.getDefocusV()) for ctf in ctfSeries]
            ctfDict[ctfSeries.getTsId()] = sum(ctfValues) / len(ctfValues)

        matchIds = tsDict.keys() & ctfDict.keys()
        matchTs = [tsDict[i] for i in matchIds]
        mismatchIds = tsDict.keys() - ctfDict.keys()

        if mismatchIds:
            self.warning("No CTFs found for tomograms with tsId: "
                         f"{mismatchIds}")

        # Iterate over tomos
        self._deconvolve(pix, acq, matchTs, ctfDict, keyName="_tsId")

    def createOutputStep(self):
        in_tomos = self.getInputTomos()
        out_tomos = self._createSetOfTomograms()
        out_tomos.copyInfo(in_tomos)
        out_tomos.copyItems(in_tomos, doClone=False,
                            updateItemCallback=self._updateItem)

        self._defineOutputs(**{outputs.Tomograms.name: out_tomos})
        self._defineTransformRelation(self.getInputTomos(pointer=True),
                                      out_tomos)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []

        if hasattr(self, outputs.Tomograms.name):
            summary.append(f"Deconvolved {self.getInputTomos().getSize()} "
                           "tomograms")

        return summary

    # -------------------------- UTILS functions ------------------------------
    def getInputTomos(self, pointer=False):
        return self.inputTomograms if pointer else self.inputTomograms.get()
