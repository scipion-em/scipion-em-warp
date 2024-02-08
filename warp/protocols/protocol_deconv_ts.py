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

import pyworkflow.protocol.params as params
from pyworkflow.utils.properties import Message
import pyworkflow.utils as pwutils

from tomo.protocols import ProtTomoBase
from tomo.objects import SetOfTiltSeries

from warp.protocols.protocol_base import ProtWarpBase


class outputs(Enum):
    TiltSeries = SetOfTiltSeries


class ProtWarpDeconvTS(ProtWarpBase, ProtTomoBase):
    """ Protocol to deconvolve (Wiener-like filter) a set of tilt-series.
    """
    _label = 'deconvolve tilt-series'
    _possibleOutputs = outputs

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
        form.addParam('inputTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      label='Input tilt-series',
                      important=True)
        form.addParam('inputCTFs',
                      params.PointerParam,
                      pointerClass='SetOfCTFTomoSeries',
                      important=True,
                      label='CTF tomo series',
                      help='Set of CTFs that correspond to the '
                           'input above. The matching is done using tsId.')

        form.addParallelSection(threads=8, mpi=0)

    # --------------------------- STEPS functions -----------------------------
    def deconvolveStep(self):
        tsSet = self.getInputTS()
        ctfSet = self.inputCTFs.get()

        tsIds_from_ts = self._getTsIds(tsSet)
        tsIds_from_ctfs = self._getTsIds(ctfSet)

        if not tsIds_from_ctfs.issubset(tsIds_from_ts):
            self.warning("Found CTFs with tsId that did not match "
                         "provided tilt-series: "
                         f"{set.difference(tsIds_from_ctfs, tsIds_from_ts)}")

        # Load CTFs
        ctfDict = dict()
        for ctfSeries in ctfSet.iterItems():
            tsKey = ctfSeries.getTsId()
            ctfValues = [0.5 * (ctf.getDefocusU() + ctf.getDefocusU()) for ctf in ctfSeries]
            ctfDict[tsKey] = sum(ctfValues) / len(ctfValues)

        acq = tsSet.getAcquisition()
        pix = tsSet.getSamplingRate()
        tsList = [{
            "_tsId": tsId,
            "_filename": tsSet.getTiltSeriesFromTsId(tsId).getFirstItem().getFileName()
        } for tsId in tsIds_from_ts]

        # Iterate over TS
        self._deconvolve(pix, acq, tsList, ctfDict, keyName="_tsId", isTS=True)

    def createOutputStep(self):
        in_ts = self.getInputTS()
        out_ts = self._createSetOfTiltSeries()
        out_ts.copyInfo(in_ts)
        out_ts.copyItems(in_ts, updateTiCallback=self.updateTi)

        self._defineOutputs(**{outputs.TiltSeries.name: out_ts})
        self._defineTransformRelation(self.getInputTS(pointer=True), out_ts)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []

        if self.isFinished():
            summary.append(f"Deconvolved {self.getInputTS().getSize()} "
                           "tilt-series")

        return summary

    # -------------------------- UTILS functions ------------------------------
    def getInputTS(self, pointer=False):
        if pointer:
            return self.inputTiltSeries
        else:
            return self.inputTiltSeries.get()

    def updateTi(self, j, ts, ti, tsOut, tiOut):
        fn = ti.getFileName()
        tiOut.setFileName(self._getOutputFn(fn))

    def _getOutputFn(self, micName):
        """ Overwrite base class method to output mrcs instead. """
        return self._getExtraPath(pwutils.removeBaseExt(micName) + "_deconv.mrcs")
