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

import pyworkflow.protocol.params as params
from pyworkflow.utils.properties import Message
from pwem.protocols import ProtMicrographs
from pwem.objects import SetOfMicrographs
from pwem.constants import RELATION_CTF

from .protocol_base import ProtWarpBase


class ProtWarpDeconv2D(ProtWarpBase, ProtMicrographs):
    """ Protocol to deconvolve (Wiener-like filter) a set of micrographs.
    """
    _label = 'deconvolve 2D'
    _possibleOutputs = {'outputMicrographs': SetOfMicrographs}

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

        form.addParallelSection(threads=8, mpi=0)

    # --------------------------- STEPS functions -----------------------------
    def deconvolveStep(self):
        # Load CTFs
        ctfDict = dict()
        for ctf in self.ctfRelations.get():
            micKey = ctf.getMicrograph().getMicName()
            ctfDict[micKey] = 0.5 * (ctf.getDefocusU() + ctf.getDefocusU())

        input_mics = self.getInputMicrographs()
        acq = input_mics.getAcquisition()
        pix = input_mics.getSamplingRate()
        micsList = input_mics.aggregate(["COUNT"], "_micName",
                                        ["_micName", "_filename"])

        # Iterate over mics
        self._deconvolve(pix, acq, micsList, ctfDict,
                         keyName="_micName")

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

        if self.isFinished():
            summary.append(f"Deconvolved {self.getInputMicrographs().getSize()} "
                           "micrographs")

        return summary
