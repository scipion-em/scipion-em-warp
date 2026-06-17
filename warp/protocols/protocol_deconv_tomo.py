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
    """
    Applies Wiener-like deconvolution to a set of cryo-ET tomograms in order to improve the visibility of
    structural information that has been attenuated during image formation. The protocol is intended to
    enhance interpretability of three-dimensional density maps while preserving their biological context,
    making them more suitable for inspection, segmentation, particle localization, and downstream
    subtomogram analysis.

    AI Generated:

    Deconvolve Tomograms (ProtWarpDeconvTomo) - User Manual
        Overview

        The Deconvolve Tomograms protocol performs contrast restoration on reconstructed tomographic
        volumes using previously estimated contrast transfer information. Its principal goal is to improve
        the visual accessibility of structural features that are present in the tomogram but may be
        difficult to recognize because of attenuation introduced by the microscope imaging process.

        In cryo-electron tomography, this operation is especially useful because tomograms frequently
        contain weak contrast, crowded environments, and complex three-dimensional cellular context.
        Deconvolution helps reveal membranes, macromolecular assemblies, and other structural landmarks
        more clearly, facilitating biological interpretation without changing the underlying specimen.

        Biological Context and Typical Use

        For biological users, the protocol is particularly valuable during early tomogram exploration and
        annotation. It often improves the visibility of large assemblies, organelle boundaries, membrane
        layers, or intracellular densities that would otherwise remain difficult to interpret. This makes
        it relevant for both purified specimen tomography and in situ cellular studies.

        In many practical workflows, deconvolved tomograms are used as improved visual references for
        manual inspection, particle picking, segmentation, or target selection prior to more quantitative
        analysis. The purpose is not to create new information but to reveal existing experimental signal
        in a more accessible way.

        Inputs and Matching Strategy

        The protocol requires a set of tomograms together with corresponding tomographic contrast transfer
        estimations. Reliable biological interpretation depends on correct correspondence between the
        reconstructed volumes and the optical information used for restoration. When this correspondence is
        incomplete, only the tomograms with valid associated information can be processed confidently.

        For best results, the input tomograms should already represent reasonably well reconstructed
        volumes. Severe reconstruction artifacts, strong missing wedge effects, or highly noisy datasets
        may still limit the practical gain obtained from deconvolution.

        Deconvolution Behavior

        The restoration process is designed to improve contrast in a controlled manner. This balance is
        biologically important because tomograms are especially sensitive to noise amplification. Stronger
        enhancement may make densities appear sharper, but excessive aggressiveness can also emphasize
        reconstruction noise and generate visually misleading impressions of structural detail.

        In routine practice, moderate restoration often provides the most useful outcome. Features become
        easier to distinguish while preserving the natural appearance of tomographic densities. This is
        particularly helpful when the biological objective is localization or interpretation rather than
        high-resolution structural measurement.

        Outputs and Interpretation

        The protocol produces a new tomogram set that preserves the identity, geometry, and experimental
        meaning of the original volumes while presenting a contrast-optimized representation of the same
        specimen. These outputs remain suitable for downstream cryo-ET workflows and maintain their direct
        biological relation to the original acquisition.

        Users should interpret the resulting volumes as enhanced representations rather than as improved
        reconstructions in the strict structural sense. Clearer visibility does not necessarily imply
        increased resolution, and visually stronger densities should always be considered in the context of
        the original data.

        Practical Recommendations

        In biological workflows, deconvolution is often most useful after tomographic reconstruction and
        reliable contrast transfer estimation, but before segmentation, manual annotation, or subtomogram
        target identification. It can be particularly valuable in crowded cellular tomograms where local
        contrast is often the limiting factor for interpretation.

        Comparing restored tomograms with the original volumes is generally advisable. This preserves
        biological caution and helps prevent overinterpretation of visually enhanced noise or reconstruction
        artifacts. Conservative settings usually provide the best compromise between clarity and fidelity.

        Final Perspective

        For cryo-electron tomography users, tomogram deconvolution is best understood as a visibility
        enhancement procedure that improves practical interpretability of three-dimensional experimental
        data. When supported by reliable contrast transfer information and used with appropriate biological
        caution, it can significantly improve the usefulness of tomograms for exploration and downstream
        structural analysis.
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
