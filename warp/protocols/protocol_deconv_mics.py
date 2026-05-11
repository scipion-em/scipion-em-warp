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
from pwem.protocols import ProtMicrographs
from pwem.objects import SetOfMicrographs
from pwem.constants import RELATION_CTF

from warp.protocols.protocol_base import ProtWarpBase


class outputs(Enum):
    Micrographs = SetOfMicrographs


class ProtWarpDeconvMics(ProtWarpBase, ProtMicrographs):
    """
    Applies Wiener-like deconvolution to a set of cryo-EM micrographs in order to improve image interpretability
    before downstream analysis. The protocol is intended to enhance contrast while compensating for the optical
    effects introduced during image formation, providing micrographs that are better suited for particle picking,
    inspection, and subsequent reconstruction tasks.

    AI Generated:

    Deconvolve Micrographs (ProtWarpDeconvMics) - User Manual
        Overview

        The Deconvolve Micrographs protocol performs image restoration on cryo-EM micrographs by applying a
        deconvolution procedure guided by previously estimated contrast transfer information. Its primary objective
        is to recover visually meaningful signal that has been attenuated by microscope optics while limiting the
        amplification of noise. In practical cryo-EM workflows, this step is often useful when the user wants
        clearer particle boundaries, improved low-contrast visibility, or more interpretable images before
        downstream processing.

        For biological users, this protocol is especially valuable during the early stages of data assessment.
        Deconvolved micrographs frequently make particle populations easier to recognize, improve manual inspection,
        and can facilitate automatic particle detection in difficult datasets. It is important to understand that
        the purpose is not to alter biological content, but to produce a more informative representation of the
        recorded experimental signal.

        Inputs and Biological Context

        The protocol requires a set of input micrographs together with associated contrast transfer estimations.
        These estimations provide the physical description needed to compensate for microscope-induced modulation
        of the recorded images. For best biological reliability, the contrast transfer information should come
        from the same micrograph dataset and should be of reasonable quality. Poor defocus estimation or
        mismatched metadata can reduce the usefulness of the resulting images.

        In typical cryo-EM practice, the input micrographs should already have passed basic quality control.
        Extremely contaminated images, heavily drifted acquisitions, or strongly damaged exposures may still
        produce limited improvement after deconvolution because the underlying signal may already be compromised.

        Deconvolution Behavior

        The protocol uses a controlled deconvolution strategy that attempts to enhance structural information
        while avoiding excessive amplification of high-frequency noise. This balance is biologically important.
        Excessively aggressive enhancement may create visually sharp images that appear attractive but may
        exaggerate noise and produce misleading apparent detail. More conservative settings often preserve a
        more faithful representation of the experimental data.

        In many practical datasets, moderate enhancement is sufficient to reveal particle outlines, membrane
        boundaries, or large macromolecular features that were previously difficult to inspect. The main benefit
        is therefore interpretability rather than the creation of new information.

        Processing Considerations

        The protocol can operate efficiently on large micrograph collections and is suitable both for exploratory
        processing and routine production workflows. Depending on computational resources, execution may be
        performed using either standard processor resources or graphical acceleration. From the biological
        perspective, the computational mode does not change the meaning of the output. The resulting micrographs
        remain physically tied to the same experimental acquisition and preserve their identity within the dataset.

        Because the output remains linked to the original micrographs, this protocol integrates naturally into
        common cryo-EM pipelines. Users can continue with particle picking, visual quality assessment, or other
        preprocessing steps without changing the biological interpretation of the dataset.

        Outputs and Interpretation

        The protocol produces a new set of micrographs that preserve the acquisition context of the originals
        while providing an enhanced representation of the signal. These outputs are typically easier to inspect
        visually and may improve confidence during early decision-making stages of processing.

        Biologically, the user should interpret the output as a filtered and contrast-optimized representation
        of the same underlying specimen. Apparent sharpening should not automatically be interpreted as improved
        structural resolution. The real value lies in clearer visibility of existing information rather than
        in generating additional structural content.

        Practical Recommendations

        In routine cryo-EM work, this protocol is often most useful after reliable contrast transfer estimation
        and before particle picking or manual dataset inspection. It can be particularly helpful in low-contrast
        datasets, membrane protein samples, or cases where particles are difficult to distinguish from the
        surrounding background.

        When interpreting results, it is advisable to compare deconvolved images with the original micrographs.
        This helps maintain biological caution and prevents overinterpretation of visually enhanced features.
        Conservative enhancement usually provides the best balance between visibility and reliability.

        Final Perspective

        For many cryo-EM users, micrograph deconvolution is best understood as a visibility-enhancing preprocessing
        step. Its main contribution is to make experimentally recorded information easier to inspect and exploit
        during downstream analysis. When used with accurate contrast transfer information and interpreted with
        appropriate biological caution, it can substantially improve the practical usability of raw micrograph data.
    """
    _label = 'deconvolve micrographs'
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
        
        self.defineProcessParams(form)
        form.addParallelSection(threads=8, mpi=0)

    # --------------------------- STEPS functions -----------------------------
    def deconvolveStep(self):
        """ Load input CTFs and micrograph sets, then run deconvolution. """
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
        self._deconvolve(pix, acq, micsList, ctfDict, keyName="_micName")

    def createOutputStep(self):
        in_mics = self.getInputMicrographs()
        out_mics = self._createSetOfMicrographs()
        out_mics.copyInfo(in_mics)
        out_mics.copyItems(in_mics, doClone=False,
                           updateItemCallback=self._updateItem)

        self._defineOutputs(**{outputs.Micrographs.name: out_mics})
        self._defineTransformRelation(self.getInputMicrographs(pointer=True),
                                      out_mics)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []

        if hasattr(self, outputs.Micrographs.name):
            summary.append(f"Deconvolved {self.getInputMicrographs().getSize()} "
                           "micrographs")

        return summary

    # -------------------------- UTILS functions ------------------------------
    def getInputMicrographs(self, pointer=False):
        return self.inputMicrographs if pointer else self.inputMicrographs.get()
