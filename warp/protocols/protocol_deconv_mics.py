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
    Deconvolves a set of input micrographs using a Wiener-like filtering
    strategy based on previously estimated CTF parameters.

    AI Generated:

    Deconvolve Micrographs (ProtWarpDeconvMics) — User Manual
        Overview

        The ProtWarpDeconvMics protocol applies image deconvolution to a set
        of cryo-EM micrographs. Its goal is to partially compensate for the
        contrast transfer effects introduced by the microscope optics and to
        enhance interpretable signal before downstream processing.

        In practical cryo-EM workflows, this protocol is commonly used after
        motion correction and CTF estimation, when the user wants to improve
        particle visibility, facilitate particle picking, or produce cleaner
        micrographs for inspection and further analysis.

        Biological Purpose

        Micrographs recorded in cryo-EM are affected by contrast attenuation,
        especially at particular spatial frequencies. This protocol uses the
        CTF information associated with each micrograph to perform a controlled
        deconvolution.

        The result is not a reconstruction of the true specimen density, but
        an enhancement of useful image contrast that may improve the visibility
        of particles, macromolecular boundaries, and weak structural features.

        Inputs

        The protocol requires:

            1. A set of input micrographs.
            2. A related CTF estimation associated with those micrographs.

        Each input micrograph must have a corresponding CTF entry. The
        deconvolution process uses the defocus value extracted from the CTF
        metadata together with acquisition parameters such as:

            - pixel size
            - accelerating voltage
            - spherical aberration

        If a micrograph has no matching CTF information, it is skipped.

        Processing Strategy

        For every micrograph:

            - The protocol identifies the corresponding CTF entry.
            - It computes the effective defocus.
            - It generates an output filename.
            - It performs deconvolution using the Warp-based implementation
              inherited from ProtWarpBase.

        The actual numerical filtering is performed through tomographic
        deconvolution utilities that process the micrograph image directly.

        Processing Parameters

        The protocol exposes several advanced deconvolution parameters.

        Deconvolution strength
            Controls the intensity of the filter. Larger values enhance
            contrast more aggressively but may amplify noise.

        SNR falloff
            Regulates how rapidly high-frequency components are attenuated.
            It helps stabilize the deconvolution in noisy regions.

        High-pass fraction
            Suppresses very low frequencies that would otherwise be
            excessively boosted during deconvolution.

        These parameters should usually remain close to defaults unless
        the user has a specific reason to optimize them.

        CPU and GPU Execution

        The protocol supports both CPU and GPU execution.

        GPU mode can accelerate processing considerably for large datasets.
        Only one GPU is used per execution.

        If CPU execution is selected, multiple threads can be used.

        Workflow

        Step 1 — Deconvolution

            The protocol reads all input micrographs and their associated
            CTF metadata.

            For each micrograph:

                - retrieve micrograph name
                - retrieve defocus
                - apply deconvolution
                - write the resulting micrograph to the protocol output folder

        Step 2 — Output Creation

            Once all micrographs are processed:

                - a new SetOfMicrographs is created
                - metadata from the input set is preserved
                - file paths are updated to point to the deconvolved images

            Micrographs that failed processing are excluded automatically.

        Outputs

        The protocol generates:

            Micrographs
                A new set of deconvolved micrographs.

        The output preserves the original metadata and acquisition
        information while replacing the image filenames with the newly
        generated deconvolved files.

        Output Naming Convention

        Each output micrograph is written as:

            <original_name>_deconv.mrc

        This makes it easy to distinguish processed micrographs from
        the original input data.

        Practical Recommendations

        Typical use cases include:

            - improving particle visibility before picking
            - enhancing low-contrast datasets
            - inspecting micrographs before classification

        Recommended practice:

            - use reliable CTF estimations
            - keep default deconvolution parameters initially
            - visually inspect the results before continuing

        Over-aggressive deconvolution may increase high-frequency noise,
        so visual validation remains important.

        Summary

        ProtWarpDeconvMics provides a convenient micrograph-level
        deconvolution protocol integrated into Scipion.

        It combines:

            - input micrographs
            - associated CTF metadata
            - Warp-based deconvolution

        to generate a new micrograph set with enhanced contrast that
        can be used in downstream cryo-EM processing.
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
