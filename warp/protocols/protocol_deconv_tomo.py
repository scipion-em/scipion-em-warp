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
    Deconvolves a set of tomograms using a Wiener-like filtering strategy
    based on the average CTF information associated with each tilt series.

    AI Generated:

    Deconvolve Tomograms (ProtWarpDeconvTomo) — User Manual
        Overview

        The ProtWarpDeconvTomo protocol applies deconvolution to a set of
        reconstructed tomograms. Its purpose is to partially compensate for
        contrast attenuation introduced by the microscope optics and improve
        the interpretability of tomographic volumes.

        In cryo-electron tomography workflows, this protocol is typically
        applied after tomogram reconstruction and CTF estimation, when the
        user wants to enhance structural contrast before particle extraction,
        subtomogram averaging, segmentation, or visual inspection.

        Biological Purpose

        Tomograms often exhibit reduced contrast, especially at high spatial
        frequencies. This protocol uses CTF information derived from the
        associated tilt series to apply a controlled deconvolution to the
        final reconstructed volume.

        The result is not a correction of missing-wedge artifacts or a full
        recovery of specimen density. Instead, it enhances useful contrast
        and improves the visibility of macromolecular features within the
        tomographic volume.

        Inputs

        The protocol requires:

            1. A set of input tomograms.
            2. A corresponding SetOfCTFTomoSeries.

        Matching between tomograms and CTF information is performed using
        the tilt-series identifier (tsId).

        For each tomogram:

            - the tomogram file is retrieved
            - the corresponding CTF series is retrieved
            - the average defocus is computed from all CTF entries
              in that series

        If a tomogram has no associated CTF series, it is skipped.

        Processing Strategy

        During execution:

            - all input tomograms are indexed by tsId
            - all CTF series are indexed by tsId
            - only matching identifiers are processed

        For each matched tomogram:

            - the average defocus of the full tilt series is calculated
            - acquisition parameters are collected
            - the tomogram is deconvolved using the Warp implementation
              inherited from ProtWarpBase

        If mismatches are found, the protocol reports them as warnings.

        Processing Parameters

        The protocol inherits the same advanced deconvolution parameters
        available in ProtWarpBase.

        Deconvolution strength
            Controls the strength of contrast enhancement.

        SNR falloff
            Determines the attenuation of noisy high-frequency signal.

        High-pass fraction
            Suppresses very low frequencies that would otherwise be
            excessively amplified.

        In most practical workflows, the default values are appropriate.

        CPU and GPU Execution

        The protocol supports both CPU and GPU execution.

        GPU execution can accelerate processing significantly for large
        tomograms. Only one GPU is used per execution.

        CPU execution allows multithreaded processing.

        Workflow

        Step 1 — Matching Tomograms and CTF Series

            The protocol builds internal dictionaries for:

                - tomograms
                - CTF series

            Matching is done by tilt-series identifier.

        Step 2 — Average Defocus Calculation

            For each matched CTF series:

                - all defocus values are collected
                - the mean defocus is computed

        Step 3 — Deconvolution

            For each matched tomogram:

                - the tomogram is processed
                - the deconvolved volume is written to the output folder

        Step 4 — Output Creation

            Once processing finishes:

                - a new SetOfTomograms is created
                - metadata from the input set is preserved
                - file paths are updated to point to deconvolved tomograms

            Tomograms that failed processing are excluded automatically.

        Outputs

        The protocol generates:

            Tomograms
                A new set of deconvolved tomograms.

        Output Naming Convention

        Each output tomogram is written as:

            <original_name>_deconv.mrc

        This preserves traceability between original and processed data.

        Practical Recommendations

        This protocol is especially useful when:

            - preparing tomograms for subtomogram particle picking
            - improving visual inspection of weak densities
            - enhancing contrast before segmentation or annotation

        Recommended practice:

            - ensure correct tsId matching between tomograms and CTF series
            - use reliable CTF estimations
            - visually inspect deconvolved tomograms before downstream analysis

        Because the protocol uses the average defocus of the tilt series,
        it provides a practical global correction rather than a per-tilt
        refinement.

        Summary

        ProtWarpDeconvTomo provides a convenient tomogram-level
        deconvolution workflow integrated into Scipion.

        It combines:

            - reconstructed tomograms
            - associated CTF tilt-series information
            - Warp-based deconvolution

        to generate a new tomogram set with enhanced contrast for
        downstream cryo-ET analysis.
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
