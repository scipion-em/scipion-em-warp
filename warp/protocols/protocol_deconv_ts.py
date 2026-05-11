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

from pyworkflow.constants import NEW
import pyworkflow.protocol.params as params
from pyworkflow.utils.properties import Message
import pyworkflow.utils as pwutils

from tomo.protocols import ProtTomoBase
from tomo.objects import SetOfTiltSeries

from warp.protocols.protocol_base import ProtWarpBase


class outputs(Enum):
    TiltSeries = SetOfTiltSeries


class ProtWarpDeconvTS(ProtWarpBase, ProtTomoBase):
    """
    Deconvolves a set of tilt-series using a Wiener-like filtering
    strategy based on the average CTF information associated with each
    tilt series.

    AI Generated:

    Deconvolve Tilt-Series (ProtWarpDeconvTS) — User Manual
        Overview

        The ProtWarpDeconvTS protocol applies deconvolution directly to
        aligned tilt-series stacks before tomogram reconstruction. Its
        purpose is to compensate for contrast attenuation introduced by
        microscope optics and improve the quality of the tilt images that
        will later contribute to tomographic reconstruction.

        In cryo-electron tomography workflows, this protocol is typically
        used after tilt-series alignment and CTF estimation, but before
        tomogram reconstruction.

        Biological Purpose

        Individual tilt images often suffer from low contrast and reduced
        visibility of structural features. Applying deconvolution at the
        tilt-series level enhances interpretable signal before volume
        reconstruction.

        Since reconstruction integrates information from all tilts,
        improving the contrast of the individual projections can lead to
        better visual quality and potentially more interpretable tomograms.

        The protocol does not correct missing-wedge effects or alignment
        errors. It only applies a controlled frequency-domain
        deconvolution.

        Inputs

        The protocol requires:

            1. A set of input tilt-series.
            2. A corresponding SetOfCTFTomoSeries.

        Matching between tilt-series and CTF metadata is performed using
        the tilt-series identifier (tsId).

        For each tilt-series:

            - the first tilt image file is used as the stack reference
            - the corresponding CTF series is retrieved
            - the average defocus is computed from all CTF entries

        Tilt-series without matching CTF information are skipped.

        Processing Strategy

        During execution:

            - all input tilt-series are indexed by tsId
            - all CTF series are indexed by tsId
            - only matched identifiers are processed

        For every matched tilt-series:

            - the average defocus is calculated
            - acquisition metadata is collected
            - the complete image stack is deconvolved using the stack
              processing mode inherited from ProtWarpBase

        Unlike micrograph or tomogram deconvolution, here the protocol
        processes a full stack of images rather than a single 2D image
        or 3D volume.

        Processing Parameters

        The protocol inherits the deconvolution parameters from
        ProtWarpBase.

        Deconvolution strength
            Controls the aggressiveness of the filter.

        SNR falloff
            Stabilizes high-frequency behavior in noisy regions.

        High-pass fraction
            Prevents excessive boosting of very low frequencies.

        In most cases, default values provide a good starting point.

        CPU and GPU Execution

        The protocol supports both CPU and GPU execution.

        GPU execution can significantly accelerate stack processing.
        Only one GPU is used per execution.

        CPU execution allows multithreaded processing.

        Workflow

        Step 1 — Matching Tilt-Series and CTF Metadata

            The protocol creates two internal dictionaries:

                - input tilt-series
                - input CTF series

            Matching is done using tsId.

        Step 2 — Average Defocus Calculation

            For each matched tilt-series:

                - all CTF defocus values are collected
                - the mean defocus is calculated

        Step 3 — Stack Deconvolution

            The complete tilt-image stack is processed as a single unit.

            Each tilt image in the stack is deconvolved and written into
            a new output stack.

        Step 4 — Output Creation

            After processing:

                - a new SetOfTiltSeries is created
                - original metadata is preserved
                - every tilt image is updated to reference the new
                  deconvolved stack file

        Outputs

        The protocol generates:

            TiltSeries
                A new set of deconvolved tilt-series.

        Output Naming Convention

        Each deconvolved stack is written as:

            <original_name>_deconv.mrcs

        The output uses .mrcs format because the protocol processes
        stacked tilt-image data.

        Practical Recommendations

        This protocol is useful when:

            - improving contrast before tomogram reconstruction
            - preparing tilt-series for more interpretable reconstructions
            - enhancing weak projection data in low-contrast datasets

        Recommended practice:

            - ensure correct tsId matching between tilt-series and CTF series
            - use reliable CTF estimations
            - inspect a few deconvolved tilt images before reconstruction

        Since the protocol applies one average defocus per tilt-series,
        it provides a practical global correction rather than a
        tilt-by-tilt refinement.

        Summary

        ProtWarpDeconvTS provides a tilt-series-level deconvolution
        workflow integrated into Scipion.

        It combines:

            - aligned tilt-series stacks
            - associated CTF tilt-series information
            - Warp-based stack deconvolution

        to generate deconvolved tilt-series that can be used as improved
        input for tomographic reconstruction.
    """
    _label = 'deconvolve tilt-series'
    _possibleOutputs = outputs
    _devStatus = NEW

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
        
        self.defineProcessParams(form)
        form.addParallelSection(threads=8, mpi=0)

    # --------------------------- STEPS functions -----------------------------
    def deconvolveStep(self):
        """ Load CTF and TS sets, match by tsId before processing. """
        tsSet = self.getInputTS()
        ctfSet = self.inputCTFs.get()
        acq = tsSet.getAcquisition()
        pix = tsSet.getSamplingRate()

        tsDict, ctfDict = {}, {}

        for ts in tsSet.iterItems():
            tsDict[ts.getTsId()] = {
                "_tsId": ts.getTsId(),
                "_filename": ts.getFirstItem().getFileName()
            }

        for ctfSeries in ctfSet.iterItems():
            ctfValues = [0.5 * (ctf.getDefocusU() + ctf.getDefocusU()) for ctf in ctfSeries]
            ctfDict[ctfSeries.getTsId()] = sum(ctfValues) / len(ctfValues)

        matchIds = tsDict.keys() & ctfDict.keys()
        matchTs = [tsDict[i] for i in matchIds]
        mismatchIds = tsDict.keys() - ctfDict.keys()

        if mismatchIds:
            self.warning("No CTFs found for tilt-series with tsId: "
                         f"{mismatchIds}")

        # Iterate over TS
        self._deconvolve(pix, acq, matchTs, ctfDict, keyName="_tsId", isTS=True)

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

        if hasattr(self, outputs.TiltSeries.name):
            summary.append(f"Deconvolved {self.getInputTS().getSize()} "
                           "tilt-series")

        return summary

    # -------------------------- UTILS functions ------------------------------
    def getInputTS(self, pointer=False):
        return self.inputTiltSeries if pointer else self.inputTiltSeries.get()

    def updateTi(self, j, ts, ti, tsOut, tiOut):
        fn = ti.getFileName()
        tiOut.setFileName(self._getOutputFn(fn))

    def _getOutputFn(self, micName):
        """ Overwrite base class method to output mrcs instead. """
        return self._getExtraPath(pwutils.removeBaseExt(micName) + "_deconv.mrcs")
