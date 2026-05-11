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
    Applies Wiener-like deconvolution to a set of tilt-series in order to
    improve contrast and recover interpretable signal before downstream
    tomographic analysis.

    AI Generated:

    Deconvolve Tilt-Series (ProtWarpDeconvTS) - User Manual
        Overview

        The Deconvolve Tilt-Series protocol is designed to improve the
        interpretability of cryo-electron tomography tilt-series before
        reconstruction or subsequent processing. Its main purpose is to
        compensate for contrast attenuation introduced during image
        formation so that structural information becomes more visible
        across the individual projections that compose each tilt-series.

        In practical cryo-ET workflows, this operation is especially
        useful when the original projections appear weak in contrast or
        when downstream alignment, reconstruction, denoising, or particle
        picking benefit from clearer signal. By enhancing the visibility
        of meaningful structural features while preserving the geometric
        relationship between projections, the protocol prepares tilt-series
        for more stable and biologically informative analysis.

        Inputs and Data Consistency

        The protocol requires a set of experimental tilt-series together
        with a corresponding set of contrast transfer function estimates.
        These inputs are expected to represent the same acquisition
        collection so that every tilt-series can be associated with its
        corresponding optical characterization.

        For biological applications, consistency between the imaging data
        and the optical estimations is particularly important. Reliable
        deconvolution assumes that the experimental conditions represented
        by the contrast transfer information accurately reflect the
        acquisition of the projections. When this correspondence is not
        available, processing may remain incomplete for some entries.

        Biological Motivation

        In cryo-electron tomography, biological structures are often
        embedded in crowded and noisy cellular environments. As a result,
        subtle molecular boundaries, membrane contours, and low-contrast
        macromolecular assemblies can become difficult to recognize
        directly in raw projections.

        Deconvolution helps recover visually meaningful information before
        tomographic reconstruction. For cellular studies, this can improve
        the interpretability of organelle boundaries, filament systems,
        membrane-associated complexes, and large molecular assemblies. For
        purified or subtomogram-oriented workflows, improved projection
        contrast can facilitate later alignment and reconstruction steps.

        Relationship to Tomographic Reconstruction

        Although this protocol does not reconstruct a tomogram itself, it
        plays an important preparatory role. The quality of a tomographic
        reconstruction is strongly influenced by the quality of the
        underlying projections. Better contrast at the tilt-series level
        often translates into cleaner reconstructed densities and more
        reliable downstream interpretation.

        Biological users should understand that deconvolution does not
        create new structural information. Instead, it improves the
        recoverability of information already present in the experimental
        images. For this reason, the protocol is best viewed as an
        enhancement stage that supports later biological interpretation
        rather than as a substitute for careful reconstruction or
        validation.

        Practical Use in Cryo-ET Workflows

        This protocol is particularly valuable before tomogram generation,
        especially when working with low-dose acquisitions or challenging
        cellular specimens. In exploratory projects, it can improve visual
        inspection of projection data and provide a clearer first estimate
        of specimen quality.

        In more advanced workflows, deconvolved tilt-series can become a
        stronger starting point for alignment refinement, tomographic
        reconstruction, and subtomogram analysis. When handling large
        datasets, the protocol also fits naturally into automated
        preprocessing pipelines intended to standardize projection quality
        before downstream computation.

        Outputs and Interpretation

        The protocol produces a new set of tilt-series that preserve the
        original experimental organization while providing projections with
        enhanced interpretability. The output remains directly compatible
        with subsequent tomographic processing steps.

        From a biological perspective, the resulting projections should be
        interpreted as improved representations of the same underlying
        specimen rather than as transformed biological states. Relative
        geometry, tilt ordering, and dataset identity remain unchanged,
        which allows the processed data to be used naturally in later
        stages of analysis.

        Final Perspective

        For most cryo-electron tomography users, deconvolution at the
        tilt-series level is an early but meaningful preprocessing step.
        It can substantially improve the visibility of weak structural
        features and provide a stronger basis for reconstruction and
        interpretation. Careful matching between imaging data and optical
        characterization remains the key requirement for obtaining
        biologically reliable results.
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
