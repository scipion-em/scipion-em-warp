# **************************************************************************
# *
# * Authors: Yunior C. Fonseca Reyna    (cfonseca@cnb.csic.es)
# *
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

import os.path
from pyworkflow import BETA
import pyworkflow.protocol.params as params
from pyworkflow.object import Set
import pyworkflow.utils as pwutils
import tomo.objects as tomoObj
from tomo.protocols import ProtTomoBase

from warp.constants import (TILTSERIE_SETTINGS, TILTSERIES_FOLDER, TS_CTF,
                            OUTPUT_CTF_SERIE, TS_RECONSTRUCTION, MRC_EXT, OUTPUT_TOMOGRAMS_NAME,
                            RECONSTRUCTION_FOLDER, RECONSTRUCTION_ODD_FOLDER, RECONSTRUCTION_EVEN_FOLDER,
                            TILTIMAGES_FOLDER, SETTINGS_FOLDER, TS_IMPORT_ALIGNMENTS, WARP_TOOLS)
from warp.protocols.protocol_base import ProtWarpBase
from warp.utils import updateCtFXMLFile


class ProtWarpTomoReconstruct(ProtWarpBase, ProtTomoBase):
    """
    Performs tomographic reconstruction from aligned cryo-electron tilt
    series, combining contrast transfer information with geometric
    alignment in order to generate tomograms suitable for downstream
    structural interpretation and subtomogram analysis.
    More info:
        https://warpem.github.io/warp/user_guide/warptools/quick_start_warptools_tilt_series/#tilt-series-ctf-estimation
        https://warpem.github.io/warp/user_guide/warptools/quick_start_warptools_tilt_series/#tilt-series-reconstruct-tomograms

    AI Generated:

    Tomogram Reconstruction (ProtWarpTomoReconstruct) - User Manual
        Overview

        The Tomogram Reconstruction protocol is designed to transform a
        collection of aligned tilt series into reconstructed tomographic
        volumes that can be directly used for biological interpretation.
        Its main objective is to convert two-dimensional tilt images into
        three-dimensional representations of the specimen while preserving
        the structural information required for cellular and molecular
        analysis.

        In cryo-electron tomography workflows, this stage is one of the
        most important because it creates the volumetric context in which
        macromolecular complexes, membranes, organelles, and intracellular
        organization become visible. The protocol is therefore especially
        relevant when the user intends to perform particle picking,
        subtomogram averaging, segmentation, or direct visual inspection
        of native biological environments.

        Inputs and Biological Context

        The protocol requires a set of tilt series that already contain
        geometric alignment information. This means that the images should
        represent a coherent angular acquisition of the same specimen area.
        In addition, a corresponding set of contrast transfer information
        can be supplied so that the reconstructed tomograms preserve more
        reliable structural detail across spatial frequencies.

        From a biological perspective, the quality of the resulting
        tomograms depends strongly on the quality of the input data.
        Stable acquisition, accurate tilt geometry, and well estimated
        optical parameters all contribute directly to the interpretability
        of the final three-dimensional volume.

        Reconstruction Strategy

        The protocol reconstructs one tomogram for each input tilt series.
        Each output volume represents the original specimen region in a
        form that can be explored in three dimensions. This allows users
        to examine cellular landscapes, locate molecular assemblies, and
        prepare regions of interest for more focused downstream analyses.

        The reconstructed voxel size can be chosen according to the
        biological objective. Coarser sampling is often appropriate for
        rapid inspection of large cellular regions, whereas finer sampling
        is generally preferred when the reconstructed volume will later be
        used for subtomogram extraction or template matching.

        Tomogram Dimensions and Spatial Coverage

        An important practical consideration is the size of the final
        reconstruction. The protocol allows the user to define the
        thickness of the tomogram as well as its lateral dimensions.

        Biologically, these parameters determine how much of the specimen
        is represented in the final volume. A reconstruction that is too
        small may truncate meaningful structural features, while a volume
        that is unnecessarily large may increase computational cost
        without adding useful information.

        When the specimen is relatively thin, moderate thickness values
        are often sufficient. In thicker cellular samples or lamellae,
        larger reconstruction depth may be necessary to preserve the full
        three-dimensional context.

        Half Tomograms and Validation Workflows

        The protocol can optionally generate two independent half
        tomograms reconstructed from complementary subsets of the tilt
        images. These half reconstructions are particularly valuable for
        validation and denoising workflows.

        In practical biological applications, half tomograms are useful
        when training denoising procedures that must avoid introducing
        artificial correlations. They also provide a principled way to
        evaluate whether observed structural features are reproducible
        rather than reconstruction artifacts.

        Contrast and Image Conditioning

        Several optional conditioning choices affect the biological
        appearance of the reconstructed density. Contrast inversion can be
        particularly relevant for template matching workflows, where the
        expected sign of density matters. Normalization helps produce more
        stable intensity behavior across tilt images, which often improves
        consistency in the final tomogram.

        The protocol can also generate a deconvolved version of the
        reconstruction. For biological users, deconvolution may improve
        visual sharpness and make boundaries or compact macromolecular
        features easier to recognize. However, deconvolution should always
        be interpreted with caution because enhanced contrast does not
        necessarily imply improved biological truth.

        Outputs and Interpretation

        The primary output is a set of tomograms, one for each input tilt
        series. Each reconstructed volume preserves the identity of its
        originating acquisition while making the sample accessible in
        three dimensions.

        When half reconstruction is enabled, each tomogram is accompanied
        by corresponding half volumes. These additional outputs can be
        important for denoiser training, reproducibility assessment, and
        careful interpretation of weak structural signals.

        For biological interpretation, the tomogram should be viewed as an
        experimentally constrained representation of the specimen rather
        than as a perfectly faithful model. Contrast variations, missing
        wedge effects, and local thickness differences remain important
        factors that can influence visibility of structural features.

        Practical Recommendations

        In routine cryo-electron tomography practice, it is usually best
        to begin with conservative reconstruction settings and inspect the
        resulting tomograms visually. If the goal is exploratory cellular
        analysis, moderate sampling and standard normalization are often
        sufficient.

        When the tomograms will feed into subtomogram averaging or
        particle extraction workflows, closer attention should be paid to
        sampling, contrast convention, and the consistency of the
        reconstructed dimensions. If denoising or validation is expected
        downstream, generating half tomograms from the start is generally
        a good strategy.

        Final Perspective

        Tomogram reconstruction is not merely a computational conversion
        of tilt images into a three-dimensional volume. It defines the
        structural landscape from which downstream biological conclusions
        will be drawn. Careful choice of reconstruction sampling, volume
        dimensions, and conditioning options directly affects the
        interpretability of macromolecular organization inside the native
        biological specimen.
    """

    _label = 'tomo reconstruction'
    _possibleOutputs = {OUTPUT_CTF_SERIE: tomoObj.SetOfCTFTomoSeries,
                        OUTPUT_TOMOGRAMS_NAME: tomoObj.SetOfTomograms}
    _devStatus = BETA

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputSet',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      label='Input set of tilt-series',
                      help='Input set of tilt-series')
        form.addParam('inputSetOfCtfTomoSeries',
                      params.PointerParam,
                      label="Input CTF estimation",
                      pointerClass='SetOfCTFTomoSeries',
                      help='Select the CTF estimation for the set '
                           'of tilt-series.')

        form.addSection(label="Reconstruction")

        form.addParam('angpix', params.IntParam, default=10,
                      label='Pixel size (Å)',
                      help='Pixel size of the reconstructed tomograms in Angstrom')

        form.addParam('halfmap_tilts', params.BooleanParam, default=False,
                      label='Produce two half-tomograms?',
                      help='Produce two half-tomograms, each reconstructed from half of the tilts')

        form.addParam('deconv', params.BooleanParam, default=False,
                      label='Produce a deconvolved version',
                      help='Produce a deconvolved version; all half-tomograms, if requested, will also be deconvolved')

        form.addParam('invert', params.BooleanParam, default=False,
                      label='Invert contrast?',
                      help='Invert the contrast; contrast inversion is needed for template matching on cryo '
                           'data, i.e. when the density is dark in original images')

        form.addParam('normalize', params.BooleanParam, default=True,
                      label='Normalize the tilt images?',
                      help='Normalize the tilt images')

        form.addParam('tomo_thickness', params.IntParam, default='1000',
                      important=True,
                      label='Tomogram thickness unbinned (pixels)',
                      help="Z height of the reconstructed volume in unbinned pixels.")

        form.addParam('x_dimension', params.IntParam, default=None,
                      allowsNull=True,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Tomogram x dimension unbinned (pixels)',
                      help="X width of the reconstructed volume in unbinned pixels. If the value is None or 0, "
                           "the dimension of the tiltseries will be taken into account. ")

        form.addParam('y_dimension', params.IntParam, default=None,
                      allowsNull=True,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Tomogram Y dimension unbinned (pixels)',
                      help="Y height of the reconstructed volume in unbinned pixels. If the value is None or 0, "
                           "the dimension of the tiltseries will be taken into account.")

        # form.addParam('tomo_dimensions', params.IntParam, default='1000',
        #               condition='reconstruct==True',
        #               label='Tomogram thickness unbinned (pixels)',
        #               help="Z height of the reconstructed volume in unbinned pixels.")

        """
       --deconv_strength         Default: 1. Strength of the deconvolution filter, if requested

       --deconv_falloff          Default: 1. Fall-off of the deconvolution filter, if requested

       --deconv_highpass         Default: 300. High-pass value (in Angstrom) of the deconvolution filter, if requested

       --keep_full_voxels        Mask out voxels that aren't contained in some of the tilt images (due to excessive sample shifts); don't use if you intend to run template matching

       --dont_mask               Don't apply a mask to each tilt image if available; otherwise, masked areas will be filled with Gaussian noise

       --dont_overwrite          Don't overwrite existing tomograms in output directory

       --subvolume_size          Default: 64. Reconstruction is performed locally using sub-volumes of this size in pixel

       --subvolume_padding       Default: 3. Padding factor for the reconstruction sub-volumes (helps with aliasing effects at sub-volume borders)

               """

        form.addHidden(params.GPU_LIST, params.StringParam, default='0',
                       label='Choose GPU IDs:', validators=[params.NonEmpty],
                       help="Space-separated list of GPU IDs to use for processing. Default: all GPUs in the system."
                            " Warp can use multiple GPUs - in that case"
                            " set to i.e. *0 1 2*.")

    def _insertAllSteps(self):
        inputTs = self.inputSet.get()
        for ts in inputTs.iterItems(iterate=False):
            if not ts.isEnabled():
                continue
            self._insertFunctionStep(self.tomoReconstructionStep, ts, needsGPU=True)
            self._insertFunctionStep(self.createOutput, ts, needsGPU=False)
            self._insertFunctionStep(self.cleanIntermediateResults, needsGPU=False)

        self._insertFunctionStep(self._closeOutputSet, needsGPU=False)

    def tsCtfEstimation(self, ts):
        """CTF estimation"""
        self.info(">>> Generating ctf estimation file fo %s ..." % ts.getTsId())
        settingFile = self._getExtraPath(SETTINGS_FOLDER, ts.getTsId() + '_' + TILTSERIE_SETTINGS)
        tsSr = ts.getSamplingRate()
        argsDict = {
            "--settings": os.path.abspath(settingFile),
            "--range_high": tsSr * 3,
            "--range_low": tsSr * 4,

        }
        try:
            self.runProgram(argsDict, WARP_TOOLS, TS_CTF)
        except Exception:
            self.info(">>> Error generating ctf estimation file...")
        ctfTomoSeries = self.inputSetOfCtfTomoSeries.get().getItem('_tsId', ts.getTsId())
        processingFolder = os.path.abspath(self._getExtraPath(TILTSERIES_FOLDER))
        defocusFilePath = os.path.join(processingFolder, ts.getTsId() + '.xml')
        updateCtFXMLFile(defocusFilePath, ctfTomoSeries)

    def tsImportAligments(self, ts):
        processingFolder = os.path.abspath(self._getExtraPath(TILTSERIES_FOLDER))
        tiltstackFolder = os.path.join(processingFolder, 'tiltstack', ts.getTsId())
        pwutils.makePath(tiltstackFolder)
        factor = self.angpix.get() / ts.getSamplingRate()
        ts.writeImodFiles(tiltstackFolder, delimiter=' ', factor=factor)
        self.info(">>> Starting import aligments...")
        settingFile = self._getExtraPath(SETTINGS_FOLDER, ts.getTsId() + '_' + TILTSERIE_SETTINGS)
        argsDict = {
            "--settings": os.path.abspath(settingFile),
            '--alignments': os.path.abspath(tiltstackFolder),
            "--alignment_angpix": self.angpix.get(),
        }
        cmd = ' '.join(['%s %s' % (k, v) for k, v in argsDict.items()])
        self.runJob(self.getPlugin().getProgram(WARP_TOOLS, TS_IMPORT_ALIGNMENTS), cmd, executable='/bin/bash')

    def tomoReconstructionStep(self, ts):
        """Tomo Reconstruction"""
        self.createTiltSeriesSetting(ts)
        self.tsDataPrepare(ts)
        self.tsCtfEstimation(ts)
        self.tsImportAligments(ts)
        self.info(">>> Starting tomogram reconstruction...")
        angpix = self.angpix.get()
        settingFile = self._getExtraPath(SETTINGS_FOLDER, ts.getTsId() + '_' + TILTSERIE_SETTINGS)
        argsDict = {
            "--settings": os.path.abspath(settingFile),
            "--angpix": angpix,
        }

        cmd = ''
        if self.halfmap_tilts.get():
            cmd += " --halfmap_tilts"
        if self.deconv.get():
            cmd += " --deconv"
        if not self.invert.get():
            cmd += " --dont_invert"
        if not self.normalize.get():
            cmd += " --dont_normalize"

        self.runProgram(argsDict, WARP_TOOLS, TS_RECONSTRUCTION, othersCmds=cmd)

    def createOutput(self, ts):
        self.info(">>> Generating outputs...")
        processingFolder = os.path.abspath(self._getExtraPath(TILTSERIES_FOLDER))
        tomogramFolder = os.path.join(processingFolder, RECONSTRUCTION_FOLDER)
        generateOutput = False
        tsId = ts.getTsId()
        tomoLocation = os.path.join(tomogramFolder, self.getOutFile(tsId, ext=MRC_EXT))
        if os.path.exists(tomoLocation):
            generateOutput = True

        if generateOutput:
            outputSetOfTomograms = self.getOutputSetOfTomograms(OUTPUT_TOMOGRAMS_NAME)
            outputSetOfTomograms.setSamplingRate(self.getAngPix())
            newTomogram = tomoObj.Tomogram(tsId=tsId)
            newTomogram.copyInfo(ts)
            newTomogram.setSamplingRate(self.getAngPix())
            newTomogram.setLocation(tomoLocation)

            if self.halfmap_tilts.get():
                halfMapsList = [os.path.join(tomogramFolder, RECONSTRUCTION_ODD_FOLDER,
                                self.getOutFile(tsId, ext=MRC_EXT)),
                                os.path.join(tomogramFolder, RECONSTRUCTION_EVEN_FOLDER,
                                self.getOutFile(tsId, ext=MRC_EXT))]
                newTomogram.setHalfMaps(halfMapsList)

            # Set default tomogram origin
            newTomogram.setOrigin(newOrigin=None)
            newTomogram.fixMRCVolume(True)
            outputSetOfTomograms.append(newTomogram)
            outputSetOfTomograms.updateDim()
            outputSetOfTomograms.update(newTomogram)
            outputSetOfTomograms.write()
            self._store(outputSetOfTomograms)

        else:
            self.error(">>> Some error occurred in the reconstruction process. Please go to the "
                       "process log(.../extra/warp_tiltseries/logs)")

    def getOutputSetOfTomograms(self, outputSetName):
        outputSetOfTomograms = getattr(self, outputSetName, None)
        if outputSetOfTomograms:
            outputSetOfTomograms.enableAppend()
        else:
            tsSet = self.inputSet.get()
            outputSetOfTomograms = tomoObj.SetOfTomograms.create(self._getExtraPath(), template='tomograms%s.sqlite')
            outputSetOfTomograms.setAcquisition(tsSet.getAcquisition())
            outputSetOfTomograms.setSamplingRate(tsSet.getSamplingRate())
            outputSetOfTomograms.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(**{outputSetName: outputSetOfTomograms})
            self._defineSourceRelation(tsSet, outputSetOfTomograms)

        return outputSetOfTomograms

    def _summary(self):
        summary = []
        if self.hasAttribute(OUTPUT_CTF_SERIE) and self.hasAttribute(OUTPUT_TOMOGRAMS_NAME):
            summary.append(f"Input tilt-series: {self.inputSet.get().getSize()}\n"
                           f"CTF Estimation: {self.CTFTomoSeries.getSize()}\n"
                           f"Tomograms: {self.Tomograms.getSize()}")
        else:
            summary.append("Outputs are not ready yet.")
        return summary

    def getAngPix(self):
        return self.angpix.get()

    def getOutFile(self, tsId, ext) -> str:
        angpix = self.getAngPix()
        suffix = str(f"{angpix:.2f}") + 'Apx'
        return f'{tsId}_{suffix}.{ext}'

    def cleanIntermediateResults(self):
        self.info(">>> Cleaning intermediate results...")
        imagesFolder = self._getExtraPath(TILTIMAGES_FOLDER)
        pwutils.cleanPath(imagesFolder)
