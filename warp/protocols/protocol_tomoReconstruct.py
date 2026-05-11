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
    Reconstructs tomograms from tilt-series using WarpTools after CTF
    estimation and alignment import. Optionally produces half-tomograms
    for denoiser training and deconvolved reconstructions.
    More info:
        https://warpem.github.io/warp/user_guide/warptools/quick_start_warptools_tilt_series/#tilt-series-ctf-estimation
        https://warpem.github.io/warp/user_guide/warptools/quick_start_warptools_tilt_series/#tilt-series-reconstruct-tomograms

    AI Generated:

    Tomo Reconstruction (ProtWarpTomoReconstruct) — User Manual
        Overview

        The Tomo Reconstruction protocol processes a set of tilt-series
        into reconstructed tomograms using the Warp tilt-series workflow.
        It combines CTF handling, alignment import, and 3D reconstruction
        into a sequential pipeline suitable for cryo-electron tomography.

        The protocol is intended for users who already have aligned
        tilt-series and associated CTF estimations. From these inputs,
        it generates tomograms that can later be used for visualization,
        subtomogram averaging, particle picking, template matching, or
        denoising workflows.

        Biological Purpose

        In cryo-ET, tomogram reconstruction transforms a stack of
        2D tilted projections into a 3D representation of the specimen.
        This step is fundamental because all downstream biological
        interpretation depends directly on the quality of the resulting
        volume.

        Accurate reconstruction requires both correct CTF information
        and reliable alignment parameters. Errors introduced here often
        propagate to later stages such as segmentation, particle
        extraction, and subtomogram refinement.

        Inputs

        The protocol requires two main inputs:

        1. A set of tilt-series.
           These define the projection images and geometric acquisition
           parameters.

        2. A set of CTF estimations.
           These are used to update Warp-compatible CTF files before
           reconstruction.

        The protocol assumes that the tilt-series are valid, accessible,
        and that each series has matching metadata.

        Reconstruction Workflow

        For every enabled tilt-series, the protocol executes the
        following sequence:

        1. Create Warp settings for the individual tilt-series.
        2. Prepare tilt images and internal metadata.
        3. Generate or update CTF estimation files.
        4. Import alignment information from IMOD-compatible files.
        5. Launch Warp tomogram reconstruction.
        6. Register reconstructed tomograms as Scipion outputs.
        7. Remove temporary intermediate image files.

        Each tilt-series is processed independently, which makes the
        workflow robust for batch processing and easier to debug.

        Reconstruction Parameters

        Pixel size (angpix)
            Defines the voxel size of the reconstructed tomogram.
            This directly determines the sampling of the final 3D map.

        Tomogram thickness
            Defines the Z-size of the reconstruction volume in
            unbinned pixels.

        X and Y dimensions
            Optional advanced controls for the lateral size of the
            reconstruction. If omitted, the original tilt-series
            dimensions are used.

        These parameters should be chosen according to specimen size,
        expected biological context, and available computational memory.

        Half-Tomograms

        The protocol can optionally reconstruct two half-tomograms,
        each generated from half of the tilts.

        This is especially useful for:

        - denoiser training,
        - validation workflows,
        - consistency analysis between independent halves.

        When enabled, both half volumes are stored together with the
        main tomogram.

        Deconvolution and Contrast Options

        Deconvolution
            Produces a deconvolved tomogram that may improve visibility
            of structural features.

        Invert contrast
            Controls whether density contrast is inverted.
            For cryo-data, disabling inversion is often appropriate
            unless template matching requires opposite contrast.

        Normalize tilt images
            Controls whether input projections are normalized before
            reconstruction.

        These options affect interpretability and should be chosen
        according to the intended downstream application.

        CTF Handling

        Before reconstruction, the protocol generates Warp-compatible
        CTF estimation files and updates them using the provided
        CTF metadata.

        This ensures that the reconstruction uses corrected optical
        parameters while preserving previously estimated defocus
        information.

        Alignment Import

        Alignment parameters are written as IMOD-compatible files and
        imported into Warp before reconstruction.

        The imported alignment is rescaled according to the requested
        reconstruction pixel size, ensuring geometric consistency
        between tilt images and reconstructed tomograms.

        Outputs

        The protocol produces:

        - A SetOfTomograms containing one tomogram per input tilt-series.
        - Optional half-maps if half-tomogram generation is enabled.

        Each output tomogram contains:

        - file location,
        - sampling rate,
        - acquisition metadata,
        - default tomogram origin.

        The resulting tomograms are immediately suitable for further
        Scipion-based tomographic analysis.

        Output Registration

        Tomograms are appended incrementally to an output set.
        This allows progressive streaming-like behavior during
        processing of multiple tilt-series.

        Output dimensions are updated automatically after each new
        tomogram is created.

        Practical Recommendations

        For exploratory work:
            Use moderate pixel size and default normalization.

        For denoiser preparation:
            Enable half-tomogram generation.

        For template matching:
            Carefully check contrast inversion settings.

        For large specimens:
            Increase tomogram thickness to fully capture specimen depth.

        In practice, the most critical factor is consistency between
        acquisition metadata, CTF estimation, and alignment geometry.

        Cleanup Strategy

        After each tilt-series is processed, intermediate tilt-image
        files are removed automatically.

        This helps reduce storage usage, especially when processing
        large tomography datasets.

        Final Perspective

        Tomogram reconstruction is the bridge between raw tilt-series
        data and biologically interpretable 3D volumes.

        This protocol automates the main Warp reconstruction workflow
        while preserving compatibility with Scipion output objects,
        making it appropriate for both routine tomography processing
        and large-scale cryo-ET pipelines.
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
        argsDict = {
            "--settings": os.path.abspath(settingFile),
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
