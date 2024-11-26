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
import tomo.objects as tomoObj
from tomo.protocols import ProtTomoBase

from warp.constants import (TILTSERIE_SETTINGS, TILTSERIES_FOLDER, TS_CTF,
                            OUTPUT_CTF_SERIE, TS_RECONSTRUCTION, MRC_EXT, OUTPUT_TOMOGRAMS_NAME,
                            RECONSTRUCTION_FOLDER, RECONSTRUCTION_ODD_FOLDER, RECONSTRUCTION_EVEN_FOLDER)
from warp.protocols.protocol_base import ProtWarpBase


class ProtWarpTomoReconstruct(ProtWarpBase, ProtTomoBase):
    """
    CTF estimation of a set of input tilt-series and reconstruct tomograms for various tasks and, optionally,
    half-tomograms for denoiser training using the Warp procedure.
    More info:
        https://warpem.github.io/warp/user_guide/warptools/quick_start_warptools_tilt_series/#tilt-series-ctf-estimation
        https://warpem.github.io/warp/user_guide/warptools/quick_start_warptools_tilt_series/#tilt-series-reconstruct-tomograms
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

        form.addSection(label="Reconstruction")

        form.addParam('binFactor', params.IntParam,
                      default=4, label='Binning factor', important=True,
                      help='Binning factor of the reconstructed tomograms')

        # form.addParam('angpix', params.IntParam, default=10,
        #               condition='reconstruct==True',
        #               label='Pixel size (Å)', help='Pixel size of the reconstructed tomograms in Angstrom')

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
        self._insertFunctionStep(self.dataPrepare, self.inputSet.get())
        self._insertFunctionStep(self.tsCtfEstimationStep)
        self._insertFunctionStep(self.tomoReconstructionStep)
        self._insertFunctionStep(self.createOutputStep)

    def tsCtfEstimationStep(self):
        """CTF estimation"""
        self.info(">>> Starting ctf estimation...")
        inputTSAdquisition = self.inputSet.get().getFirstItem().getAcquisition()
        argsDict = {
            "--settings": os.path.abspath(self._getExtraPath(TILTSERIE_SETTINGS)),
            "--window": 512,
            "--range_low": 30,
            "--range_high": 8,
            # "--range_high": self.inputSet.get().getSamplingRate() * 2 + 0.1,
            "--defocus_min": 0.5,
            "--defocus_max": 5,
            "--voltage": int(inputTSAdquisition.getVoltage()),
            "--cs": inputTSAdquisition.getSphericalAberration(),
            "--amplitude": inputTSAdquisition.getAmplitudeContrast(),
        }
        self.runProgram(argsDict, TS_CTF)

    def tomoReconstructionStep(self):
        """Tomo Reconstruction"""
        self.info(">>> Starting tomogram reconstruction...")
        angpix = self.getAngPix()
        argsDict = {
            "--settings": os.path.abspath(self._getExtraPath(TILTSERIE_SETTINGS)),
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

        self.runProgram(argsDict, TS_RECONSTRUCTION, othersCmds=cmd)

    def createOutputStep(self):
        self.info(">>> Generating outputs...")
        processingFolder = os.path.abspath(self._getExtraPath(TILTSERIES_FOLDER))
        tomogramFolder = os.path.join(processingFolder, RECONSTRUCTION_FOLDER)
        generateOutput = False
        tsSet = self.inputSet.get()

        for ts in tsSet.iterItems():
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
                outputSetOfTomograms.append(newTomogram)
                outputSetOfTomograms.updateDim()
                outputSetOfTomograms.update(newTomogram)
                outputSetOfTomograms.write()
                self._store(outputSetOfTomograms)

            else:
                self.error(">>> Some error occurred in the reconstruction process. Please go to the "
                           "process log(.../extra/warp_tiltseries/logs)")

        self._closeOutputSet()

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
        return round(self.inputSet.get().getSamplingRate() * self.binFactor.get())

    def getOutFile(self, tsId, ext) -> str:
        angpix = self.getAngPix()
        suffix = str(f"{angpix:.2f}") + 'Apx'
        return f'{tsId}_{suffix}.{ext}'

    def getBinFactor(self):
        import math
        return math.floor(math.log2(self.binFactor.get()))
