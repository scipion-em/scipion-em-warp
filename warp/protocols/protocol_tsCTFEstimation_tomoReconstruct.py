# **************************************************************************
# *
# * Authors: Yunior C. Fonseca Reyna    (cfonseca@cnb.csic.es)
# *
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
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
import tomo.objects as tomoObj
from pyworkflow.object import Set
from tomo.protocols import ProtTomoBase
from warp import TILTSERIE_SETTINGS, TILTSERIES_FOLDER
from warp.constants import TS_CTF, OUTPUT_CTF_SERIE, TS_RECONSTRUCTION, MRC_EXT, OUTPUT_TOMOGRAMS_NAME, \
    RECONSTRUCTION_FOLDER, RECONSTRUCTION_ODD_FOLDER, RECONSTRUCTION_EVEN_FOLDER
from warp.protocols.protocol_base import ProtWarpBase
from warp.utils import parseCtfXMLFile


class ProtWarpTSCtfEstimationTomoReconstruct(ProtWarpBase, ProtTomoBase):
    """
    CTF estimation of a set of input tilt-series and reconstruct tomograms for various tasks and, optionally,
    half-tomograms for denoiser training using the Warp procedure.
    More info:
        https://warpem.github.io/warp/user_guide/warptools/quick_start_warptools_tilt_series/#tilt-series-ctf-estimation
        https://warpem.github.io/warp/user_guide/warptools/quick_start_warptools_tilt_series/#tilt-series-reconstruct-tomograms
    """

    _label = 'CTF estimation and tomo reconstruction'
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

        form.addSection(label="CTF")
        form.addParam('window', params.IntParam, default=512,
                      label='Windows', help='Patch size for CTF estimation in binned pixels')

        line = form.addLine('Resolution (Å)',
                            help='Resolution in Angstrom to consider in fit.')

        line.addParam('range_low', params.FloatParam, default=30,
                      label='Min', help='Lowest (worst) resolution in Angstrom to consider in fit')

        line.addParam('range_high', params.FloatParam, default=4,
                      label="Max",
                      help="Highest (best) resolution in Angstrom to consider in fit")

        line = form.addLine('Defocus search range (Å)',
                            help='Defocus values in um to explore during fitting (positive = underfocus)')
        line.addParam('defocus_min', params.FloatParam, default=0.5,
                      label='Min', help='Minimum defocus value in um to explore during fitting (positive = underfocus)')
        line.addParam('defocus_max', params.FloatParam, default=5,
                      label='Max', help='Maximum defocus value in um to explore during fitting (positive = underfocus)')

        form.addParam('fit_phase', params.BooleanParam, default=True,
                      label='Fit phase', help='Fit the phase shift of a phase plate')

        form.addSection(label="Reconstruction")

        form.addParam('reconstruct', params.BooleanParam, default=False,
                      label='Tomogram reconstruction ?',
                      help='reconstruct tomograms for various tasks and, optionally'
                           'half-tomograms for denoiser training using the Warp procedure')

        form.addParam('angpix', params.IntParam, default=10,
                      condition='reconstruct==True',
                      label='Pixel size (Å)', help='Pixel size of the reconstructed tomograms in Angstrom')

        form.addParam('halfmap_tilts', params.BooleanParam, default=False,
                      condition='reconstruct==True',
                      label='Produce two half-tomograms?',
                      help='Produce two half-tomograms, each reconstructed from half of the tilts')

        form.addParam('deconv', params.BooleanParam, default=False,
                      condition='reconstruct==True',
                      label='Produce a deconvolved version',
                      help='Produce a deconvolved version; all half-tomograms, if requested, will also be deconvolved')

        form.addParam('dont_invert', params.BooleanParam, default=False,
                      condition='reconstruct==True',
                      label='Invert contrast?',
                      help='Invert the contrast; contrast inversion is needed for template matching on cryo '
                           'data, i.e. when the density is dark in original images')

        form.addParam('dont_normalize', params.BooleanParam, default=False,
                      condition='reconstruct==True',
                      label='Normalize the tilt images?',
                      help='Normalize the tilt images')

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
        if self.reconstruct.get():
            self._insertFunctionStep(self.tomoReconstructionStep)
        self._insertFunctionStep(self.createOutputStep)

    def tsCtfEstimationStep(self):
        """CTF estimation"""
        self.info(">>> Starting ctf estimation...")
        inputTSAdquisition = self.inputSet.get().getAcquisition()
        argsDict = {
            "--settings": os.path.abspath(self._getExtraPath(TILTSERIE_SETTINGS)),
            "--window": self.window.get(),
            "--range_low": self.range_low.get(),
            "--range_high": self.range_high.get(),
            "--defocus_min": self.defocus_min.get(),
            "--defocus_max": self.defocus_max.get(),
            "--voltage": int(inputTSAdquisition.getVoltage()),
            "--cs": inputTSAdquisition.getSphericalAberration(),
            "--amplitude": inputTSAdquisition.getAmplitudeContrast(),
            "--fit_phase": self.fit_phase.get(),
            "--auto_hand": 0,
        }
        self.runProgram(argsDict, TS_CTF)

    def tomoReconstructionStep(self):
        """Tomo Reconstruction"""
        self.info(">>> Starting tomogram reconstruction...")
        argsDict = {
            "--settings": os.path.abspath(self._getExtraPath(TILTSERIE_SETTINGS)),
            "--angpix": self.angpix.get(),
            "--halfmap_tilts": self.halfmap_tilts.get(),
            "--deconv": self.deconv.get(),
            "--dont_invert": not self.dont_invert.get(),
            "--dont_normalize": not self.dont_normalize.get(),
        }
        self.runProgram(argsDict, TS_RECONSTRUCTION)

    def createOutputStep(self):
        self.info(">>> Generating outputs...")
        processingFolder = os.path.abspath(self._getExtraPath(TILTSERIES_FOLDER))
        tomogramFolder = os.path.join(processingFolder, RECONSTRUCTION_FOLDER)
        generateOutput = False
        tsSet = self.inputSet.get()

        for ts in tsSet.iterItems():
            tsId = ts.getTsId()
            outputSetOfCTFTomoSeries = self.getOutputSetOfCTFTomoSeries(OUTPUT_CTF_SERIE)

            # CTF outputs
            newCTFTomoSeries = tomoObj.CTFTomoSeries(tsId=tsId)
            newCTFTomoSeries.copyInfo(ts)
            newCTFTomoSeries.setTiltSeries(ts)
            outputSetOfCTFTomoSeries.append(newCTFTomoSeries)
            defocusFilePath = os.path.join(processingFolder, ts.getTsId() + '.xml')
            ctfData, gridCtfData = parseCtfXMLFile(defocusFilePath)
            defocusDelta = float(ctfData['DefocusDelta']) * 1e4

            for ti in ts.iterItems():
                tiObjId = ti.getObjId()
                newCTFTomo = tomoObj.CTFTomo()
                newCTFTomo.setAcquisitionOrder(ti.getAcquisitionOrder())
                newCTFTomo.setIndex(ti.getIndex())
                newCTFTomo.setObjId(tiObjId)
                newCTFTomo.setDefocusU(gridCtfData["Nodes"][tiObjId] + defocusDelta)
                newCTFTomo.setDefocusV(gridCtfData["Nodes"][tiObjId] - defocusDelta)
                newCTFTomoSeries.append(newCTFTomo)

            outputSetOfCTFTomoSeries.update(newCTFTomoSeries)
            outputSetOfCTFTomoSeries.write()
            self._store(outputSetOfCTFTomoSeries)

            # Reconstruction outputs
            if self.reconstruct.get():
                tomoLocation = os.path.join(tomogramFolder, self.getOutFile(tsId, ext=MRC_EXT))
                if os.path.exists(tomoLocation):
                    generateOutput = True

                if generateOutput:
                    outputSetOfTomograms = self.getOutputSetOfTomograms(OUTPUT_TOMOGRAMS_NAME)
                    outputSetOfTomograms.setSamplingRate(self.angpix.get())
                    newTomogram = tomoObj.Tomogram(tsId=tsId)
                    newTomogram.copyInfo(ts)
                    newTomogram.setSamplingRate(self.angpix.get())
                    newTomogram.setLocation(tomoLocation)

                    if self.halfmap_tilts.get():
                        halfMapsList = [os.path.join(tomogramFolder, RECONSTRUCTION_ODD_FOLDER, self.getOutFile(tsId, ext=MRC_EXT)),
                                        os.path.join(tomogramFolder, RECONSTRUCTION_EVEN_FOLDER, self.getOutFile(tsId, ext=MRC_EXT))]
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

    def getOutputSetOfCTFTomoSeries(self, outputSetName):
        outputSetOfCTFTomoSeries = getattr(self, outputSetName, None)

        if outputSetOfCTFTomoSeries:
            outputSetOfCTFTomoSeries.enableAppend()
        else:
            outputSetOfCTFTomoSeries = tomoObj.SetOfCTFTomoSeries.create(self._getPath(),
                                                                  template='CTFmodels%s.sqlite')
            tsSet = self.inputSet.get()
            outputSetOfCTFTomoSeries.setSetOfTiltSeries(tsSet)
            outputSetOfCTFTomoSeries.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(**{outputSetName: outputSetOfCTFTomoSeries})
            self._defineCtfRelation(outputSetOfCTFTomoSeries, tsSet)

        return outputSetOfCTFTomoSeries

    def getOutputSetOfTomograms(self, outputSetName):
        outputSetOfTomograms = getattr(self, outputSetName, None)
        if outputSetOfTomograms:
            outputSetOfTomograms.enableAppend()
        else:
            tsSet = self.inputSet.get()
            outputSetOfTomograms = tomoObj.SetOfTomograms.create(self._getExtraPath(), template='tomograms%s.sqlite')
            outputSetOfTomograms.setAcquisition(tsSet.getAcquisition())
            outputSetOfTomograms.setSamplingRate(tsSet.getSamplingRate())

            self._defineOutputs(**{outputSetName: outputSetOfTomograms})
            self._defineSourceRelation(tsSet, outputSetOfTomograms)

        return outputSetOfTomograms

    def _summary(self):
        summary = []
        if self.hasAttribute('CTFTomoSeries'):
            summary.append(f"Input tilt-series: {self.inputSet.get().getSize()}\n"
                           f"CTF Estimation: {self.CTFTomoSeries.getSize()}")
        else:
            summary.append("Outputs are not ready yet.")
        return summary

    def _validate(self):
        rangeHigh = self.range_high.get()
        if rangeHigh is None or rangeHigh < self.inputSet.get().getSamplingRate() * 2:
            return ["Resolution parameter(Max) can't be higher than the binned data's Nyquist resolution"]

    def getOutFile(self, tsId, ext):
        angpix = self.angpix.get()
        suffix = str(f"{angpix:.2f}") + 'Apx'
        return f'{tsId}_{suffix}.{ext}'









