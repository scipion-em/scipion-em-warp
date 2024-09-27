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
import glob
import os.path

from pwem.emlib.image.image_readers import ImageStack, ImageReadersRegistry
from pyworkflow import BETA
import pyworkflow.protocol.params as params
import tomo.objects as tomoObj
import pyworkflow.utils as pwutils
from tomo.protocols import ProtTomoBase
from warp.constants import CREATE_SETTINGS, TS_CTF, FS_MOTION_AND_CTF
from warp.protocols.protocol_base import ProtWarpBase
from warp.utils import parseCtfXMLFile


class ProtWarpTSCtfEstimation(ProtWarpBase, ProtTomoBase):
    """
    CTF estimation of a set of input tilt-series using the Warp procedure.
    More info:
        https://warpem.github.io/warp/user_guide/warptools/quick_start_warptools_tilt_series/#tilt-series-ctf-estimation
    """

    _label = 'CTF estimation'
    _devStatus = BETA

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputSet',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      label='Input set of tilt-series',
                      help='Input set of tilt-series')

        form.addParam('window', params.IntParam, default=512,
                      label='Windows', help='Patch size for CTF estimation in binned pixels')

        line = form.addLine('Resolution (Å)',
                            help='Resolution in Angstrom to consider in fit.')

        line.addParam('range_low', params.IntParam, default=30,
                      label='Min', help='Lowest (worst) resolution in Angstrom to consider in fit')

        line.addParam('range_high', params.IntParam, default=4,
                      label="Max",
                      help="Highest (best) resolution in Angstrom to consider in fit")

        line = form.addLine('Defocus search range (Å)',
                            help='Defocus values in um to explore during fitting (positive = underfocus)')
        line.addParam('defocus_min', params.FloatParam, default=0.5,
                      label='Min', help='Minimum defocus value in um to explore during fitting (positive = underfocus)')
        line.addParam('defocus_max', params.FloatParam, default=5,
                      label='Max', help='Maximum defocus value in um to explore during fitting (positive = underfocus)')

        form.addParam('voltage', params.IntParam, default=300,
                      label='Voltage', help=' Acceleration voltage of the microscope in kV')
        form.addParam('cs', params.FloatParam, default=2.7,
                      label=' Spherical aberration', help='Spherical aberration of the microscope in mm')

        form.addParam('amplitude', params.FloatParam, default=0.07,
                      label='Amplitude', help='Amplitude contrast of the sample, usually 0.07-0.10 for cryo')

        form.addParam('fit_phase', params.BooleanParam, default=True,
                      label='Fit phase', help='Fit the phase shift of a phase plate')

        form.addParam('auto_hand', params.IntParam, default=0,
                      label='Auto hand', help='Run defocus handedness estimation based on this many tilt series '
                                              '(e.g. 10), then estimate CTF with the correct handedness')

        form.addParam(params.GPU_LIST, params.StringParam, default='0',
                      label='Choose GPU IDs:', validators=[params.NonEmpty],
                      help='This argument is necessary. By default, the '
                           'protocol will attempt to launch on GPU 0. You can '
                           'override the default allocation by providing a '
                           'list of which GPUs (0,1,2,3, etc) to use. '
                           'GPU are separated by ",". For example: "0,1,5"')

    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.tsCtfEstimationStep)
        self._insertFunctionStep(self.createOutputStep)

    def convertInputStep(self):
        setOfTiltSeries = self.inputSet.get()
        starFolder = self._getExtraPath('tomostar')
        pwutils.makePath(starFolder)
        tiltimagesFolder = self._getExtraPath('tiltimages')
        pwutils.makePath(tiltimagesFolder)

        # 1. Extract all tilt images from the tiltseries
        for ts in setOfTiltSeries.iterItems():
            tsId = ts.getTsId()
            properties = {"sr": ts.getSamplingRate()}
            index = 1
            tiValues = {}
            for ti in ts.iterItems():
                newBinaryName = tsId + '_%s_%s.mrc' % (index, ti.getTiltAngle())
                newFrame = ImageStack(properties=properties)
                newFrame.append(ImageReadersRegistry.open(str(index) + '@' + ti.getFileName()))
                ImageReadersRegistry.write(newFrame, os.path.join(self._getExtraPath('tiltimages'), newBinaryName))
                # Taking angleTilt axisAngle Dose AverageIntensity MaskedFraction
                dose = 0
                if ts.hasAcquisition():
                    axisAngle = ts._acquisition.getTiltAxisAngle()
                else:
                    axisAngle = 0
                amplitudeContrast = 0
                maskedFraction = 0
                if ti.getAcquisition():
                    amplitudeContrast = ti.getAcquisition().getAmplitudeContrast()
                    dose = ti.getAcquisition().getDosePerFrame()

                tiValues[newBinaryName] = [ti.getTiltAngle(), axisAngle, dose, amplitudeContrast, maskedFraction]

                index += 1

            self.tomoStarGenerate(tsId, tiValues, starFolder)

        # 3. Create symbolic links from tiltseries folder to average folder
        #    We need to do this because warp needs both folders: The tiltimages folder to get
        #    the header of the images (it only needs the first tiltimage of each tiltseries),
        #    and the averages folder to read them.
        tiltSeriesFiles = glob.glob(os.path.join(tiltimagesFolder, '*'))
        averagesFolder = os.path.join(tiltimagesFolder, 'average')
        pwutils.makePath(averagesFolder)
        for file in tiltSeriesFiles:
            fileName = os.path.basename(file)
            destFolder = os.path.join(averagesFolder, fileName)
            os.symlink('../' + fileName, destFolder)

        # 4. Create warp_tiltseries.settings file
        self.processingFolder = os.path.abspath(self._getExtraPath('warp_tiltseries'))
        pwutils.makePath(self.processingFolder)
        argsDict = {
            "--folder_data": os.path.abspath(self._getExtraPath('tomostar')),
            "--extension": "*.tomostar",
            "--folder_processing":  self.processingFolder,
            '--angpix': 0.7894,
            '--exposure': 2.64,
            "--output": os.path.abspath(self._getExtraPath("warp_tiltseries.settings")),
        }

        cmd = ' '.join(['%s %s' % (k, v) for k, v in argsDict.items()])
        self.runJob(self.getPlugin().getProgram(CREATE_SETTINGS), cmd, executable='/bin/bash')

    def tomoStarGenerate(self, tsId, tiValues, otputFolder):
        """Generate the .tomostar files from TS"""
        _fileName = os.path.abspath(otputFolder) + '/%s.tomostar' % tsId
        _file = open(_fileName, 'a+')
        header = """
data_

loop_
_wrpMovieName #1
_wrpAngleTilt #2
_wrpAxisAngle #3
_wrpDose #4
_wrpAverageIntensity #5
_wrpMaskedFraction #6
"""
        _file.write(header)

        for key, value in tiValues.items():
            tiPath = '../tiltimages/' + key
            angleTilt = value[0]
            axisAngle = value[1]
            dose = value[2]
            averageIntensity = value[3]
            maskedFraction = value[4]
            _file.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (
                    tiPath, angleTilt, axisAngle, dose, averageIntensity, maskedFraction))

        _file.close()

    def motionAndCtfStep(self):
        """Estimate 2D sample motion and contrast transfer function"""
        argsDict = {
            "--settings": os.path.abspath(self._getExtraPath("warp_frameseries.settings")),
            "--m_grid": '1x1x3',
            "--c_grid": '2x2x1',
            "--c_range_max": self.range_high.get(),
            "--c_defocus_max": self.defocus_max.get(),
        }
        cmd = ' '.join(['%s %s' % (k, v) for k, v in argsDict.items()])
        cmd += ' --c_use_sum --out_averages --out_average_halves'

        self.runJob(self.getPlugin().getProgram(FS_MOTION_AND_CTF), cmd, executable='/bin/bash')

    def tsCtfEstimationStep(self):
        """CTF estimation"""
        argsDict = {
            "--settings": os.path.abspath(self._getExtraPath("warp_tiltseries.settings")),
            "--window": self.window.get(),
            "--range_low": self.range_low.get(),
            "--range_high": self.range_high.get(),
            "--defocus_min": self.defocus_min.get(),
            "--defocus_max": self.defocus_max.get(),
            "--voltage": self.voltage.get(),
            "--cs": self.cs.get(),
            "--amplitude": self.amplitude.get(),
            "--fit_phase": self.fit_phase.get(),
            "--auto_hand": self.auto_hand.get(),
        }
        cmd = ' '.join(['%s %s' % (k, v) for k, v in argsDict.items()])
        self.runJob(self.getPlugin().getProgram(TS_CTF), cmd, executable='/bin/bash')

    def createOutputStep(self):
        tsSet = self.inputSet.get()
        outputSetOfCTFTomoSeries = tomoObj.SetOfCTFTomoSeries.create(self._getPath(),
                                                                     template='CTFmodels%s.sqlite')
        outputSetOfCTFTomoSeries.setSetOfTiltSeries(tsSet)
        for ts in tsSet.iterItems():
            newCTFTomoSeries = tomoObj.CTFTomoSeries(tsId=ts.getTsId())
            newCTFTomoSeries.copyInfo(ts)
            newCTFTomoSeries.setTiltSeries(ts)
            outputSetOfCTFTomoSeries.append(newCTFTomoSeries)
            defocusFilePath = os.path.join(self.processingFolder, ts.getTsId() + '.xml')
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

        self._defineOutputs(**{'CTFTomoSeries': outputSetOfCTFTomoSeries})
        self._defineCtfRelation(outputSetOfCTFTomoSeries, tsSet)










