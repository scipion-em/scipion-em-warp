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
from tomo.protocols import ProtTomoBase
from warp.constants import TS_CTF
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

        form.addHidden(params.GPU_LIST, params.StringParam, default='0',
                      label='Choose GPU IDs:', validators=[params.NonEmpty],
                      help='This argument is necessary. By default, the '
                           'protocol will attempt to launch on GPU 0. You can '
                           'override the default allocation by providing a '
                           'list of which GPUs (0,1,2,3, etc) to use. '
                           'GPU are separated by ",". For example: "0,1,5"')

    def _insertAllSteps(self):
        self._insertFunctionStep(self.dataPrepare, self.inputSet.get())
        self._insertFunctionStep(self.tsCtfEstimationStep)
        self._insertFunctionStep(self.createOutputStep)

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
                # TODO We need to understand how to calculate the defocus values
                newCTFTomo.setDefocusU(gridCtfData["Nodes"][tiObjId] + defocusDelta)
                newCTFTomo.setDefocusV(gridCtfData["Nodes"][tiObjId] - defocusDelta)
                newCTFTomoSeries.append(newCTFTomo)

        self._defineOutputs(**{'CTFTomoSeries': outputSetOfCTFTomoSeries})
        self._defineCtfRelation(outputSetOfCTFTomoSeries, tsSet)










