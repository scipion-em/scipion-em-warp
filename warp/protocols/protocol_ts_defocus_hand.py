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
from pyworkflow.object import  String
from tomo.protocols import ProtTomoBase
from warp.constants import FS_CTF,  TS_DEFOCUS_HAND
from warp.protocols.protocol_base import ProtWarpBase


class ProtWarpTSDefocusHand(ProtWarpBase, ProtTomoBase):
    """
    Check and/or set defocus handedness for all tilt series
    More info:
        https://warpem.github.io/warp/user_guide/warptools/quick_start_warptools_tilt_series/#tilt-series-check-defocus-handedness
    """

    _label = 'CTF defocus hand'
    _devStatus = BETA

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputSet',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      label='Input set of tilt-series',
                      help='Input set of tilt-series')

    def _insertAllSteps(self):
        self._insertFunctionStep(self.dataPrepare, self.inputSet.get())
        self._insertFunctionStep(self.tsDefocusEstimationStep)
        self._insertFunctionStep(self.tsDefocusHandStep)
        self._insertFunctionStep(self.createOutputStep)

    def tsDefocusEstimationStep(self):
        """Estimate CTF parameters in frame series(averages)"""
        self.info(">>> Starting ctf estimation...")
        inputTS = self.inputSet.get()
        averagesPath = os.path.abspath(self._getExtraPath('tiltimages'))
        argsDict = {
            "--settings": os.path.abspath(self._getExtraPath("warp_tiltseries.settings")),
            "--input_data": os.path.join(averagesPath, '*.mrc'),
            "--grid": "2x2x1",
            "--output_processing": averagesPath,
            "--range_max": inputTS.getSamplingRate() * 2 + 0.1
        }
        cmd = ' '.join(['%s %s' % (k, v) for k, v in argsDict.items()])
        self.runJob(self.getPlugin().getProgram(FS_CTF), cmd, executable='/bin/bash')

    def tsDefocusHandStep(self):
        """Defocus handedness"""
        self.info(">>> Starting defocus handedness...")
        argsDict = {
            "--settings": os.path.abspath(self._getExtraPath("warp_tiltseries.settings")),
            '--output_processing': os.path.abspath(self._getExtraPath('warp_tiltseries'))
        }
        cmd = ' '.join(['%s %s' % (k, v) for k, v in argsDict.items()])
        cmd += ' --check'
        self.runJob(self.getPlugin().getProgram(TS_DEFOCUS_HAND), cmd, executable='/bin/bash')

    def createOutputStep(self):
        stdoutFile = os.path.abspath(os.path.join(self.getPath(), 'logs', 'run.stdout'))
        with open(stdoutFile, 'r', encoding='utf-8') as file:
            lines = file.readlines()

        for line in reversed(lines):
            if 'Average correlation:' in line:
                average = line.split()[-1]
                outputAverage = String(average + " (defocus handedness should be set to 'no flip')")
                if float(average) < 0:
                    outputAverage = String(average + " (defocus handedness should be set to 'flip')")

                self._defineOutputs(**{'Average correlation': outputAverage})
                break



