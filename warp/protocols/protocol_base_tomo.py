# *
# * Authors:     Scipion Team
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
# *  e-mail address 'scipion-users@lists.sourceforge.net'
# *
# **************************************************************************
import logging
import math
import traceback
from abc import ABC, abstractmethod
from os.path import abspath, join, exists, normpath
from typing import Tuple, Union
from pwem.protocols import EMProtocol
from pyworkflow.utils import cyanStr, makePath, redStr
from tomo.objects import SetOfTiltSeriesM, SetOfTiltSeries
from warp import Plugin, WARP_TOOLS, CREATE_SETTINGS, TOMOSTAR_EXT, TILTSERIES_SETTINGS, FRAMESERIES_FOLDER, \
    TILTSERIES_FOLDER, TOMOSTAR_FOLDER, MATCHING_FOLDER, RELION_FOLDER, TILTSTACK_FOLDER
from warp.convert import writeTsStar

logger = logging.getLogger(__name__)


class ProtWarpTomoBase(EMProtocol, ABC):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.failedTsIds = []
        self.samplingRate = None
        self.exposure = None
        self.tsmDir = None
        self.tsDir = None
        self.tsStarDir = None
        self.seriesDict = None

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        inSet = self._getInSet()
        self.samplingRate = inSet.getSamplingRate()
        # self.outSamplingRate = self._getOutputSamplingRate()
        self.exposure = inSet.getAcquisition().getDosePerFrame()
        self.seriesDict = {series.getTsId(): series.clone() for series in inSet.iterItems()}
        # # Gain and defects
        # gainFile = tsMSet.getGain()
        # darkFile = tsMSet.getDark()
        # self.gainFile = abspath(gainFile) if gainFile else None
        # self.darkFile = abspath(darkFile) if darkFile else None
        # Make paths
        self.tsmDir = self._getFrameSeriesDir()
        self.tsDir = self._getTiltSeriesDir()
        self.tsStarDir = self._getTsStarDir()
        makePath(*[self.tsmDir, self.tsDir, self.tsStarDir])

    def createTsSettingsStep(self):
        logger.info(cyanStr(">>> Creating tilt-series settings..."))
        try:
            cmd = self._genCreateSettingsTsArgs()
            self.runJob(Plugin.getProgram(WARP_TOOLS, CREATE_SETTINGS), cmd, executable='/bin/bash')
        except Exception as e:
            raise Exception(f"{WARP_TOOLS} {CREATE_SETTINGS} failed with the exception --> {e}")

    def createTsMStarStep(self, tsId: str):
        try:
            logger.info(cyanStr(f">>> {tsId} - Creating the star file..."))
            writeTsStar(self.seriesDict[tsId], self._getOutTsStarFile(tsId))

        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to create the TS star file '
                                f'with exception {e}. Skipping... '))
            logger.error(traceback.format_exc())
            self.failedTsIds.append(tsId)

    # --------------------------- UTILS functions -----------------------------
    @abstractmethod
    def _getInSet(self) -> Union[SetOfTiltSeriesM, SetOfTiltSeries]:
        pass

    def _getTiltSeriesSettingsFn(self) -> str:
        return abspath(self._getExtraPath(TILTSERIES_SETTINGS))

    def _getFrameSeriesDir(self) -> str:
        return abspath(self._getExtraPath(FRAMESERIES_FOLDER))

    def _getTiltSeriesDir(self) -> str:
        return abspath(self._getExtraPath(TILTSERIES_FOLDER))

    def _getTsStarDir(self) -> str:
        return abspath(self._getExtraPath(TOMOSTAR_FOLDER))

    def _getOutTsStarFile(self, tsId: str) -> str:
        return join(self.tsStarDir, f'{tsId}{TOMOSTAR_EXT}')

    def _getTiltStackDir(self) -> str:
        return join(self._getTiltSeriesDir(), TILTSTACK_FOLDER)

    def _getTsTiltStackDir(self, tsId: str) -> str:
        return join(self._getTiltStackDir(), tsId)

    def _getParticlesFolder(self) -> str:
        return self._getExtraPath(MATCHING_FOLDER)

    def _getRelionFolder(self) -> str:
        return self._getExtraPath(RELION_FOLDER)

    def _genCreateSettingsTsArgs(self) -> str:
        x, y, z = self._getTomoDims()
        args = [
            f'--output {self._getTiltSeriesSettingsFn()}',
            f'--folder_processing {self.tsDir}',
            f'--folder_data {self.tsStarDir}',
            f'--extension *{TOMOSTAR_EXT}',
            f'--angpix {self.samplingRate}',
            f'--exposure {self.exposure}',
            f'--tomo_dimensions {x}x{y}x{z}'
        ]
        return ' '.join(args)

    def _getTomoDims(self) -> Tuple[int, int, int]:
        inSet = self._getInSet()
        x = int(inSet.getDimensions()[0])
        y = int(inSet.getDimensions()[1])
        # The binary swap is exact only near φ∈{0,90}. Near 45°, the detector rectangle rotated into
        # the tomogram frame needs a box larger than either W×H or H×W. The exact axis-aligned bounding
        # box of the rotated rectangle covers every imaged pixel at any angle and self-pads near 45°I
        tiltAxisAngle = inSet.getAcquisition().getTiltAxisAngle()
        c, s = abs(math.cos(tiltAxisAngle)), abs(math.sin(tiltAxisAngle))
        tomoX = math.ceil(x * c + y * s)
        tomoY = math.ceil(x * s + y * c)
        tomoZ = int(max(tomoX, tomoY) / 3)
        return tomoX, tomoY, tomoZ

    def normalizeParticlesPath(self, particleFile: str) -> str:
        if not exists(particleFile):
            particleFile = normpath(join(self._getTiltSeriesDir(), particleFile))
        return particleFile
    
    def normalizeTomogramsPath(self, tomogramFile: str) -> str:
        if not exists(tomogramFile):
            tomogramFile = join(self._getRelionFolder(), tomogramFile)
        return tomogramFile

    def normalizeOptimizationPath(self, optimizationValue:str) -> str:
        optimizationFile = optimizationValue
        if not exists(optimizationFile):
            optimizationFile = join(self._getExtraPath(RELION_FOLDER), optimizationValue)
        return optimizationFile