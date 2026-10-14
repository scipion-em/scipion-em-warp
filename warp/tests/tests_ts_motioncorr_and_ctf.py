# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Scipion Team
# *
# * National Center of Biotechnology, CSIC, Spain
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
from typing import Optional, Tuple
from pyworkflow.tests import setupTestProject, DataSet
from pyworkflow.utils import cyanStr, magentaStr
from tomo.objects import SetOfCTFTomoSeries, SetOfTiltSeries
from tomo.protocols import ProtImportTsMovies
from tomo.tests import EMPIAR_10491, DataSetEmpiar10491, TS_1, TS_32
from tomo.tests.test_base_centralized_layer import TestBaseCentralizedLayer
from warp.protocols import ProtWarpTSMotionCorr
from warp.protocols.protocol_ts_motion_ctf import BINNING_FACTOR, FLIP_Y


class TestWarpMCandCTF(TestBaseCentralizedLayer):
    importedTsM = None

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet(EMPIAR_10491)
        cls.testAcq = DataSetEmpiar10491.acq.value
        cls.unbinnedApix = DataSetEmpiar10491.unbinnedApix
        cls.runPrevProtocols()

    @classmethod
    def runPrevProtocols(cls):
        print(cyanStr('--------------------------------- RUNNING PREVIOUS PROTOCOLS ---------------------------------'))
        cls._runPreviousProtocols()
        print(
            cyanStr('\n-------------------------------- PREVIOUS PROTOCOLS FINISHED ---------------------------------'))

    @classmethod
    def _runPreviousProtocols(cls):
        cls.importedTsM = cls._runImportTiltSeriesM()

    @classmethod
    def _runImportTiltSeriesM(cls):
        print(magentaStr("\n==> Importing data - tilt-series movies:"))
        protImport = cls.newProtocol(ProtImportTsMovies,
                                     filesPath=cls.dataset.getFile(DataSetEmpiar10491.tsMDir.value),
                                     filesPattern='*.mdoc',
                                     dosePerFrame=cls.testAcq.getDosePerFrame(),
                                     gainFile=cls.dataset.getFile(DataSetEmpiar10491.gainFile.value)
                                     )
        cls.launchProtocol(protImport)
        return getattr(protImport, protImport.OUTPUT_NAME, None)

    @classmethod
    def _runMotioncorrAndCtfEst(cls,
                                estimateCtf: bool = True,
                                doOddEven: bool = False,
                                binningFactorMode: int = BINNING_FACTOR,
                                binningtargetVal: int = 2) \
            -> Tuple[Optional[SetOfTiltSeries], Optional[SetOfCTFTomoSeries]]:

        if binningFactorMode == BINNING_FACTOR:
            binningModeStr = 'Binning factor'
            binModeLabel = ' Bin'
        else:
            binningModeStr = 'Target sampling rate (Angst/pix)'
            binModeLabel = ' Sr'

        print(magentaStr("\n==> Running the tilt-series motion correction:"
                         f"\n\t- Estimate CTF = {estimateCtf}"
                         f"\n\t- Split odd/even = {doOddEven}"
                         f"\n\t- Binning mode = {binningModeStr}"
                         f"\n\t- Binning target value = {binningtargetVal}"))

        prot = cls.newProtocol(ProtWarpTSMotionCorr,
                               inputTSMovies=cls.importedTsM,
                               binFactorMode=binningFactorMode,
                               binFactor=binningtargetVal,
                               binTarget=binningtargetVal,
                               average_halves=doOddEven,
                               gainFlip=FLIP_Y,
                               estimateCTF=estimateCtf,
                               )
        ctfLabel = ' CTF' if estimateCtf else ''
        oeLabel = ' OE' if doOddEven else ''
        label = f'MC{ctfLabel}{oeLabel}{binModeLabel}'
        prot.setObjLabel(label)
        cls.launchProtocol(prot)
        return (getattr(prot, prot._possibleOutputs.tiltSeries.value, None),
                getattr(prot, prot._possibleOutputs.ctfs.value, None))

    def _checkTsSet(self,
                    tsSet: SetOfTiltSeries,
                    testSRate: float,
                    hasOddEven: bool) -> None:
        testNAngles = DataSetEmpiar10491.nAngles.value
        self.checkTiltSeries(tsSet,
                             expectedSetSize=DataSetEmpiar10491.numTs.value,
                             expectedSRate=testSRate,
                             imported=True,
                             expectedDimensions=DataSetEmpiar10491.tsDims.value,
                             testSetAcqObj=self.testAcq,
                             testAcqObj=self.testAcq,
                             hasOddEven=hasOddEven,
                             anglesCount=testNAngles,
                             anglesCountSet=testNAngles,
                             isHeterogeneousSet=False,
                             presentTsIds=[TS_1, TS_32])

    def _checkCtfSet(self,
                     ctfSet: SetOfCTFTomoSeries, ):
        self.checkCTFs(ctfSet,
                       expectedSetSize=DataSetEmpiar10491.numTs.value,
                       expectedPsdFile=True,
                       handedness=1)
