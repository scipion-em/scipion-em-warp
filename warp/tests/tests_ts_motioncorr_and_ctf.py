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
from typing import Optional, Tuple, List
from pyworkflow.tests import setupTestProject, DataSet
from pyworkflow.utils import cyanStr, magentaStr
from tomo.objects import SetOfCTFTomoSeries, SetOfTiltSeries
from tomo.protocols import ProtImportTsMovies
from tomo.tests import EMPIAR_10491, DataSetEmpiar10491, TS_1, TS_32
from tomo.tests.test_base_centralized_layer import TestBaseCentralizedLayer
from warp.protocols import ProtWarpTSMotionCorr
from warp.protocols.protocol_ts_motion_ctf import BINNING_FACTOR, FLIP_Y, TARGET_SAMPLING_RATE


class TestWarpMCandCTF(TestBaseCentralizedLayer):
    importedTsM = None

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet(EMPIAR_10491)
        cls.testAcq = DataSetEmpiar10491.acq.value
        cls.unbinnedApix = DataSetEmpiar10491.unbinnedApix.value
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
    def getWarpExpectedTsDimsTargetSr(cls,
                                      inDims: List[int],
                                      targetSr: float) -> List[int]:
        downsample = targetSr / cls.unbinnedApix  # == 2**bintimes
        binnedDims = [(round(d / downsample) // 2) * 2 for d in inDims[:2]]
        binnedDims.append(inDims[-1])
        return binnedDims

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
                                binningtargetVal: float = 1.,
                                ctfRangeMax: float = 4.) \
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
                               range_max=ctfRangeMax
                               )
        ctfLabel = ' CTF' if estimateCtf else ''
        oeLabel = ' OE' if doOddEven else ''
        label = f'MC{ctfLabel}{oeLabel}{binModeLabel}'
        prot.setObjLabel(label)
        cls.launchProtocol(prot)
        return (getattr(prot, prot._possibleOutputs.tiltSeries.name, None),
                getattr(prot, prot._possibleOutputs.ctfs.name, None))

    def _checkTsSet(self,
                    tsSet: SetOfTiltSeries,
                    estimateCtf: bool,
                    testSRate: float,
                    hasOddEven: bool,
                    targetSRateMode: bool) -> None:
        testNAngles = DataSetEmpiar10491.nAngles.value
        expectedDims = DataSetEmpiar10491.tsDims.value
        if targetSRateMode:
            expectedDims = self.getWarpExpectedTsDimsTargetSr(expectedDims, testSRate)
        self.checkTiltSeries(tsSet,
                             expectedSetSize=DataSetEmpiar10491.numTs.value,
                             expectedSRate=testSRate,
                             imported=True,
                             expectedDimensions=expectedDims,
                             testSetAcqObj=self.testAcq,
                             testAcqObj=self.testAcq,
                             hasOddEven=hasOddEven,
                             anglesCount=testNAngles,
                             anglesCountSet=testNAngles,
                             isHeterogeneousSet=False,
                             presentTsIds=[TS_1, TS_32],
                             defocusHandedness=1 if estimateCtf else None,
                             tiltAnglesTolDeg=0.02)

    def _checkCtfSet(self,
                     ctfSet: SetOfCTFTomoSeries):
        self.checkCTFs(ctfSet,
                       expectedSetSize=DataSetEmpiar10491.numTs.value,
                       expectedPsdFile=True)

    def _checkResults(self,
                      tsSet: Optional[SetOfTiltSeries],
                      ctfSet: Optional[SetOfCTFTomoSeries],
                      estimateCtf: bool,
                      testSRate: float = -1.,
                      hasOddEven: bool = False,
                      targetSRateMode: bool = False) -> None:
        if isinstance(tsSet, SetOfTiltSeries):
            self._checkTsSet(tsSet,
                             estimateCtf=estimateCtf,
                             testSRate=testSRate,
                             hasOddEven=hasOddEven,
                             targetSRateMode=targetSRateMode)
        else:
            raise Exception(f'Output SetOfTiltSeries is of type {tsSet}')

        if estimateCtf and ctfSet is None:
            raise Exception(f'Output SetOfCTF is unexpectedly empty.')
        elif isinstance(ctfSet, SetOfCTFTomoSeries):
            self._checkCtfSet(ctfSet)

    def test_01(self):
        estimateCtf = True
        doOddEven = False
        tsSet, ctfSet = self._runMotioncorrAndCtfEst(
            estimateCtf=estimateCtf,
            doOddEven=doOddEven,
            binningFactorMode=BINNING_FACTOR,
            binningtargetVal=1
        )
        self._checkResults(tsSet,
                           ctfSet,
                           estimateCtf=estimateCtf,
                           testSRate=self.unbinnedApix,
                           hasOddEven=doOddEven)

    def test_02(self):
        estimateCtf = False
        doOddEven = True
        binningTargetVal = 4
        tsSet, ctfSet = self._runMotioncorrAndCtfEst(
            estimateCtf=estimateCtf,
            doOddEven=doOddEven,
            binningFactorMode=TARGET_SAMPLING_RATE,
            binningtargetVal=binningTargetVal
        )
        self._checkResults(tsSet,
                           ctfSet,
                           estimateCtf=estimateCtf,
                           testSRate=binningTargetVal,
                           hasOddEven=doOddEven,
                           targetSRateMode=True)

    def test_03(self):
        estimateCtf = True
        doOddEven = True
        binningTargetVal = 3.12
        tsSet, ctfSet = self._runMotioncorrAndCtfEst(
            estimateCtf=estimateCtf,
            doOddEven=doOddEven,
            binningFactorMode=TARGET_SAMPLING_RATE,
            binningtargetVal=binningTargetVal,
            ctfRangeMax=8
        )
        self._checkResults(tsSet,
                           ctfSet,
                           estimateCtf=estimateCtf,
                           testSRate=binningTargetVal,
                           hasOddEven=doOddEven,
                           targetSRateMode=True)