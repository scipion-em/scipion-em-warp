# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology (MRC-LMB)
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

from pwem import Domain
from pyworkflow.tests import BaseTest, setupTestProject
from pyworkflow.utils import magentaStr

from tomo.protocols import ProtImportTs, ProtImportTsCTF, ProtImportTsMovies
from tomo.tests import DataSet
from warp.protocols import (ProtWarpTomoReconstruct, ProtWarpDeconvTomo,
                            ProtWarpDeconvTS, ProtWarpTSMotionCorr)


class TestWarpBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.inputDataSet = DataSet.getDataSet('empiar10064')
        cls.ts_path = cls.inputDataSet.getPath()
        cls.inputTSM = DataSet.getDataSet('empiar_10491')
        cls.tsm_path = cls.inputTSM.getPath()

    @classmethod
    def runImportTiltSeries(cls, **kwargs):
        cls.protImportTS = cls.newProtocol(ProtImportTs, **kwargs)
        cls.launchProtocol(cls.protImportTS)
        cls.assertIsNotNone(cls.protImportTS.outputTiltSeries,
                            "SetOfTiltSeries has not been produced.")
        return cls.protImportTS

    @classmethod
    def runImportTiltSeriesM(cls, **kwargs):
        cls.protImportTSM = cls.newProtocol(ProtImportTsMovies, **kwargs)
        cls.launchProtocol(cls.protImportTSM)
        cls.assertIsNotNone(cls.protImportTSM.outputTiltSeriesM,
                            "SetOfTiltSeriesM has not been imported.")
        return cls.protImportTSM

    @classmethod
    def runMotioncorrTSMovieAligment(cls, **kwargs):
        cls.motioncorrTSMovieAligment = cls.newProtocol(ProtWarpTSMotionCorr, **kwargs)
        cls.launchProtocol(cls.motioncorrTSMovieAligment)
        cls.assertIsNotNone(cls.motioncorrTSMovieAligment.TiltSeries,
                            "SetOfTiltSeries has not been produced.")
        cls.assertIsNotNone(cls.motioncorrTSMovieAligment.CTFTomoSeries,
                            "CTFTomoSeries has not been produced.")
        return cls.motioncorrTSMovieAligment

    @classmethod
    def runImodImportTMatrix(cls, **kwargs):
        imodProts = Domain.importFromPlugin('imod.protocols', doRaise=True)
        cls.importTransformationMatrix = cls.newProtocol(imodProts.ProtImodImportTransformationMatrix, **kwargs)
        cls.launchProtocol(cls.importTransformationMatrix)
        cls.assertIsNotNone(cls.importTransformationMatrix.TiltSeries,
                            "Transformation matrix has not been produced.")
        return cls.importTransformationMatrix

    @classmethod
    def runWarpCTFEstimationTomoReconstruction(cls, **kwargs):
        cls.tomoReconstruct = cls.newProtocol(ProtWarpTomoReconstruct, **kwargs)
        cls.launchProtocol(cls.tomoReconstruct)
        cls.assertIsNotNone(cls.tomoReconstruct.Tomograms, "SetOfTomograms has not been produced.")
        return cls.tomoReconstruct

    @classmethod
    def runImportCtf(cls, **kwargs):
        cls.protImportCtf = cls.newProtocol(ProtImportTsCTF, **kwargs)
        cls.launchProtocol(cls.protImportCtf)
        ctf = getattr(cls.protImportCtf, cls.protImportCtf._possibleOutputs.CTFs.name)
        cls.assertIsNotNone(ctf, "SetOfCTFTomoSeries has not been produced.")
        return cls.protImportCtf

    @classmethod
    def runTomoReconstruct(cls, **kwargs):
        imodProts = Domain.importFromPlugin('imod.protocols', doRaise=True)
        cls.protRecon = cls.newProtocol(
            imodProts.ProtImodTomoReconstruction, **kwargs)
        cls.launchProtocol(cls.protRecon)
        cls.assertIsNotNone(cls.protRecon.Tomograms,
                            "SetOfTomograms has not been produced.")

        return cls.protRecon


class TestDeconvolveTomo(TestWarpBase):

    def test_run(self):
        print(magentaStr("\n==> Importing data - tilt series:"))
        protImportTS = self.runImportTiltSeries(filesPath=self.ts_path,
                                                filesPattern="{TS}.mrcs",
                                                anglesFrom=2,  # tlt file
                                                voltage=300,
                                                sphericalAberration=2.7,
                                                amplitudeContrast=0.07,
                                                magnification=50000,
                                                samplingRate=20.96,
                                                tiltAxisAngle=85.0,
                                                dosePerFrame=1.67)

        print(magentaStr("\n==> Importing data - tomo CTFs:"))
        protImportCtf = self.runImportCtf(filesPath=self.ts_path,
                                          filesPattern="mixedCTEM_tomo*.defocus",
                                          importFrom=1,  # imod
                                          inputSetOfTiltSeries=protImportTS.outputTiltSeries)

        print(magentaStr("\n==> Running imod - tomo reconstruct:"))
        protImportTomo = self.runTomoReconstruct(
            tomoThickness=110.0,
            inputSetOfTiltSeries=protImportTS.outputTiltSeries)

        print(magentaStr("\n==> Testing warp - deconvolve tomograms:"))
        ctf = getattr(protImportCtf, protImportCtf._possibleOutputs.CTFs.name)
        protDeconvTomo = self.newProtocol(
            ProtWarpDeconvTomo,
            inputTomograms=protImportTomo.Tomograms,
            inputCTFs=ctf)

        self.launchProtocol(protDeconvTomo)
        outputTomos = getattr(protDeconvTomo, protDeconvTomo._possibleOutputs.Tomograms.name)
        self.assertIsNotNone(outputTomos, "Warp deconvolve tomograms has failed")
        self.assertSetSize(outputTomos, 2)

        print(magentaStr("\n==> Testing warp - deconvolve tilt-series:"))
        protDeconvTS = self.newProtocol(
            ProtWarpDeconvTS,
            inputTiltSeries=protImportTS.outputTiltSeries,
            inputCTFs=ctf)

        self.launchProtocol(protDeconvTS)
        outputTS = getattr(protDeconvTS, protDeconvTS._possibleOutputs.TiltSeries.name)
        self.assertIsNotNone(outputTS, "Warp deconvolve tilt-series has failed")
        self.assertSetSize(outputTS, 2)


class TestWarpEstimateCTFTomoReconstruction(TestWarpBase):
    excludedViewsDict = {'TS_1': [1, 2, 3, 37, 38, 39]}

    def test_warpCTFEstimationTomoReconstruction(self):
        print(magentaStr("\n==> Importing data - tilt series movies:"))
        protImportTSM = self.runImportTiltSeriesM(filesPath=self.tsm_path,
                                                  filesPattern="*/*.mdoc",
                                                  voltage=300,
                                                  samplingRate=0.79,
                                                  tiltAxisAngle=-85.6,
                                                  dosePerFrame=2.64,
                                                  gainFile=os.path.join(self.tsm_path, 'gain_ref.mrc'))

        print(magentaStr("\n==> Running Warp - align tiltseries movies and ctf estimation "))
        protAlignAndCtf = self.runMotioncorrTSMovieAligment(inputTSMovies=protImportTSM.outputTiltSeriesM,
                                                            binFactor=1,
                                                            x=1, y=1, z=3,
                                                            c_x=2, c_y=2, c_z=1,
                                                            range_max=7,
                                                            defocus_max=8,
                                                            gainFlip=2)

        print(magentaStr("\n==> Running Imod - import transformation matrix "))
        protImportTM = self.runImodImportTMatrix(inputSetOfTiltSeries=protAlignAndCtf.TiltSeries,
                                                 filesPath=os.path.join(self.tsm_path, 'tiltstack/TS_1'),
                                                 filesPattern='*.xf',
                                                 binningTM=13)

        print(magentaStr("\n==> Running Warp - Tomo Reconstruction "))
        ctfEstimationTomoReconstruct = self.runWarpCTFEstimationTomoReconstruction(inputSet=protImportTM.TiltSeries,
                                                                                   binFactor=13,
                                                                                   tomo_thickness=1000,
                                                                                   x_dimension=4400,
                                                                                   y_dimension=6000)
        setOfTomogram = ctfEstimationTomoReconstruct.Tomograms
        self.assertSetSize(setOfTomogram, 1)
        self.assertTrue(setOfTomogram.getSamplingRate() == 10.00)
        self.assertTrue(setOfTomogram.getDim() == (348, 474, 80))

        print(magentaStr("\n==> Running Warp - align tiltseries movies and ctf estimation "))
        protAlignAndCtf = self.runMotioncorrTSMovieAligment(inputTSMovies=protImportTSM.outputTiltSeriesM,
                                                            binFactor=1,
                                                            x=1, y=1, z=3,
                                                            c_x=2, c_y=2, c_z=1,
                                                            range_max=7,
                                                            defocus_max=8,
                                                            gainFlip=2)

        self._excludeTsSetViews(protAlignAndCtf.TiltSeries)

        print(magentaStr("\n==> Running Imod - import transformation matrix (Excluding views)"))
        protImportTM = self.runImodImportTMatrix(inputSetOfTiltSeries=protAlignAndCtf.TiltSeries,
                                                 filesPath=os.path.join(self.tsm_path, 'tiltstack/TS_1'),
                                                 filesPattern='*.xf',
                                                 binningTM=13)

        print(magentaStr("\n==> Running Warp - Tomo Reconstruction (Excluding views) "))
        ctfEstimationTomoReconstruct = self.runWarpCTFEstimationTomoReconstruction(inputSet=protImportTM.TiltSeries,
                                                                                   binFactor=13,
                                                                                   tomo_thickness=1000,
                                                                                   x_dimension=4400,
                                                                                   y_dimension=6000)

    @classmethod
    def _excludeTsSetViews(cls, tsSet):
        tsList = [ts.clone(ignoreAttrs=[]) for ts in tsSet]
        for ts in tsList:
            cls._excludeTsViews(tsSet, ts, cls.excludedViewsDict[ts.getTsId()])

    @staticmethod
    def _excludeTsViews(tsSet, ts, excludedViewsList):
        tiList = [ti.clone() for ti in ts]
        for i, ti in enumerate(tiList):
            if i in excludedViewsList:
                ti._objEnabled = False
                ts.update(ti)
        ts.write()
        tsSet.update(ts)
        tsSet.write()
