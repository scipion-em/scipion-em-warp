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
from pwem.protocols import ProtImportMovies

from tomo.protocols import ProtImportTs, ProtImportTsCTF, ProtImportTsMovies
from tomo.tests import DataSet
from warp.protocols import (ProtWarpTSCtfEstimationTomoReconstruct, ProtWarpTSDefocusHand, ProtWarpDeconvTomo,
                            ProtWarpDeconvTS)


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
        motioncorrProts = Domain.importFromPlugin('motioncorr.protocols', doRaise=True)
        cls.motioncorrTSMovieAligment = cls.newProtocol(motioncorrProts.ProtTsMotionCorr, **kwargs)
        cls.launchProtocol(cls.motioncorrTSMovieAligment)
        cls.assertIsNotNone(cls.motioncorrTSMovieAligment.outputTiltSeries,
                            "SetOfTiltSeries has not been produced.")
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
        cls.ctfEstimationTomoReconstruct = cls.newProtocol(ProtWarpTSCtfEstimationTomoReconstruct, **kwargs)
        cls.launchProtocol(cls.ctfEstimationTomoReconstruct)
        cls.assertIsNotNone(cls.ctfEstimationTomoReconstruct.CTFTomoSeries, "CTFTomoSeries has not been produced.")
        cls.assertIsNotNone(cls.ctfEstimationTomoReconstruct.Tomograms, "SetOfTomograms has not been produced.")
        return cls.ctfEstimationTomoReconstruct

    @classmethod
    def runWarpDefocusHand(cls, **kwargs):
        cls.defocusHand = cls.newProtocol(ProtWarpTSDefocusHand, **kwargs)
        cls.launchProtocol(cls.defocusHand)
        cls.assertIsNotNone(cls.defocusHand, "Warp defocus handedness for all tilt series has failed")
        return cls.defocusHand


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
    def test_warpCTFEstimationTomoReconstruction(self):
        print(magentaStr("\n==> Importing data - tilt series movies:"))
        protImportTSM = self.runImportTiltSeriesM(filesPath=self.tsm_path,
                                                  filesPattern="*/*.mdoc",
                                                  samplingRate=0.79,
                                                  gainFile=os.path.join(self.tsm_path, 'gain_ref.mrc'))

        print(magentaStr("\n==> Running Motioncorr - tiltseries movies "))
        protImportCtf = self.runMotioncorrTSMovieAligment(inputTiltSeriesM=protImportTSM.outputTiltSeriesM,
                                                          binFactor=4.0)

        print(magentaStr("\n==> Running Imod - import transformation matrix "))
        protImportTM = self.runImodImportTMatrix(inputSetOfTiltSeries=protImportCtf.outputTiltSeries,
                                                  filesPath=os.path.join(self.tsm_path, 'tiltstack/TS_1'),
                                                  filesPattern='*.xf')

        print(magentaStr("\n==> Running Warp - Defocus Hand "))
        self.runWarpDefocusHand(inputSet=protImportTM.TiltSeries)

        print(magentaStr("\n==> Running Warp - CTF estimation and Tomo Reconstruction "))
        ctfEstimationTomoReconstruct = self.runWarpCTFEstimationTomoReconstruction(inputSet=protImportTM.TiltSeries,
                                                                                   reconstruct=True,
                                                                                   binFactor=13,
                                                                                   range_high=6.42,
                                                                                   tomo_thickness=1000)
        self.assertSetSize(ctfEstimationTomoReconstruct.CTFTomoSeries, 1)
        setOfTomogram = ctfEstimationTomoReconstruct.Tomograms
        self.assertSetSize(setOfTomogram, 1)
        self.assertTrue(setOfTomogram.getSamplingRate() == 10.00)
        self.assertTrue(setOfTomogram.getDim() == (1172, 1172, 442))

