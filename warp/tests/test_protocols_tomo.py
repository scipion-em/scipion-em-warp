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

from pwem import Domain
from pyworkflow.tests import BaseTest, setupTestProject
from pyworkflow.utils import magentaStr

from tomo.protocols import ProtImportTs, ProtImportTsCTF
from tomo.tests import DataSet

from warp.protocols.protocol_deconv_tomo import ProtWarpDeconvTomo
from warp.protocols.protocol_deconv_ts import ProtWarpDeconvTS


class TestDeconvolveTomo(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.inputDataSet = DataSet.getDataSet('empiar10064')
        cls.ts_path = cls.inputDataSet.getPath()

    @classmethod
    def runImportTiltSeries(cls, **kwargs):
        cls.protImportTS = cls.newProtocol(ProtImportTs, **kwargs)
        cls.launchProtocol(cls.protImportTS)
        cls.assertIsNotNone(cls.protImportTS.outputTiltSeries,
                            "SetOfTiltSeries has not been produced.")
        return cls.protImportTS

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
