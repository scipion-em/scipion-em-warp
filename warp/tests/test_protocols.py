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

from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from pyworkflow.utils import magentaStr
from pwem.protocols import ProtImportMicrographs, ProtImportCTF

from warp.protocols.protocol_deconv_mics import ProtWarpDeconvMics, outputs


class TestDeconvolveMics(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsXmipp = DataSet.getDataSet('xmipp_tutorial')
        cls.micsFn = cls.dsXmipp.getFile('allMics')
        cls.dsGrigorieff = DataSet.getDataSet('grigorieff')
        cls.ctfFn = cls.dsGrigorieff.getFile('ctffind3/ctfs.sqlite')

    @classmethod
    def runImportMics(cls, filesPath, samplingRate):
        """ Import micrographs. """
        print(magentaStr("\n==> Importing data - micrographs:"))
        protImport = cls.newProtocol(ProtImportMicrographs,
                                     filesPath=filesPath,
                                     samplingRate=samplingRate)
        cls.launchProtocol(protImport)
        output = getattr(protImport, protImport._possibleOutputs.outputMicrographs.name)
        cls.assertIsNotNone(output,
                            "SetOfMicrographs has not been produced.")

        return protImport

    @classmethod
    def runImportCTFs(cls, filesPath, inputMics):
        print(magentaStr("\n==> Importing data - ctfs:"))
        protCTF = cls.newProtocol(ProtImportCTF,
                                  importFrom=ProtImportCTF.IMPORT_FROM_SCIPION,
                                  filesPath=filesPath,
                                  inputMicrographs=inputMics)
        cls.launchProtocol(protCTF)
        cls.assertIsNotNone(protCTF.outputCTF,
                            "There was a problem when importing CTFs.")
        return protCTF

    def test_run(self):
        protImport = self.runImportMics(self.micsFn, 1.237)
        mics = getattr(protImport, protImport._possibleOutputs.outputMicrographs.name)
        protCTF = self.runImportCTFs(self.ctfFn, mics)

        print(magentaStr("\n==> Testing warp - deconvolve micrographs:"))
        protDeconv2D = self.newProtocol(
            ProtWarpDeconvMics,
            inputMicrographs=mics,
            ctfRelations=protCTF.outputCTF)
        self.launchProtocol(protDeconv2D)
        outputMics = getattr(protDeconv2D, outputs.Micrographs.name)
        self.assertIsNotNone(outputMics, "Warp deconvolve micrographs has failed")
        self.assertSetSize(outputMics, 3)
