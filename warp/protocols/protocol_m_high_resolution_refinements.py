# ******************************************************************************
# *
# * Authors:     Yunior C. Fonseca Reyna
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
# ******************************************************************************

import os
import time

import mrcfile

from pyworkflow import BETA
import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from pyworkflow.object import Integer
from pyworkflow.protocol import GPU_LIST

from warp.constants import *
from warp.protocols.protocol_base import ProtWarpBase
from warp.utils import updateCtFXMLFile


class ProtWarpMHigResolutionRefinement(ProtWarpBase):
    """
    High Resolution Refinements in M
    More info:
        https://warpem.github.io/warp/user_guide/warptools/quick_start_warptools_tilt_series/#initial-3d-refinement-in-relion
    """

    _label = 'M high resolution refinement'
    _devStatus = BETA

    def __init__(self, **kwargs):
        ProtWarpBase.__init__(self, **kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection("Import")
        form.addParam('inputSet',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      label='Input set of tilt-series',
                      help='Input set of tilt-series')
        form.addParam('inputSetOfCtfTomoSeries',
                      params.PointerParam,
                      label="Input CTF estimation",
                      pointerClass='SetOfCTFTomoSeries',
                      help='Select the CTF estimation for the set '
                           'of tilt-series.')
        form.addParam('relionRefineProt', params.PointerParam,
                      pointerClass='ProtRelionRefineSubtomograms',
                      label='Relion 3D refinement',
                      help='Relion 3D refinement protocol')
        form.addParam('refMask', params.PointerParam, pointerClass='VolumeMask, Volume',
                      default=None,
                      label='Mask',
                      help='Path to a tight binary mask file. M will automatically expand and smooth it based on '
                           'current resolution')

        form.addSection("Species")
        form.addParam('diameter', params.IntParam, default=None,
                      important=True,
                      label='Molecule diameter (Å)',
                      help=' Molecule diameter in Angstrom.')
        form.addParam('symmetry', params.StringParam, default='C1',
                      label='Symmetry',
                      help='Default: C1. Symmetry of the template, e.g. C1, D7, O')
        form.addParam('temporal_samples', params.IntParam, default=1,
                      label='Temporal samples',
                      help="Number of temporal samples in each particle pose's trajectory.")
        form.addParam('angpix_resample', params.FloatParam, default=None,
                      allowsNull=True,
                      label='Angpix resample',
                      help="Resample half-maps and masks to this pixel size.")
        form.addParam('lowpass', params.FloatParam, default=None,
                      label='Gaussian low-pass filter',
                      help="Optional low-pass filter (in Å), applied to both half-maps")

        form.addSection("Refinement")

        form.addHidden(GPU_LIST, params.StringParam, default='0',
                       label='Choose GPU IDs:', validators=[params.NonEmpty],
                       help='Space-separated list of GPU IDs to use for processing. '
                            'Default: all GPUs in the system For example: "0 1 5"')

    def _insertAllSteps(self):
        self._insertFunctionStep(self.prepareDataStep, needsGPU=True)
        self._insertFunctionStep(self.createPopulationStep, needsGPU=False)
        self._insertFunctionStep(self.createSourcesStep, needsGPU=False)
        self._insertFunctionStep(self.createSpeciesStep, needsGPU=True)
        self._insertFunctionStep(self.refinementStep, needsGPU=True)

    def prepareDataStep(self):
        inputTs = self.inputSet.get()
        coordSet = self.relionRefineProt.get().inReParticles.get().getCoordinates3D()
        tsSr = inputTs.getSamplingRate()
        tomoSr = coordSet.getSamplingRate()
        tomoDim = coordSet.getPrecedents().getDim()
        scaleFactor = tomoSr / tsSr
        self.tomo_thickness = Integer(round(tomoDim[2] * scaleFactor))
        self.x_dimension = Integer(round(tomoDim[0] * scaleFactor))
        self.y_dimension = Integer(round(tomoDim[1] * scaleFactor))

        self.createTiltSeriesSetting(None)
        for ts in inputTs.iterItems(iterate=False):
            self.tsDataPrepare(ts)
            self.createImodFiles(ts)
        self.tsCtfEstimation()
        self.updateCTFValues()
        self.tsImportAligments()

    def createPopulationStep(self):
        self.info(">>> Creating population...")
        tsResults = os.path.join(self._getExtraPath(), 'm')
        argsDict = {
            "--directory": tsResults,
            "--name": 'processing'
        }
        self.runProgram(argsDict, MTOOLS, CREATE_POPULATION)

    def createSourcesStep(self):
        self.info(">>> Creating sources...")
        populationPath = os.path.join(self._getExtraPath('m'))
        settingFile = self._getExtraPath(TILTSERIE_SETTINGS)
        argsDict = {
            "--name": 'processing',
            "--population": os.path.join(populationPath, 'processing.population'),
            "--processing_settings": os.path.abspath(settingFile)
        }
        self.runProgram(argsDict, MTOOLS, CREATE_SOURCES)

    def createSpeciesStep(self):
        relionProt = self.relionRefineProt.get()
        extraPath = relionProt._getExtraPath()
        averageSubTomoHalfMaps = relionProt.average.getHalfMaps(asList=True)
        populationPath = os.path.join(self._getExtraPath('m'))
        mask = self.refMask.get().getFileName()
        voxelSize = relionProt.average.getSamplingRate()

        argsDict = {
            "--population": os.path.join(populationPath, 'processing.population'),
            "--name": 'processing_species',
            "--diameter": self.diameter.get(),
            "--sym": self.symmetry.get(),
            "--temporal_samples": self.temporal_samples.get(),
            "--half1": averageSubTomoHalfMaps[0],
            "--half2": averageSubTomoHalfMaps[1],
            "--mask": mask,
            "--particles_relion": os.path.join(extraPath, '_data.star'),
            "--lowpass": self.lowpass.get(),
            "--angpix": voxelSize

        }
        if self.angpix_resample.get() is not None:
            argsDict["--angpix_resample"] = self.angpix_resample.get()

        self.runProgram(argsDict, MTOOLS, CREATE_SPECIES)

    def refinementStep(self):
        populationPath = os.path.join(self._getExtraPath('m'))
        self.info(">>> Running M to check setup...")
        self.checkSetup(populationPath)
        self.info(">>> First M Refinement with 2D Image Warp, Particle Poses Refinement and CTF Refinement")
        self.firstRefinement(populationPath)
        self.info(">>> Stage Angle Refinement")
        self.stageAngleRefinement(populationPath)

    def checkSetup(self, populationPath):
        argsDict = {
            "--population": os.path.join(populationPath, 'processing.population'),
            "--iter": 0
        }
        self.runProgram(argsDict, MCORE, None)

    def firstRefinement(self, populationPath):
        argsDict = {
            "--population": os.path.join(populationPath, 'processing.population'),
            "--refine_imagewarp": '6x4'
        }
        cmd = '--refine_particles --ctf_defocus --ctf_defocusexhaustive'
        self.runProgram(argsDict, MCORE, None, othersCmds=cmd)

    def stageAngleRefinement(self, populationPath):
        argsDict = {
            "--population": os.path.join(populationPath, 'processing.population'),
            "--refine_imagewarp": '6x4'
        }
        cmd = '--refine_particles --refine_stageangles '
        self.runProgram(argsDict, MCORE, None, othersCmds=cmd)

    def tsCtfEstimation(self):
        """CTF estimation"""

        self.info(">>> Generating ctf estimation...")
        settingFile = self._getExtraPath(TILTSERIE_SETTINGS)
        argsDict = {
            "--settings": os.path.abspath(settingFile)
        }
        try:
            self.runProgram(argsDict, WARP_TOOLS, TS_CTF)
        except Exception:
            self.info(">>> Error generating ctf estimation file...")

    def updateCTFValues(self):
        for ts in self.inputSet.get().iterItems():
            ctfTomoSeries = self.inputSetOfCtfTomoSeries.get().getItem('_tsId', ts.getTsId())
            processingFolder = os.path.abspath(self._getExtraPath(TILTSERIES_FOLDER))
            defocusFilePath = os.path.join(processingFolder, ts.getTsId() + '.xml')
            updateCtFXMLFile(defocusFilePath, ctfTomoSeries)

    def tsImportAligments(self):
        self.info(">>> Starting import aligments...")
        processingFolder = os.path.abspath(self._getExtraPath(TILTSERIES_FOLDER))
        tiltstackFolder = os.path.join(processingFolder, 'tiltstack')
        angpix = self.inputSet.get().getSamplingRate()
        settingFile = self._getExtraPath(TILTSERIE_SETTINGS)
        argsDict = {
            "--settings": os.path.abspath(settingFile),
            '--alignments': os.path.abspath(tiltstackFolder),
            "--alignment_angpix": angpix
        }
        cmd = ' '.join(['%s %s' % (k, v) for k, v in argsDict.items()])
        self.runJob(self.getPlugin().getProgram(WARP_TOOLS, TS_IMPORT_ALIGNMENTS), cmd, executable='/bin/bash')

    def createImodFiles(self, ts):
        processingFolder = os.path.abspath(self._getExtraPath(TILTSERIES_FOLDER))
        tiltstackFolder = os.path.join(processingFolder, 'tiltstack', ts.getTsId())
        pwutils.makePath(tiltstackFolder)
        ts.writeImodFiles(tiltstackFolder, delimiter=' ')
