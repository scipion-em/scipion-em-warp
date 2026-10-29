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
import logging
import os
import traceback
from enum import Enum
from os.path import join

from emtable import Table
from pyworkflow import BETA
import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from pyworkflow.object import Integer
from pyworkflow.utils import Message, redStr, makePath
from reliontomo.convert import readSetOfPseudoSubtomograms
from reliontomo.objects import createSetOfRelionPSubtomograms, RelionSetOfPseudoSubtomograms
from tomo.constants import BOTTOM_LEFT_CORNER
from tomo.objects import SetOfTiltSeries
from warp import MATCHING_FOLDER, TILTSERIES_SETTINGS, RELION_FOLDER, WARP_TOOLS, TS_EXPORT_PARTICLES, \
    OPTIMISATION_SET_STAR, MATCHING_PARTICLES_STAR, MATCHING_TOMOGRAMS_STAR, TS_CTF, TILTSERIES_FOLDER, TS_DEFOCUS_HAND, \
    TS_IMPORT_ALIGNMENTS, TILTIMAGES_FOLDER, tomoStarFields, TILTSTACK_FOLDER
from warp.protocols.protocol_base_tomo import ProtWarpTomoBase
from warp.utils import updateCtFXMLFile, getTransformInfoFromCoordOrSubtomo, modifyStarFileMultiTable

logger = logging.getLogger(__name__)


class outputObjects(Enum):
    relionParticles = RelionSetOfPseudoSubtomograms


class ProtWarpExportParticles(ProtWarpTomoBase):
    """
    Export particles as 3D volumes or 2D image series
    More info:
       https://warpem.github.io/warp/user_guide/warptools/quick_start_warptools_tilt_series/#export-particles
    """

    _label = 'Export particles'
    _devStatus = BETA

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('coordinates', params.PointerParam,
                      pointerClass='SetOfCoordinates3D',
                      label="3D Coordinates",
                      important=True,
                      allowsNull=False,
                      help='3D coordinates')
        form.addParam('inputSet', params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      label="Tilt series",
                      important=True,
                      help='Tilt series with alignment (non interpolated) used in the tomograms reconstruction.')
        form.addParam('inputSetOfCtfTomoSeries',
                      params.PointerParam,
                      label="Input CTF estimation",
                      pointerClass='SetOfCTFTomoSeries',
                      help='Select the CTF estimation for the set '
                           'of tilt-series.')

        form.addSection(label='Reconstruct')
        form.addParam('output_angpix', params.FloatParam,
                      label='Output pixel size',
                      important=True,
                      default=None,
                      help='Pixel size at which to export particles')

        form.addParam('box', params.IntParam,
                      label='Output box size',
                      important=True,
                      default=None,
                      help='Output has this many pixels/voxels on each side')

        form.addParam('diameter', params.IntParam,
                      label='Particles diameter (Å)',
                      important=True,
                      default=None,
                      help='Particle diameter in angstroms')

        form.addParam('nTilts', params.IntParam,
                      label='Number of tilts per particle',
                      default=None, allowsNull=True,
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Number of tilt images to include per particle in the output; images '
                           'with the lowest overall exposure are included first. Leave empty to '
                           'use all tilts (WarpTools --n_tilts).')
        form.addParam('maxMissingTilts', params.IntParam,
                      label='Max. missing tilts',
                      default=5,
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Particles not visible in more than this number of tilts are excluded '
                           '(WarpTools --max_missing_tilts; only applies to 2D particle series, '
                           'which is what this protocol exports).')
        form.addParam('dontNormalizeInput', params.BooleanParam,
                      label="Don't normalize input images",
                      default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Skip normalization of the input tilt images (WarpTools '
                           '--dont_normalize_input).')
        form.addParam('dontPremultiply', params.BooleanParam,
                      label="Don't premultiply by CTF / weights",
                      default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      help="Don't premultiply the 2D particle series by CTFs or RELION weights "
                           '(WarpTools --dont_premultiply). Set this if your downstream STA tool '
                           'expects un-weighted particles.')

        form.addSection(label='Geometry / handedness')
        form.addParam('flipZHand', params.BooleanParam, default=False,
                      label='Invert Z coordinate (handedness)',
                      help="Mirror particle Z coordinates through the tomogram mid-plane "
                           "(z -> Nz-1-z) to reconcile the geometric handedness between the "
                           "picking tomogram and Warp/M.\n\n"
                           "NOTE: this is an improper (chirality-inverting) transform. It is "
                           "correct for picking-stage coordinates whose orientations are the "
                           "identity (the common case). If the input coordinates carry refined "
                           "orientations, those Euler angles must be mirror-transformed "
                           "consistently, which this option does NOT do.")
        form.addParam('invertDefocusHand', params.BooleanParam, default=False,
                      label='Invert defocus handedness',
                      help="Run 'WarpTools ts_defocus_hand --set_flip' before export so the "
                           "exported RELION tomogram hand (rlnTomoHand, derived from Warp's "
                           "AreAnglesInverted flag) is flipped. Use this if M/RELION refinement "
                           "converges to a mirror-image map. Independent of the Z-coordinate "
                           "option above. Unchecked leaves Warp's default handedness untouched.")

        # form.addParam('writeStacks', params.EnumParam,
        #               label='Export type',
        #               default=0,
        #               choices=['2D', '3D'],
        #               display=params.EnumParam.DISPLAY_HLIST,
        #               help='If set to 2D, this program will write output particles as 2d image series centered on '
        #                    'the particle (particle series). If set to 3D, this program will write output '
        #                    'particles as 3d images (subtomograms)')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._initialize()
        pId = self._insertFunctionStep(self.createTsSettingsStep,
                                       prerequisites=[],
                                       needsGPU=False)

        for tsId, tsM in self.seriesDict.items():
            pId = self._insertFunctionStep(self.createTsMStarStep, tsId,
                                           prerequisites=pId,
                                           needsGPU=False)
            pId = self._insertFunctionStep(self.createTsAliFilesStep, tsId,
                                           prerequisites=pId,
                                           needsGPU=False)



        # self._insertFunctionStep(self.prepareDataStep, needsGPU=False)
        self._insertFunctionStep(self.exportParticlesStep, needsGPU=True)
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)
        self._insertFunctionStep(self.cleanIntermediateResults, needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        super()._initialize()
        makePath(self._getTiltStackDir())

    def createTsAliFilesStep(self, tsId: str):
        """Take advantage of Scipion's method to generate IMOD files (compatible with Warp)."""
        if tsId in self.failedTsIds:
            return
        try:
            ts = self.seriesDict[tsId]
            self.createImodFiles(ts)
        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to create the TS star file '
                                f'with exception {e}. Skipping... '))
            logger.error(traceback.format_exc())
            self.failedTsIds.append(tsId)

    def prepareDataStep(self):
        inputTs = self.inputSet.get()
        coordSet = self.coordinates.get()
        tsSr = inputTs.getSamplingRate()
        tomoSr = coordSet.getSamplingRate()
        tomoDim = coordSet.getPrecedents().getDim()
        scaleFactor = tomoSr / tsSr
        self.tomo_thickness = Integer(round(tomoDim[2] * scaleFactor))
        self.x_dimension = Integer(round(tomoDim[0] * scaleFactor))
        self.y_dimension = Integer(round(tomoDim[1] * scaleFactor))

        # Create the warp settings file and retrieve ctf and alignments values
        self.info(">>> Creating the warp settings file and retrieve ctf and alignments values...")
        self.createTiltSeriesSetting(None)
        for ts in inputTs.iterItems(iterate=False):
            self.tsDataPrepare(ts)
            self.createImodFiles(ts)
        self.tsCtfEstimation(tsSr)
        self.updateCTFValues()
        self.tsImportAligments()
        if self.invertDefocusHand.get():
            self.tsDefocusHand()

        sRate = coordSet.getSamplingRate()
        outPath = self._getExtraPath(MATCHING_FOLDER)
        pwutils.makePath(outPath)

        currentTomoId = None
        particlesTable = None
        f = None

        # Generate .star file per tomogram
        self.info(">>> Generate .star file per tomogram...")
        for coord in coordSet.iterCoordinates(orderBy='_tomoId'):
            tomoId = coord.getTomoId()

            if tomoId != currentTomoId:
                if f is not None:
                    particlesTable.writeStar(f, tableName='')
                    f.close()

                currentTomoId = tomoId
                tomoDim = coord.getVolume().getDim()
                filePath = os.path.join(outPath, f'{currentTomoId}.star')
                f = open(filePath, 'w')
                particlesTable = Table(columns=tomoStarFields)

            # Orientations: use the coordinate's transform when present; otherwise (pure picks
            # with no orientation, e.g. a de-novo STA start) write zeroed Euler angles so Warp
            # gets valid rlnAngleRot/Tilt/Psi instead of choking on a missing transform.
            if coord.hasTransform():
                angles, _ = getTransformInfoFromCoordOrSubtomo(coord, sRate)
            else:
                angles = (0., 0., 0.)
            zCoord = coord.getZ(BOTTOM_LEFT_CORNER)
            if self.flipZHand.get():
                # Mirror Z through the tomogram mid-plane (corner origin, 0-indexed
                # voxels) to reconcile Warp/M handedness with the picking tomogram.
                # Improper transform: only valid for identity-orientation coordinates
                # (see the 'Invert Z coordinate' help). Normalisation by tomoDim[2]
                # below is kept identical to X/Y to preserve the existing convention.
                zCoord = (tomoDim[2] - 1) - zCoord
            particlesTable.addRow(
                coord.getX(BOTTOM_LEFT_CORNER) / tomoDim[0],
                coord.getY(BOTTOM_LEFT_CORNER) / tomoDim[1],
                zCoord / tomoDim[2],
                angles[0],
                angles[1],
                angles[2],
                f'{currentTomoId}.tomostar',
                coord.getScore()
            )
        # Close the last open file
        if f is not None:
            particlesTable.writeStar(f, tableName='')
            f.close()

    def exportParticlesStep(self):
        self.info(">>> Exporting particles...")
        settingFile = self._getExtraPath(TILTSERIES_SETTINGS)
        matchinFolder = self._getExtraPath(MATCHING_FOLDER)
        output = self._getExtraPath(RELION_FOLDER)
        pwutils.makePath(output)
        argsDict = {
            "--settings": os.path.abspath(settingFile),
            "--input_directory": matchinFolder,
            "--input_pattern": "*.star",
            "--output_star": os.path.join(output, 'matching.star'),
            "--output_angpix": self.output_angpix.get(),
            "--box": self.box.get(),
            "--diameter": self.diameter.get(),
        }
        cmd = '--relative_output_paths --normalized_coords --2d'
        # Tilt-sampling controls
        if self.nTilts.get() is not None:
            cmd += f' --n_tilts {self.nTilts.get()}'
        cmd += f' --max_missing_tilts {self.maxMissingTilts.get()}'
        # Normalization / premultiply flags
        if self.dontNormalizeInput.get():
            cmd += ' --dont_normalize_input'
        if self.dontPremultiply.get():
            cmd += ' --dont_premultiply'
        self.runProgram(argsDict, WARP_TOOLS, TS_EXPORT_PARTICLES, othersCmds=cmd)

    def writeOptimisationSetStar(self, relionFolder):
        optFile = os.path.join(relionFolder, OPTIMISATION_SET_STAR)

        particlesFile = self.normalizeOptimizationPath(
            os.path.join(relionFolder, MATCHING_PARTICLES_STAR)
        )
        tomogramsFile = self.normalizeOptimizationPath(
            os.path.join(relionFolder, MATCHING_TOMOGRAMS_STAR)
        )

        self.info(">>> Rewriting optimisation set STAR file...")

        with open(optFile, 'w') as f:
            f.write('data_\n\n')
            f.write(f'_rlnTomoParticlesFile {particlesFile}\n')
            f.write(f'_rlnTomoTomogramsFile {tomogramsFile}\n')

        return optFile

    def createOutputStep(self):
        coordSet = self.coordinates.get()
        tsSet = self.inputSet.get()
        tsSRate = tsSet.getSamplingRate()
        boxSize = self.box.get()
        acq = tsSet.getAcquisition()
        relionFolder = self._getExtraPath(RELION_FOLDER)
        are2dStacks = True

        optFile = self.writeOptimisationSetStar(relionFolder)
        psubtomoSet = createSetOfRelionPSubtomograms(self._getPath(),
                                                     optFile,
                                                     coordSet,
                                                     template='pseudosubtomograms%s.sqlite',
                                                     tsSamplingRate=tsSRate,
                                                     relionBinning=self.output_angpix.get() / tsSet.getSamplingRate(),
                                                     boxSize=boxSize,
                                                     are2dStacks=are2dStacks,
                                                     acquisition=acq)

        modifyStarFileMultiTable(os.path.join(relionFolder, MATCHING_PARTICLES_STAR),
                                 '_rlnImageName', lambda v: self.normalizeParticlesPath(v))
        modifyStarFileMultiTable(os.path.join(relionFolder, MATCHING_TOMOGRAMS_STAR),
                                 '_rlnTomoTiltSeriesName', lambda v: self.normalizeTomogramsPath(v))

        readSetOfPseudoSubtomograms(psubtomoSet)
        outDict = {outputObjects.relionParticles.name: psubtomoSet}
        self._defineOutputs(**outDict)
        self._defineSourceRelation(self.coordinates, psubtomoSet)
        self._defineSourceRelation(self.inputSetOfCtfTomoSeries, psubtomoSet)
        self._defineSourceRelation(self.inputSet, psubtomoSet)

    def tsCtfEstimation(self, tsSr):
        """CTF estimation"""

        self.info(">>> Generating ctf estimation...")
        settingFile = self._getExtraPath(TILTSERIES_SETTINGS)
        argsDict = {
            "--settings": os.path.abspath(settingFile),
            "--range_high": tsSr * 3,
            "--range_low": tsSr * 4,
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
            updateCtFXMLFile(defocusFilePath, ctfTomoSeries, ts)

    def tsDefocusHand(self):
        """Flip the per-tilt-series defocus/geometric handedness (Warp's
        AreAnglesInverted flag), which ts_export_particles writes out as the RELION
        tomogram hand (rlnTomoHand) for M. '--set_flip' is a pure metadata setter and
        does not require the per-movie 2x2x1 defocus grid that '--check'/'--set_auto'
        need."""
        self.info(">>> Flipping defocus handedness (ts_defocus_hand --set_flip)...")
        settingFile = self._getExtraPath(TILTSERIES_SETTINGS)
        argsDict = {"--settings": os.path.abspath(settingFile)}
        self.runProgram(argsDict, WARP_TOOLS, TS_DEFOCUS_HAND, othersCmds='--set_flip')

    def tsImportAligments(self):
        self.info(">>> Starting import alignments...")
        processingFolder = os.path.abspath(self._getExtraPath(TILTSERIES_FOLDER))
        tiltstackFolder = os.path.join(processingFolder, 'tiltstack')
        angpix = self.inputSet.get().getSamplingRate()
        settingFile = self._getExtraPath(TILTSERIES_SETTINGS)
        argsDict = {
            "--settings": os.path.abspath(settingFile),
            '--alignments': os.path.abspath(tiltstackFolder),
            "--alignment_angpix": angpix,
        }
        cmd = ' '.join(['%s %s' % (k, v) for k, v in argsDict.items()])
        self.runJob(self.getPlugin().getProgram(WARP_TOOLS, TS_IMPORT_ALIGNMENTS), cmd, executable='/bin/bash')

    @staticmethod
    def _invertTiltAngles(tltPath):
        with open(tltPath, 'r') as f:
            angles = [
                float(line.strip())
                for line in f
                if line.strip()
            ]

        with open(tltPath, 'w') as f:
            for angle in angles:
                f.write(f'{-angle}\n')

    def cleanIntermediateResults(self):
        self.info(">>> Cleaning intermediate results...")
        imagesFolder = self._getExtraPath(TILTIMAGES_FOLDER)
        pwutils.cleanPath(imagesFolder)

    # --------------------------- UTILS functions -----------------------------
    def _getInSet(self) -> SetOfTiltSeries:
        return self.inputSet.get()

    def _getTiltStackDir(self) -> str:
        return join(self._getTiltSeriesDir(), TILTSTACK_FOLDER)

    def _getTsTiltStackDir(self, tsId: str) -> str:
        return join(self._getTiltStackDir(), tsId)

    def createImodFiles(self, ts):
        tsId = ts.getTsId()
        tiltstackFolder = makePath(self._getTsTiltStackDir(tsId))
        ts.writeImodFiles(tiltstackFolder, delimiter=' ', factor=1)
        tltPath = os.path.join(tiltstackFolder, f'{ts.getTsId()}.tlt')
        # The .tomostar is generated with NEGATED tilt angles (convert method writeTsStar use -tiltAngle).
        # ts_import_alignments must receive the alignment
        # in the SAME sign convention, so the native IMOD .tlt angles are inverted here.
        # If the tomostar negation ever changes, this inversion must change in lockstep.
        self._invertTiltAngles(tltPath)




