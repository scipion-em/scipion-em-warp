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
import traceback
from enum import Enum
from os.path import join
from emtable import Table
from pyworkflow import BETA
import pyworkflow.protocol.params as params
from pyworkflow.utils import Message, redStr, cyanStr, makePath, cleanPath
from reliontomo.convert import readSetOfPseudoSubtomograms
from reliontomo.objects import createSetOfRelionPSubtomograms, RelionSetOfPseudoSubtomograms
from tomo.constants import BOTTOM_LEFT_CORNER
from tomo.objects import SetOfTiltSeries
from tomo.utils import getTsIdsIntersection, getTsIdsDicts
from warp import WARP_TOOLS, TS_EXPORT_PARTICLES, \
    OPTIMISATION_SET_STAR, MATCHING_PARTICLES_STAR, MATCHING_TOMOGRAMS_STAR, TS_DEFOCUS_HAND, \
    TS_IMPORT_ALIGNMENTS, TILTIMAGES_FOLDER, tomoStarFields, Plugin
from warp.protocols.protocol_base_tomo import ProtWarpTomoBase
from warp.utils import writeCtfXml, getTransformInfoFromCoordOrSubtomo, modifyStarFileMultiTable

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
        self.ctfDict = None
        self.tomoDict = None

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
        settingsId = self._insertFunctionStep(self.createTsSettingsStep,
                                              prerequisites=[],
                                              needsGPU=False)

        # Per tilt-series: tomostar -> IMOD alignment files -> CTF XML.
        ctfSteps = []
        for tsId in self.seriesDict:
            sId = self._insertFunctionStep(self.createTsMStarStep, tsId,
                                           prerequisites=settingsId,
                                           needsGPU=False)

            sId = self._insertFunctionStep(self.createTsAliFilesStep, tsId,
                                           prerequisites=sId,
                                           needsGPU=False)

            sId = self._insertFunctionStep(self.createCtfXmlStep, tsId,
                                           prerequisites=sId,
                                           needsGPU=False)
            ctfSteps.append(sId)

        # Particle .star files only depend on the input coordinates / tomoDict
        # (prepared in _initialize), not on any prior step.
        partId = self._insertFunctionStep(self.createParticlesStarStep,
                                          prerequisites=[],
                                          needsGPU=False)

        # Import alignments needs every per-series .tomostar + CTF XML in place.
        prereq = self._insertFunctionStep(self.tsImportAligmentsStep,
                                          prerequisites=ctfSteps,
                                          needsGPU=False)
        if self.invertDefocusHand.get():
            prereq = self._insertFunctionStep(self.tsDefocusHandStep,
                                              prerequisites=prereq,
                                              needsGPU=False)

        expId = self._insertFunctionStep(self.exportParticlesStep,
                                         prerequisites=[partId, prereq],
                                         needsGPU=True)

        outId = self._insertFunctionStep(self.createOutputStep,
                                         prerequisites=expId,
                                         needsGPU=False)

        self._insertFunctionStep(self.cleanIntermediateResults,
                                 prerequisites=outId,
                                 needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        super()._initialize()
        # Get the intersection of tsIds
        tsSet = self._getInSet()
        ctfSet = self.inputSetOfCtfTomoSeries.get()
        coordSet = self.coordinates.get()
        tomoSet = coordSet.getPrecedents()
        commonTsIds = getTsIdsIntersection(coordSet, tsSet, ctfSet)
        self.seriesDict, self.ctfDict, self.tomoDict = getTsIdsDicts(tsSet, ctfSet, tomoSet,
                                                                     present_ts_ids=commonTsIds)
        makePath(*[self._getTiltStackDir(), self._getParticlesFolder()])

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

    def createCtfXmlStep(self, tsId: str):
        """Build the Warp per-series CTF XML directly from Scipion's
        previously estimated CTF metadata. ts_import_alignments later
        re-serialises this file to the full schema."""
        if tsId in self.failedTsIds:
            return
        try:
            logger.info(cyanStr(f">>> {tsId} - Writing CTF XML from Scipion CTF..."))
            ts = self.seriesDict[tsId]
            # Use the intersected, pre-cloned CTF dict (tsId -> CTFTomoSeries clone)
            # rather than querying the whole input set, which may have been subsetted.
            ctfTomoSeries = self.ctfDict.get(tsId)
            if ctfTomoSeries is None:
                raise ValueError(f'No input CTF estimation found for tilt-series {tsId}')
            acq = ts.getAcquisition()
            x, y, _ = ts.getDim()
            xmlPath = join(self._getTiltSeriesDir(), f'{tsId}.xml')
            writeCtfXml(xmlPath, ctfTomoSeries, ts,
                        apix=ts.getSamplingRate(),
                        voltage=acq.getVoltage(),
                        cs=acq.getSphericalAberration(),
                        amplitude=acq.getAmplitudeContrast(),
                        imageDimPix=(x, y))
        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to write the CTF XML '
                                f'with exception {e}. Skipping... '))
            logger.error(traceback.format_exc())
            self.failedTsIds.append(tsId)

    def createParticlesStarStep(self):
        logger.info(cyanStr(">>> Generating the particles .star file per tomogram..."))
        coordSet = self.coordinates.get()
        outPath = self._getParticlesFolder()
        
        # Iterate only the intersected, pre-cloned tomograms (the input sets may have
        # been subsetted); iterCoordinates(volume=tomo) then yields just that
        # tomogram's coordinates, so no .star references a missing .tomostar / CTF XML.
        for tomoId, tomo in self.tomoDict.items():
            if tomoId in self.failedTsIds:
                continue

            # Dimensions and sampling rate from the tomogram element: the set may be
            # heterogeneous, so these are read per tomogram, not from the set.
            tomoDim = tomo.getDim()
            sRate = tomo.getSamplingRate()
            particlesTable = Table(columns=tomoStarFields)

            for coord in coordSet.iterCoordinates(volume=tomo):
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
                    f'{tomoId}.tomostar',
                    coord.getScore()
                )

            # Only write a .star for tomograms that actually have coordinates.
            if len(particlesTable):
                with open(join(outPath, f'{tomoId}.star'), 'w') as f:
                    particlesTable.writeStar(f, tableName='')

    def tsImportAligmentsStep(self):
        logger.info(cyanStr(">>> Starting import alignments..."))
        # ts_import_alignments runs once over every series and takes a single global
        # --alignment_angpix, so it is read from the set.
        args = [
            f"--settings {self._getTiltSeriesSettingsFn()}",
            f"--alignments {self._getTiltStackDir()}",
            f"--alignment_angpix {self._getInSet().getSamplingRate()}",
        ]
        args = ' '.join(args)
        self.runJob(Plugin.getProgram(WARP_TOOLS, TS_IMPORT_ALIGNMENTS), args, executable='/bin/bash')

    def tsDefocusHandStep(self):
        """Flip the per-tilt-series defocus/geometric handedness (Warp's
        AreAnglesInverted flag), which ts_export_particles writes out as the RELION
        tomogram hand (rlnTomoHand) for M. '--set_flip' is a pure metadata setter and
        does not require the per-movie 2x2x1 defocus grid that '--check'/'--set_auto'
        need."""
        logger.info(cyanStr(">>> Flipping defocus handedness (ts_defocus_hand --set_flip)..."))
        args = f"--settings {self._getTiltSeriesSettingsFn()} --set_flip"
        self.runJob(Plugin.getProgram(WARP_TOOLS, TS_DEFOCUS_HAND), args, executable='/bin/bash')

    def exportParticlesStep(self):
        logger.info(cyanStr(">>> Exporting particles..."))
        settingsFile = self._getTiltSeriesSettingsFn()
        particlesDir = self._getParticlesFolder()
        output = self._getRelionFolder()
        # 1. Initialize the list directly with the base arguments
        cmd_parts = [
            f"--settings {settingsFile}",
            f"--input_directory {particlesDir}",
            "--input_pattern *.star",
            f"--output_star {join(output, 'matching.star')}",
            f"--output_angpix {self.output_angpix.get()}",
            f"--box {self.box.get()}",
            f"--diameter {self.diameter.get()}",
            # Split base flags into individual list items for cleaner structure
            "--relative_output_paths",
            "--normalized_coords",
            "--2d",
        ]

        # 2. Conditionally append tilt-sampling controls
        nTilts = self.nTilts.get()
        if nTilts is not None:
            cmd_parts.append(f"--n_tilts {nTilts}")

        cmd_parts.append(f"--max_missing_tilts {self.maxMissingTilts.get()}")

        # 3. Conditionally append normalization flags
        if self.dontNormalizeInput.get():
            cmd_parts.append("--dont_normalize_input")

        if self.dontPremultiply.get():
            cmd_parts.append("--dont_premultiply")

        # 4. Perform an O(N) join
        args =  " ".join(cmd_parts)
        self.runJob(Plugin.getProgram(WARP_TOOLS, TS_EXPORT_PARTICLES), args, executable='/bin/bash')

    def createOutputStep(self):
        coordSet = self.coordinates.get()
        # The output is a single RelionSetOfPseudoSubtomograms carrying one sampling
        # rate and one acquisition, so these global values are read from the set.
        tsSet = self._getInSet()
        tsSRate = tsSet.getSamplingRate()
        boxSize = self.box.get()
        acq = tsSet.getAcquisition()
        relionFolder = self._getRelionFolder()
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

        modifyStarFileMultiTable(join(relionFolder, MATCHING_PARTICLES_STAR),
                                 '_rlnImageName', lambda v: self.normalizeParticlesPath(v))
        modifyStarFileMultiTable(join(relionFolder, MATCHING_TOMOGRAMS_STAR),
                                 '_rlnTomoTiltSeriesName', lambda v: self.normalizeTomogramsPath(v))

        readSetOfPseudoSubtomograms(psubtomoSet)
        outDict = {outputObjects.relionParticles.name: psubtomoSet}
        self._defineOutputs(**outDict)
        self._defineSourceRelation(self.coordinates, psubtomoSet)
        self._defineSourceRelation(self.inputSetOfCtfTomoSeries, psubtomoSet)
        self._defineSourceRelation(self.inputSet, psubtomoSet)

    # --------------------------- UTILS functions -----------------------------
    def _getInSet(self) -> SetOfTiltSeries:
        return self.inputSet.get()

    def createImodFiles(self, ts):
        tsId = ts.getTsId()
        tiltstackFolder = makePath(self._getTsTiltStackDir(tsId))
        ts.writeImodFiles(tiltstackFolder, delimiter=' ', factor=1)
        tltPath = join(tiltstackFolder, f'{ts.getTsId()}.tlt')
        # The .tomostar is generated with NEGATED tilt angles (convert method writeTsStar use -tiltAngle).
        # ts_import_alignments must receive the alignment
        # in the SAME sign convention, so the native IMOD .tlt angles are inverted here.
        # If the tomostar negation ever changes, this inversion must change in lockstep.
        self._invertTiltAngles(tltPath)

    def writeOptimisationSetStar(self, relionFolder):
        optFile = join(relionFolder, OPTIMISATION_SET_STAR)

        particlesFile = self.normalizeOptimizationPath(
            join(relionFolder, MATCHING_PARTICLES_STAR)
        )
        tomogramsFile = self.normalizeOptimizationPath(
            join(relionFolder, MATCHING_TOMOGRAMS_STAR)
        )

        logger.info(cyanStr(">>> Rewriting optimisation set STAR file..."))

        with open(optFile, 'w') as f:
            f.write('data_\n\n')
            f.write(f'_rlnTomoParticlesFile {particlesFile}\n')
            f.write(f'_rlnTomoTomogramsFile {tomogramsFile}\n')

        return optFile

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
        logger.info(cyanStr(">>> Cleaning intermediate results..."))
        imagesFolder = self._getExtraPath(TILTIMAGES_FOLDER)
        cleanPath(imagesFolder)
