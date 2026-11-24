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
from glob import glob
from os.path import join
from typing import Set, List
from emtable import Table
from pyworkflow import BETA
from pyworkflow.protocol import PointerParam, FloatParam, IntParam, LEVEL_ADVANCED, BooleanParam
from pyworkflow.utils import Message, redStr, cyanStr, yellowStr, makePath, cleanPath, createLink
from pwem.emlib.image.image_readers import ImageStack, ImageReadersRegistry
from reliontomo.convert import readSetOfPseudoSubtomograms
from reliontomo.objects import createSetOfRelionPSubtomograms, RelionSetOfPseudoSubtomograms
from tomo.constants import BOTTOM_LEFT_CORNER
from tomo.objects import SetOfTiltSeries, SetOfCTFTomoSeries, SetOfCoordinates3D, TiltSeries
from tomo.utils import getTsIdsIntersection, getTsIdsDicts, getPresentAcqOrdersTsAndCtf
from warp import WARP_TOOLS, TS_EXPORT_PARTICLES, \
    OPTIMISATION_SET_STAR, MATCHING_PARTICLES_STAR, MATCHING_TOMOGRAMS_STAR, TS_DEFOCUS_HAND, \
    TS_IMPORT_ALIGNMENTS, TOMOSTAR_EXT, tomoStarFields, Plugin
from warp.convert import writeTiltImagesTsStar
from warp.protocols.protocol_base_tomo import ProtWarpTomoBase
from warp.utils import writeCtfXml, getTransformInfoFromCoordOrSubtomo, modifyStarFileMultiTable, \
    getWarpTiltImages

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
        form.addParam('coordinates', PointerParam,
                      pointerClass='SetOfCoordinates3D',
                      label="3D Coordinates",
                      important=True,
                      allowsNull=False)

        form.addParam('inTsSet', PointerParam,
                      pointerClass='SetOfTiltSeries',
                      label="Tilt series",
                      important=True,
                      allowsNull=False,
                      help='Tilt series with alignment (non interpolated) used in the tomograms reconstruction.')

        form.addParam('inCtfSet', PointerParam,
                      pointerClass='SetOfCTFTomoSeries',
                      label="Input CTF estimation",
                      important=True,
                      allowsNull=False,
                      help='CTF estimation for the set of tilt-series.')

        form.addSection(label='Reconstruct')
        form.addParam('output_angpix', FloatParam,
                      label='Output pixel size',
                      important=True,
                      default=None,
                      help='Pixel size at which to export particles')

        form.addParam('box', IntParam,
                      label='Output box size',
                      important=True,
                      default=None,
                      help='Output has this many pixels/voxels on each side')

        form.addParam('diameter', IntParam,
                      label='Particles diameter (Å)',
                      important=True,
                      default=None,
                      help='Particle diameter in angstroms')

        form.addParam('nTilts', IntParam,
                      label='Number of tilts per particle',
                      default=None, allowsNull=True,
                      expertLevel=LEVEL_ADVANCED,
                      help='Number of tilt images to include per particle in the output; images '
                           'with the lowest overall exposure are included first. Leave empty to '
                           'use all tilts (WarpTools --n_tilts).')
        form.addParam('maxMissingTilts', IntParam,
                      label='Max. missing tilts',
                      default=5,
                      expertLevel=LEVEL_ADVANCED,
                      help='Particles not visible in more than this number of tilts are excluded '
                           '(WarpTools --max_missing_tilts; only applies to 2D particle series, '
                           'which is what this protocol exports).')
        form.addParam('dontNormalizeInput', BooleanParam,
                      label="Don't normalize input images",
                      default=False,
                      expertLevel=LEVEL_ADVANCED,
                      help='Skip normalization of the input tilt images (WarpTools '
                           '--dont_normalize_input).')
        form.addParam('dontPremultiply', BooleanParam,
                      label="Don't premultiply by CTF / weights",
                      default=False,
                      expertLevel=LEVEL_ADVANCED,
                      help="Don't premultiply the 2D particle series by CTFs or RELION weights "
                           '(WarpTools --dont_premultiply). Set this if your downstream STA tool '
                           'expects un-weighted particles.')

        form.addSection(label='Geometry / handedness')
        form.addParam('flipZHand', BooleanParam, default=False,
                      label='Invert Z coordinate (handedness)',
                      help="Mirror particle Z coordinates through the tomogram mid-plane "
                           "(z -> Nz-1-z) to reconcile the geometric handedness between the "
                           "picking tomogram and Warp/M.\n\n"
                           "NOTE: this is an improper (chirality-inverting) transform. It is "
                           "correct for picking-stage coordinates whose orientations are the "
                           "identity (the common case). If the input coordinates carry refined "
                           "orientations, those Euler angles must be mirror-transformed "
                           "consistently, which this option does NOT do.")

        form.addParam('invertDefocusHand', BooleanParam, default=False,
                      label='Invert defocus handedness',
                      help="Run 'WarpTools ts_defocus_hand --set_flip' before export so the "
                           "exported RELION tomogram hand (rlnTomoHand, derived from Warp's "
                           "AreAnglesInverted flag) is flipped. Use this if M/RELION refinement "
                           "converges to a mirror-image map. Independent of the Z-coordinate "
                           "option above. Unchecked leaves Warp's default handedness untouched.")

        # form.addParam('writeStacks', EnumParam,
        #               label='Export type',
        #               default=0,
        #               choices=['2D', '3D'],
        #               display=EnumParam.DISPLAY_HLIST,
        #               help='If set to 2D, this program will write output particles as 2d image series centered on '
        #                    'the particle (particle series). If set to 3D, this program will write output '
        #                    'particles as 3d images (subtomograms)')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        commonTsIds = self._init()
        settingsId = self._insertFunctionStep(self.createTsSettingsStep,
                                              prerequisites=[],
                                              needsGPU=False)

        # Per tilt-series: tilt images + tomostar -> IMOD alignment files -> CTF XML.
        ctfSteps = []
        for tsId in commonTsIds:
            ts = self.seriesDict[tsId]
            ctf = self.ctfDict[tsId]
            presentAcqOrders = getPresentAcqOrdersTsAndCtf(ts, ctf)
            if not presentAcqOrders:
                # Warp loads every .tomostar in the folder and aborts the whole command on
                # an empty one, so such a series must not get any step (nor a .tomostar).
                logger.warning(yellowStr(f'{tsId}: no views are enabled in both the tilt-series '
                                         f'and the CTF. Skipping...'))
                # Recorded so createParticlesStarStep doesn't write a .star for it either.
                self.failedTsIds.append(tsId)
                continue
            sId = self._insertFunctionStep(self.prepareTsDataStep,
                                           tsId,
                                           presentAcqOrders,
                                           prerequisites=settingsId,
                                           needsGPU=False)

            sId = self._insertFunctionStep(self.createTsAliFilesStep,
                                           tsId,
                                           presentAcqOrders,
                                           prerequisites=sId,
                                           needsGPU=False)

            sId = self._insertFunctionStep(self.createCtfXmlStep,
                                           tsId,
                                           presentAcqOrders,
                                           prerequisites=sId,
                                           needsGPU=False)
            ctfSteps.append(sId)

        # Particle .star files are written after every per-series step so that the
        # series that failed there (failedTsIds) are already known and skipped.
        partId = self._insertFunctionStep(self.createParticlesStarStep,
                                          prerequisites=ctfSteps,
                                          needsGPU=False)

        # Import alignments needs every per-series .tomostar + CTF XML in place.
        prereq = self._insertFunctionStep(self.tsImportAligmentsStep,
                                          prerequisites=partId,
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
    def _init(self) -> Set[str]:
        super()._initialize()
        # Get the intersection of tsIds
        tsSet = self._getInTsSet()
        ctfSet = self._getInCtfSet()
        coordSet = self._getInCoords()
        tomoSet = coordSet.getPrecedents()
        commonTsIds = getTsIdsIntersection(coordSet, tsSet, ctfSet)
        self.seriesDict, self.ctfDict, self.tomoDict = getTsIdsDicts(tsSet, ctfSet, tomoSet,
                                                                     present_ts_ids=commonTsIds)
        makePath(*[self._getTiltStackDir(), self._getParticlesFolder(), self._getAveragesDir()])
        return commonTsIds

    def prepareTsDataStep(self, tsId: str, presentAcqOrders: List[int]):
        """Give Warp the per-tilt images it reads in place of fs_motion output, and
        the .tomostar that points to them.

        Warp opens each tilt as a "movie" at <tomostar dir>/<wrpMovieName>: it reads the
        header of that file, then the image itself from the sibling average/ folder
        (TiltSeries.LoadMovieData / LoadMovieSizes). So each enabled tilt is written as
        a single-image .mrc in tiltimages/, linked from tiltimages/average/, and the
        .tomostar rows point to ../tiltimages/<name>. No Warp frame-series metadata is
        needed; a missing per-movie XML is ignored by Warp."""
        if tsId in self.failedTsIds:
            return
        try:
            logger.info(cyanStr(f">>> {tsId} - Extracting tilt images and writing the tomostar file..."))
            ts = self.seriesDict[tsId]
            tiltImages = getWarpTiltImages(ts, presentAcqOrders)
            imagesDir = self._getTiltImagesDir()
            averagesDir = self._getAveragesDir()
            # The pixel size goes to the MRC header: Warp takes the physical image size from it.
            properties = {'sr': ts.getSamplingRate()}
            imageNames = []
            for ti in tiltImages:
                imageName = f'{tsId}_{ti.getIndex():03d}.mrc'
                tiltImage = ImageStack(properties=properties)
                tiltImage.append(ImageReadersRegistry.open(f'{ti.getIndex()}@{ti.getFileName()}'))
                ImageReadersRegistry.write(tiltImage, join(imagesDir, imageName))
                linkedFile = join(averagesDir, imageName)
                # createLink resolves the source against the current directory, not
                # against the link's folder, so it must get the image's real path; it
                # then writes the link as relative (../<imageName>).
                createLink(join(imagesDir, imageName), linkedFile)
                imageNames.append(imageName)
            writeTiltImagesTsStar(ts, tiltImages, imageNames, self._getOutTsStarFile(tsId))
        except Exception as e:
            self._markFailed(tsId, f'Unable to prepare the tilt images and TS star file '
                                   f'with exception {e}')

    def createTsAliFilesStep(self, tsId: str, presentAcqOrders: List[int]):
        """Write the IMOD .xf and .tlt files that ts_import_alignments reads."""
        if tsId in self.failedTsIds:
            return
        try:
            ts = self.seriesDict[tsId]
            self.createImodFiles(ts, presentAcqOrders)
        except Exception as e:
            self._markFailed(tsId, f'Unable to create the IMOD alignment files '
                                   f'with exception {e}')

    def createCtfXmlStep(self, tsId: str, presentAcqOrders: List[int]):
        """Build the Warp per-series CTF XML directly from Scipion's
        previously estimated CTF metadata. ts_import_alignments later
        re-serialises this file to the full schema."""
        if tsId in self.failedTsIds:
            return
        try:
            logger.info(cyanStr(f">>> {tsId} - Writing CTF XML from Scipion CTF..."))
            ts = self.seriesDict[tsId]
            ctf = self.ctfDict[tsId]
            acq = ts.getAcquisition()
            x, y, _ = ts.getDim()
            xmlPath = join(self._getTiltSeriesDir(), f'{tsId}.xml')
            writeCtfXml(xmlPath, ctf, ts, presentAcqOrders,
                        apix=ts.getSamplingRate(),
                        voltage=acq.getVoltage(),
                        cs=acq.getSphericalAberration(),
                        amplitude=acq.getAmplitudeContrast(),
                        imageDimPix=(x, y))
        except Exception as e:
            self._markFailed(tsId, f'Unable to write the CTF XML with exception {e}')

    def createParticlesStarStep(self):
        logger.info(cyanStr(">>> Generating the particles .star file per tomogram..."))
        coordSet = self._getInCoords()
        outPath = self._getParticlesFolder()
        # Iterate only the intersected, pre-cloned tomograms (the input sets may have
        # been subsetted); iterCoordinates(volume=tomo) then yields just that
        # tomogram's coordinates, so no .star references a missing .tomostar / CTF XML.
        for tomoId, tomo in self.tomoDict.items():
            if tomoId in self.failedTsIds:
                continue
            starFile = join(outPath, f'{tomoId}.star')
            try:
                self._writeTomoParticlesStar(coordSet, tomoId, tomo, starFile)
            except Exception as e:
                # Drop the possibly half-written .star, and the .tomostar (_markFailed) so
                # Warp doesn't process the series without its particles either.
                cleanPath(starFile)
                self._markFailed(tomoId, f'Unable to write the particles .star file '
                                         f'with exception {e}')
        if not glob(join(outPath, '*.star')):
            raise Exception('No particles .star file could be written: no tomogram with '
                            'coordinates is left. Check the log for the per-series errors.')

    def tsImportAligmentsStep(self):
        # Failed series have had their .tomostar removed (_markFailed); fail with a clear
        # message if none is left, instead of letting WarpTools run on an empty input.
        if not glob(join(self._getTsStarDir(), f'*{TOMOSTAR_EXT}')):
            raise Exception('No tilt-series could be prepared for Warp. Check the log '
                            'for the per-series errors.')
        logger.info(cyanStr(">>> Starting import alignments..."))
        # ts_import_alignments runs once over every series and takes a single global
        # --alignment_angpix, so it is read from the set.
        args = [
            f"--settings {self._getTiltSeriesSettingsFn()}",
            f"--alignments {self._getTiltStackDir()}",
            f"--alignment_angpix {self._getInTsSet().getSamplingRate()}",
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
        args = " ".join(cmd_parts)
        self.runJob(Plugin.getProgram(WARP_TOOLS, TS_EXPORT_PARTICLES), args, executable='/bin/bash')

    def createOutputStep(self):
        logger.info(cyanStr(">>> Creating the outputs..."))
        coordSet = self._getInCoords()
        # The output is a single RelionSetOfPseudoSubtomograms carrying one sampling
        # rate and one acquisition, so these global values are read from the set.
        tsSet = self._getInTsSet()
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
        if psubtomoSet.getSize() == 0:
            raise Exception('Warp exported no particles. Please check the logs run.stdout '
                            'and run.stderr.')
        outDict = {outputObjects.relionParticles.name: psubtomoSet}
        self._defineOutputs(**outDict)
        self._defineSourceRelation(self.coordinates, psubtomoSet)
        self._defineSourceRelation(self.inCtfSet, psubtomoSet)
        self._defineSourceRelation(self.inTsSet, psubtomoSet)

    # --------------------------- UTILS functions -----------------------------
    def _markFailed(self, tsId: str, msg: str):
        """Exclude a tilt-series from the rest of the run. Recording it in failedTsIds
        only stops the next Scipion steps: WarpTools loads every .tomostar in the
        tomostar folder, so a series with no CTF XML would still be processed (with
        default CTF values, silently). Removing its .tomostar keeps it out of Warp."""
        logger.error(redStr(f'tsId = {tsId} -> {msg}. Skipping... '))
        logger.error(traceback.format_exc())
        self.failedTsIds.append(tsId)
        cleanPath(self._getOutTsStarFile(tsId))

    def _writeTomoParticlesStar(self, coordSet: SetOfCoordinates3D, tomoId: str, tomo,
                                starFile: str):
        """Write the Warp particles .star of one tomogram (normalized coordinates).
        Nothing is written if the tomogram has no coordinates."""
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
            with open(starFile, 'w') as f:
                particlesTable.writeStar(f, tableName='')

    def _getInCoords(self) -> SetOfCoordinates3D:
        return self.coordinates.get()

    def _getInTsSet(self) -> SetOfTiltSeries:
        return self.inTsSet.get()

    def _getInCtfSet(self) -> SetOfCTFTomoSeries:
        return self.inCtfSet.get()

    def createImodFiles(self, ts: TiltSeries, presentAcqOrders: List[int]):
        """Write <tsId>.xf and <tsId>.tlt with one line per tilt of getWarpTiltImages,
        the same tilts and order as the .tomostar and the CTF XML. ts_import_alignments
        assigns the lines to tilts by index and rejects an .xf whose line count differs
        from the tomostar, so Scipion's ts.writeImodFiles (all tilts, including excluded
        ones, in stack order) can't be used here."""
        tsId = ts.getTsId()
        tiltstackFolder = self._getTsTiltStackDir(tsId)
        makePath(tiltstackFolder)
        tiltImages = getWarpTiltImages(ts, presentAcqOrders)
        xfLines = []
        tltLines = []
        missingTransform = False
        for ti in tiltImages:
            transform = ti.getTransform()
            if transform is None:
                missingTransform = True
                xfLines.append('1.0000000 0.0000000 0.0000000 1.0000000 0.000 0.000')
            else:
                m = transform.getMatrix()
                xfLines.append(f'{m[0, 0]:.7f} {m[0, 1]:.7f} {m[1, 0]:.7f} {m[1, 1]:.7f} '
                               f'{m[0, 2]:.3f} {m[1, 2]:.3f}')
            # The .tomostar stores NEGATED tilt angles (writeTiltImagesTsStar), and
            # ts_import_alignments must receive the alignment in the SAME sign convention.
            # If the tomostar negation ever changes, this must change in lockstep.
            tltLines.append(f'{-ti.getTiltAngle()}')
        if missingTransform:
            logger.warning(yellowStr(f'{tsId}: some tilt images have no alignment; '
                                     f'the identity transform was written for them.'))
        with open(join(tiltstackFolder, f'{tsId}.xf'), 'w') as f:
            f.write('\n'.join(xfLines) + '\n')
        with open(join(tiltstackFolder, f'{tsId}.tlt'), 'w') as f:
            f.write('\n'.join(tltLines) + '\n')

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

    def cleanIntermediateResults(self):
        # The exported particles and STAR files don't reference the tilt images
        # (Warp writes rlnTomoTiltSeriesName as dummy.mrc), so they can go.
        logger.info(cyanStr(">>> Cleaning intermediate results..."))
        cleanPath(self._getTiltImagesDir())
