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
import glob
import logging
import shutil
import traceback
from collections import defaultdict
from os.path import join, abspath, exists, basename, realpath
import emtools.metadata
from pwem.convert.headers import fixVolume
from pwem.objects import VolumeMask
from pyworkflow import BETA
import pyworkflow.protocol.params as params
from pyworkflow.object import Integer, Set, Float
from pyworkflow.protocol import GPU_LIST
from pyworkflow.utils import Message, removeBaseExt, redStr, makePath, copyFile, createLink, cleanPath, cyanStr, \
    yellowStr
from reliontomo.convert import convert50_tomo, readSetOfPseudoSubtomograms
from reliontomo.objects import createSetOfRelionPSubtomograms
from tomo.objects import CTFTomoSeries, CTFTomo, SetOfCTFTomoSeries, AverageSubTomogram, TiltSeries
from warp.constants import MTOOLS, TILTSERIES_SETTINGS, CREATE_POPULATION, CREATE_SOURCE, M_RESULT_FOLDER, \
    CREATE_SPECIES, MCORE, WARP_TOOLS, TILTSERIES_FOLDER, TOMOSTAR_EXT, \
    MATCHING_PROCESSING_SPECIES_PATTERN, PROCESSING_SPECIES_PARTICLES, MATCHING_FOLDER, RLN_COORDINATE_X, \
    RLN_COORDINATE_Y, RLN_COORDINATE_Z, RLN_ANGLE_ROT, RLN_ANGLE_TILT, RLN_ANGLE_PSI, RLN_RAMDOM_SUBSET, \
    RLN_MICROGRAPH_NAME, RLN_SOURCE_HASH, RELION_FOLDER, OPTIMISATION_SET_STAR, MATCHING_PARTICLES_STAR, \
    MATCHING_TOMOGRAMS_STAR, OUTPUT_RELION_PARTICLES, TS_EXPORT_PARTICLES, PROCESSING_SPECIES_AVERAGE, \
    PROCESSING_SPECIES_HALF1, PROCESSING_SPECIES_HALF2, OUPUT_AVERAGE_SUBTOMOGRAM, PROCESSING_SPECIES_SPECIES, \
    PROCESSING_SPECIES_MASK, OUTPUT_MASK_SUBTOMOGRAM, SPECIES_FOLDER, OUTPUT_CTF_SERIES, TOMOSTAR_FOLDER, \
    TILTIMAGES_FOLDER
from warp.protocols.protocol_base_tomo import ProtWarpTomoBase
from warp.utils import parseCtfXMLFile, extractGlobalResolution, modifyStarFileMultiTable, copyWarpSettings, \
    copyWarpTsXml, readWarpSettingsTomoDims, readWarpTsMoviePaths, readWarpTsGeometry, writeWarpTiltImages, \
    getWarpTiltImageName

logger = logging.getLogger(__name__)
PROCESSING_POPULATION = 'processing.population'


class ProtWarpMHigResolutionRefinement(ProtWarpTomoBase):
    """
    High Resolution Refinements in M
    More info:
        https://warpem.github.io/warp/user_guide/warptools/quick_start_warptools_tilt_series/#initial-3d-refinement-in-relion
    """

    _label = 'M high resolution refinement'
    _devStatus = BETA

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(Message.LABEL_INPUT)
        form.addParam("inputFromMProtocol", params.BooleanParam, default=False,
                      label="Input from M protocol ?",
                      help='Select a M high resolution refinement protocol or an input data')

        form.addParam("inputToMProtocol", params.PointerParam, default=None,
                      pointerClass="ProtWarpMHigResolutionRefinement",
                      label="M protocol",
                      help='M high resolution refinement protocol',
                      condition='inputFromMProtocol')

        form.addParam('exportProt',
                      params.PointerParam,
                      pointerClass='ProtWarpExportParticles',
                      condition='not inputFromMProtocol',
                      label='Warp export particles run',
                      help='Run of "Export particles" (Warp) that exported the particles refined '
                           'here. M uses its Warp project: the settings (pixel size and the volume '
                           'the particles were placed in), and each tilt-series XML with the CTF, '
                           'the defocus handedness and the alignment. Its input tilt-series are '
                           'used to write the tilt images M reads.')

        form.addParam('inReParticles', params.PointerParam,
                      condition='not inputFromMProtocol',
                      pointerClass='RelionSetOfPseudoSubtomograms',
                      label='Relion pseudosubtomograms',
                      help='Relion set of pseudoSubtomograms')

        form.addParam('averageSubtomogram', params.PointerParam, pointerClass='Volume',
                      pointerCondition='hasHalfMaps',
                      condition='not inputFromMProtocol',
                      default=None,
                      label='Reference volume',
                      help='Average subtomogram with half maps')

        form.addParam('refMask', params.PointerParam, pointerClass='VolumeMask, Volume',
                      condition='not inputFromMProtocol',
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
        form.addParam('iter', params.IntParam, default=None,
                      allowsNull=True,
                      label='Refinement sub-iterations',
                      help='Number of refinement sub-iterations')

        form.addParam('refine_particles', params.BooleanParam, default=True,
                      label='Refine particle poses',
                      help="Refine particle poses")

        form.addParam('refine_imagewarp', params.StringParam, default=None,
                      allowsNull=True,
                      label='Image warp grid',
                      help="Refine image warp with a grid of XxY dimensions. "
                           "Examples: leave blank = don't refine, '1x1', '6x4'")

        form.addParam('refine_stageangles', params.BooleanParam, default=True,
                      label='Refine stage angles',
                      help="Refine stage angles (tilt series only)")

        form.addParam('refine_mag', params.BooleanParam, default=False,
                      label='Refine anisotropic magnification',
                      help="Refine anisotropic magnification")

        form.addParam('ctf_defocus', params.BooleanParam, default=True,
                      label='Refine defocus',
                      help="Refine defocus using a local search")

        form.addParam('ctf_defocusexhaustive', params.BooleanParam, default=True,
                      label='Refine defocus more exhaustive',
                      help="Refine defocus using a more exhaustive grid search in the first "
                           "sub-iteration; only works in combination with ctf_defocus")

        form.addParam('ctf_phase', params.BooleanParam, default=False,
                      label='Refine phase shift',
                      help="Refine phase shift (phase plate data only)")

        form.addParam('ctf_cs', params.BooleanParam, default=False,
                      label='Refine spherical aberration',
                      help="Refine spherical aberration, which is also a proxy for pixel size")

        form.addParam('ctf_zernike2', params.BooleanParam, default=False,
                      label='Refine Zernike polynomials of 2nd order',
                      help="Refine Zernike polynomials of 2nd order (slow)")

        form.addParam('ctf_zernike3', params.BooleanParam, default=False,
                      label='Refine Zernike polynomials of 3rd order',
                      help="Refine Zernike polynomials of 3rd order (beam tilt, trefoil – fast)")

        form.addParam('ctf_zernike4', params.BooleanParam, default=False,
                      label='Refine Zernike polynomials of 4th order',
                      help="Refine Zernike polynomials of 4th order (slow)")

        form.addParam('ctf_zernike5', params.BooleanParam, default=False,
                      label='Refine Zernike polynomials of 5th order',
                      help="Refine Zernike polynomials of 5th order (fast)")

        # Threads = MCore --perdevice (worker processes per GPU), as in ProtWarpTSMotionCorr.
        # It replaces the --perdevice_refine option of older Warp versions, which MCore no
        # longer accepts. The steps run serially, so the threads aren't used for anything else.
        form.addParallelSection(threads=1, mpi=0)
        form.addHidden(GPU_LIST, params.StringParam, default='0',
                       label='Choose GPU IDs:', validators=[params.NonEmpty],
                       help='Space-separated list of GPU IDs to use for processing. '
                            'Default: all GPUs in the system For example: "0 1 5"')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self.globalResolution = Float()
        # Volume the particles were placed in by the export, in tilt-series pixels. M
        # takes it from the settings copied from the export (prepareProjectStep).
        x, y, z = readWarpSettingsTomoDims(self._getExportProt()._getTiltSeriesSettingsFn())
        self.x_dimension = Integer(x)
        self.y_dimension = Integer(y)
        self.tomo_thickness = Integer(z)

        if not self.inputFromMProtocol.get():
            self._insertFunctionStep(self.prepareProjectStep, needsGPU=False)
            for tsId in self._getExportTsIds():
                self._insertFunctionStep(self.writeTiltImagesStep, tsId, needsGPU=False)
            self._insertFunctionStep(self.writeParticlesStep, needsGPU=False)
            self._insertFunctionStep(self.createPopulationStep, needsGPU=False)
            self._insertFunctionStep(self.createSourcesStep, needsGPU=False)
            self._insertFunctionStep(self.createSpeciesStep, needsGPU=True)
            self._insertFunctionStep(self.checkSetupStep, needsGPU=True)
        else:
            self._insertFunctionStep(self.prepareMDataStep, needsGPU=False)
        self._insertFunctionStep(self.refinementStep, needsGPU=True)
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)
        self._insertFunctionStep(self.closeOutputStep, needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    def prepareProjectStep(self):
        """Set up the Warp project M runs on from the export's one: the settings and,
        per tilt-series, the .tomostar and the XML.

        The XML already holds everything M needs about the series: the Scipion CTF and
        the defocus handedness (written by the export), and the alignment (written by
        ts_import_alignments). It is copied, not linked, because M rewrites it after
        each refinement. Only the paths change; the settings keep the export's pixel
        size and volume size, which the particle positions depend on. The .tomostar is
        only listed by create_source: M reads the tilts from the XML."""
        logger.info(cyanStr(f">>> Preparing Warp project..."))
        exportProt = self._getExportProt()
        tsStarDir = self._getTsStarDir()
        tsDir = self._getTiltSeriesDir()
        makePath(tsStarDir, tsDir, self._getAveragesDir())
        copyWarpSettings(exportProt._getTiltSeriesSettingsFn(), self._getTiltSeriesSettingsFn(),
                         dataFolder=tsStarDir, processingFolder=tsDir)
        for tsId in self._getExportTsIds():
            copyFile(join(exportProt._getTsStarDir(), f'{tsId}{TOMOSTAR_EXT}'),
                             join(tsStarDir, f'{tsId}{TOMOSTAR_EXT}'))
            copyWarpTsXml(join(exportProt._getTiltSeriesDir(), f'{tsId}.xml'),
                          join(tsDir, f'{tsId}.xml'),
                          dataDirectory=tsStarDir)

    def writeTiltImagesStep(self, tsId: str):
        """Write the per-tilt images M reads (tiltimages/ and tiltimages/average/), for
        the tilts of the series' XML, from the export's input tilt-series."""
        logger.info(cyanStr(f">>> {tsId} - Writing the tilt images..."))
        try:
            ts = self._getInTsSet().getItem(TiltSeries.TS_ID_FIELD, tsId)
            tiltImages = self._getWarpTiltImages(ts)
            writeWarpTiltImages(ts, tiltImages, self._getTiltImagesDir(), self._getAveragesDir())
        except Exception as e:
            self._markFailed(tsId, f'Unable to write the tilt images with exception {e}.')

    def writeParticlesStep(self):
        if not self._getProjectTsIds():
            raise Exception('Failed to write the tilt images of all the tilt-series. Check the '
                            'log for the per-series errors.')
        logger.info(cyanStr(">>> Writing the input particles..."))
        writer = convert50_tomo.Writer()
        writer.pseudoSubtomograms2Star(self.getInputSetOfReParticles(), self._getExtraPath(), isWarp=True)

    def createPopulationStep(self):
        logger.info(cyanStr(">>> Creating population..."))
        args = [
            f"--directory {self._getPopulationPath()}",
            "--name processing"
        ]
        self.runProgram(MTOOLS, CREATE_POPULATION, args)

    def createSourcesStep(self):
        logger.info(cyanStr(">>> Creating sources..."))
        args = [
            "--name", "processing",
            "--population", self._getProcessingPopulation(),
            "--processing_settings", self._getTiltSeriesSettingsFn()
        ]
        self.runProgram(MTOOLS, CREATE_SOURCE, args)

    def createSpeciesStep(self):
        logger.info(cyanStr(">>> Creating species..."))
        averageSubTomoHalfMaps = self.averageSubtomogram.get().getHalfMaps(asList=True)
        mask = self.refMask.get().getFileName()
        voxelSize = self.averageSubtomogram.get().getSamplingRate()

        args = [
            f"--population {self._getProcessingPopulation()}",
            "--name processing_species",
            f"--diameter {self.diameter.get()}",
            f"--sym {self.symmetry.get()}",
            f"--temporal_samples {self.temporal_samples.get()}",
            f"--half1 {averageSubTomoHalfMaps[0]}",
            f"--half2 {averageSubTomoHalfMaps[1]}",
            f"--mask {mask}",
            f"--particles_relion {self._getExtraPath('inParticles.star')}",
            f"--angpix {voxelSize}"
        ]
        if self.lowpass.get() is not None:
            args.append(f"--lowpass {self.lowpass.get()}")
        if self.angpix_resample.get() is not None:
            args.append(f"--angpix_resample {self.angpix_resample.get()}")
        failedTsIds = sorted(set(self._getExportTsIds()) - set(self._getProjectTsIds()))
        if failedTsIds:
            # Their particles have no data source in the population (writeTiltImagesStep)
            logger.warning(yellowStr(f'The particles of {", ".join(failedTsIds)} are skipped: '
                                     f'their tilt images could not be written.'))
            args.append("--ignore_unmatched")

        self.runProgram(MTOOLS, CREATE_SPECIES, args, useGpu=True)

    def checkSetupStep(self):
        logger.info(cyanStr(">>> Checking setup..."))
        args = [
            f"--population {self._getProcessingPopulation()}",
            "--iter 0",
            f"--perdevice {self.numberOfThreads.get()}"
        ]
        self.runProgram(MCORE, None, args, useGpu=True)

    def prepareMDataStep(self):
        """Continue from the previous M run: copy its Warp project so that this run
        refines its own population, species and per-series XMLs, and leaves the
        previous run's files (and the outputs built on them) untouched.

        Copied: the population and the species (m/; the population stores the source
        and species paths relative to itself), the source and the per-series XMLs as
        MCore left them (warp_tiltseries/), the .tomostar files and the settings. The
        settings and the XMLs are pointed to this run's folders: ts_export_particles
        reads the settings, and the XMLs locate the tilt images from the .tomostar
        folder (DataDirectory). The tilt images, read-only and as large as the
        tilt-series, are linked, not copied. What the previous run exported (relion/,
        matching/, particle series) is not copied: this run writes its own."""
        prevExtra = realpath(self.inputToMProtocol.get()._getExtraPath())
        logger.info(cyanStr(f">>> Copying the Warp project of the previous M run ({prevExtra})..."))
        # Start clean (also on a restart of this step)
        cleanPath(self._getExtraPath())
        makePath(self._getExtraPath())

        shutil.copytree(join(prevExtra, M_RESULT_FOLDER), self._getPopulationPath(), symlinks=True)

        tsStarDir = self._getTsStarDir()
        tsDir = self._getTiltSeriesDir()
        makePath(tsStarDir, tsDir)
        copyWarpSettings(join(prevExtra, TILTSERIES_SETTINGS), self._getTiltSeriesSettingsFn(),
                         dataFolder=tsStarDir, processingFolder=tsDir)
        for tomoStar in glob.glob(join(prevExtra, TOMOSTAR_FOLDER, f'*{TOMOSTAR_EXT}')):
            copyFile(tomoStar, join(tsStarDir, basename(tomoStar)))
        prevTsDir = join(prevExtra, TILTSERIES_FOLDER)
        for xmlFile in glob.glob(join(prevTsDir, '*.xml')):
            copyWarpTsXml(xmlFile, join(tsDir, basename(xmlFile)), dataDirectory=tsStarDir)
        for sourceFile in glob.glob(join(prevTsDir, '*.source')):
            copyFile(sourceFile, join(tsDir, basename(sourceFile)))
        # Earlier versions of the source (XML copies), kept by M's versioning
        prevVersions = join(prevTsDir, 'versions')
        if exists(prevVersions):
            shutil.copytree(prevVersions, join(tsDir, 'versions'), symlinks=True)

        # realpath: a chain of continuations links straight to the first run's images
        createLink(realpath(join(prevExtra, TILTIMAGES_FOLDER)), self._getTiltImagesDir())

    def refinementStep(self):
        logger.info(cyanStr(">>> Running refinement with M..."))
        self.refinement(self._getPopulationPath())

    def createOutputStep(self):
        logger.info(cyanStr(">>> Creating the outputs..."))
        for ts in self._getInTsSet().iterItems(iterate=False):
            self.createOutputCTF(ts)
        self.createOutputAverage()
        self.createOutputMask()
        self.createOutputParticles()

    # --------------------------- UTILS functions -----------------------------
    def _getPopulationPath(self) -> str:
        return self._getExtraPath(M_RESULT_FOLDER)

    def _getProcessingPopulation(self) -> str:
        return join(self._getPopulationPath(), PROCESSING_POPULATION)

    def _getExportProt(self):
        """Export particles run of the first M run of the chain."""
        prot = self
        while prot.inputFromMProtocol.get():
            prot = prot.inputToMProtocol.get()
        return prot.exportProt.get()

    def _getInTsSetPointer(self):
        return self._getExportProt().inTsSet

    def _getInTsSet(self):
        return self._getInTsSetPointer().get()

    def _getExportTsIds(self):
        """Series the export prepared for Warp: those with a .tomostar and an XML (the
        export removes the .tomostar of a series that failed)."""
        exportProt = self._getExportProt()
        tsIds = []
        for tomoStar in sorted(glob.glob(join(exportProt._getTsStarDir(), f'*{TOMOSTAR_EXT}'))):
            tsId = basename(tomoStar)[:-len(TOMOSTAR_EXT)]
            if exists(join(exportProt._getTiltSeriesDir(), f'{tsId}.xml')):
                tsIds.append(tsId)
        return tsIds

    def _getProjectTsIds(self):
        """Series of M's own Warp project (prepareProjectStep): those with a .tomostar
        and an XML. A series whose tilt images could not be written has been removed
        from it (_markFailed)."""
        tsIds = []
        for tomoStar in sorted(glob.glob(join(self._getTsStarDir(), f'*{TOMOSTAR_EXT}'))):
            tsId = basename(tomoStar)[:-len(TOMOSTAR_EXT)]
            if exists(join(self._getTiltSeriesDir(), f'{tsId}.xml')):
                tsIds.append(tsId)
        return tsIds

    def _getWarpTiltImages(self, ts):
        """Tilt images of ts in the series' XML (MoviePath), in the same order: the
        tilts and order Warp indexes the series with. The export wrote the XML, and
        named each tilt image after its index in the tilt-series stack."""
        tsId = ts.getTsId()
        xmlPath = join(self._getTiltSeriesDir(), f'{tsId}.xml')
        names = [basename(path) for path in readWarpTsMoviePaths(xmlPath)]
        tiByName = {getWarpTiltImageName(tsId, ti): ti.clone() for ti in ts.iterItems()}
        missing = [name for name in names if name not in tiByName]
        if missing:
            raise Exception(f'{tsId}: tilt images {missing} of {xmlPath} are not in the input '
                            f'tilt-series of the export run.')
        tiltImages = [tiByName[name] for name in names]
        # The XML stores the negated tilt angles (as the .tomostar)
        xmlAngles = readWarpTsGeometry(xmlPath)['angles']
        for ti, angle in zip(tiltImages, xmlAngles):
            if abs(-ti.getTiltAngle() - angle) > 0.05:
                raise Exception(f'{tsId}: the tilt angle of image {ti.getIndex()} '
                                f'({ti.getTiltAngle()}) does not match {xmlPath} ({-angle}).')
        return tiltImages

    def getInputSetOfReParticles(self):
        inReParticles = self.inReParticles.get()
        if self.inputFromMProtocol.get():
            mPrevProt = self.inputToMProtocol.get()
            inReParticles = mPrevProt.RelionParticles
        return inReParticles

    def getInputAverageSubtomogram(self):
        averageSubtomogram = self.averageSubtomogram.get()
        if self.inputFromMProtocol.get():
            mPrevProt = self.inputToMProtocol.get()
            averageSubtomogram = mPrevProt.AverageSubTomogram
        return averageSubtomogram

    def getRefMask(self):
        refMask = self.refMask.get()
        if self.inputFromMProtocol.get():
            mPrevProt = self.inputToMProtocol.get()
            refMask = mPrevProt.MaskSubTomogram
        return refMask

    def refinement(self, populationPath):
        msg = ">>>Refinement with 2D Image Warp: \n"
        cmd = [
            f"--population {join(populationPath, 'processing.population')}"
        ]

        if self.iter.get() is not None:
            cmd.append(f"--iter {self.iter.get()}")

        if self.refine_imagewarp.get():  # None or empty: don't refine
            msg += "* Refinement Image Warp\n"
            cmd.append(f"--refine_imagewarp {self.refine_imagewarp.get()}")

        if self.refine_particles.get():
            msg += "* Refinement Particles Poses\n"
            cmd.append("--refine_particles")

        if self.refine_mag.get():
            msg += "* Refine anisotropic magnification\n"
            cmd.append("--refine_mag")

        if self.refine_stageangles.get():
            msg += "* Refine stage angles (tilt series only)\n"
            cmd.append("--refine_stageangles")

        if self.ctf_defocus.get():
            msg += "* Refine defocus using a local search\n"
            cmd.append("--ctf_defocus")

        if self.ctf_defocusexhaustive.get():
            msg += ("* Refine defocus using a more exhaustive grid search in the first "
                    "sub-iteration; only works in combination with ctf_defocus\n")
            cmd.append("--ctf_defocusexhaustive")

        if self.ctf_phase.get():
            msg += "* Refine phase shift (phase plate data only)\n"
            cmd.append("--ctf_phase")

        if self.ctf_cs.get():
            msg += "* Refine spherical aberration, which is also a proxy for pixel size\n"
            cmd.append("--ctf_cs")

        if self.ctf_zernike2.get():
            msg += "* Refine Zernike polynomials of 2nd order (slow)\n"
            cmd.append("--ctf_zernike2")

        if self.ctf_zernike3.get():
            msg += "* Refine Zernike polynomials of 3rd order (beam tilt, trefoil - fast)\n"
            cmd.append("--ctf_zernike3")

        if self.ctf_zernike4.get():
            msg += "* Refine Zernike polynomials of 4th order (slow)\n"
            cmd.append("--ctf_zernike4")

        if self.ctf_zernike5.get():
            msg += "* Refine Zernike polynomials of 5th order (fast)\n"
            cmd.append("--ctf_zernike5")

        cmd.append(f"--perdevice {self.numberOfThreads.get()}")
        logger.info(cyanStr(msg))
        self.runProgram(MCORE, None, cmd, useGpu=True)

    @staticmethod
    def _getTemporalMean(rowData, prefix):
        samples = [
            (int(key[len(prefix):]), float(value))
            for key, value in rowData.items()
            if key.startswith(prefix) and key[len(prefix):].isdigit()
        ]

        if not samples:
            raise ValueError(f'No temporal samples found for {prefix}')

        samples.sort(key=lambda item: item[0])
        return sum(value for _, value in samples) / len(samples)

    def createOutputParticles(self):
        processingFolder = self.getProcessingFolder()

        if not processingFolder:
            raise FileNotFoundError(f"No folder found matching the pattern: {MATCHING_PROCESSING_SPECIES_PATTERN}")
            # Use the first matching folder (you can modify this logic if needed)

        particlesPath = join(processingFolder, PROCESSING_SPECIES_PARTICLES)
        sourceTable = emtools.metadata.StarFile(particlesPath).getTable('')
        outputFolder = self._getExtraPath(MATCHING_FOLDER)
        makePath(outputFolder)

        groupedRows = defaultdict(list)
        for row in sourceTable:
            micrographName = row.wrpSourceName
            groupedRows[micrographName].append(row)

        for micrographName, rows in groupedRows.items():
            cleanName = removeBaseExt(micrographName)
            outputFilePath = join(outputFolder, f"{cleanName}.star")

            columnNames = [
                RLN_COORDINATE_X,
                RLN_COORDINATE_Y,
                RLN_COORDINATE_Z,
                RLN_ANGLE_ROT,
                RLN_ANGLE_TILT,
                RLN_ANGLE_PSI,
                RLN_RAMDOM_SUBSET,
                RLN_MICROGRAPH_NAME,
                RLN_SOURCE_HASH
            ]

            outputTable = emtools.metadata.Table(columns=columnNames)
            outputStar = emtools.metadata.StarFile(outputFilePath, 'w')
            outputStar.writeHeader('', outputTable)
            angpix = self.angpix_resample.get()

            for row in rows:
                rowData = row._asdict()

                coordX = self._getTemporalMean(rowData, 'wrpCoordinateX')
                coordY = self._getTemporalMean(rowData, 'wrpCoordinateY')
                coordZ = self._getTemporalMean(rowData, 'wrpCoordinateZ')

                angleRot = self._getTemporalMean(rowData, 'wrpAngleRot')
                angleTilt = self._getTemporalMean(rowData, 'wrpAngleTilt')
                anglePsi = self._getTemporalMean(rowData, 'wrpAnglePsi')

                rowValues = [
                    coordX / (angpix * self.x_dimension.get()),
                    coordY / (angpix * self.y_dimension.get()),
                    coordZ / (angpix * self.tomo_thickness.get()),
                    angleRot,
                    angleTilt,
                    anglePsi,
                    rowData['wrpRandomSubset'],
                    rowData['wrpSourceName'],
                    rowData['wrpSourceHash']
                ]

                outputStar._writeRowValues(rowValues)

            outputStar.close()

        self.exportParticles()

        # Output Relion particles
        inReParticles = self.getInputSetOfReParticles()
        coordSet = inReParticles.getCoordinates3D()
        tsSet = self._getInTsSet()
        boxSize = self.AverageSubTomogram.getDim()[0]
        acq = tsSet.getAcquisition()
        relionFolder = self._getExtraPath(RELION_FOLDER)
        are2dStacks = True

        psubtomoSet = createSetOfRelionPSubtomograms(self._getPath(),
                                                     join(relionFolder, OPTIMISATION_SET_STAR),
                                                     coordSet,
                                                     template='pseudosubtomograms%s.sqlite',
                                                     tsSamplingRate=angpix,
                                                     relionBinning=angpix / tsSet.getSamplingRate(),
                                                     boxSize=boxSize,
                                                     are2dStacks=are2dStacks,
                                                     acquisition=acq)

        modifyStarFileMultiTable(join(relionFolder, MATCHING_PARTICLES_STAR),
                                 '_rlnImageName', lambda v: self.normalizeParticlesPath(v))
        modifyStarFileMultiTable(join(relionFolder, MATCHING_TOMOGRAMS_STAR),
                                 '_rlnTomoTiltSeriesName', lambda v: self.normalizeTomogramsPath(v))
        # Fill the set with the generated particles
        readSetOfPseudoSubtomograms(psubtomoSet)
        outDict = {OUTPUT_RELION_PARTICLES: psubtomoSet}
        self._defineOutputs(**outDict)

    def exportParticles(self):
        logger.info(cyanStr(">>> Exporting particles..."))
        settingFile = self._getExtraPath(TILTSERIES_SETTINGS)
        matchinFolder = self._getExtraPath(MATCHING_FOLDER)
        output = self._getExtraPath(RELION_FOLDER)
        makePath(output)
        boxSixe = self.AverageSubTomogram.getDim()[0]
        cmd = [
            f"--settings {abspath(settingFile)}",
            f"--input_directory {matchinFolder}",
            "--input_pattern *.star",
            f"--output_star {join(output, 'matching.star')}",
            f"--output_angpix {self.angpix_resample.get()}",
            f"--box {boxSixe}",
            f"--diameter {self.diameter.get()}",
            "--relative_output_paths",
            "--normalized_coords",
            "--2d"
        ]
        self.runProgram(WARP_TOOLS, TS_EXPORT_PARTICLES, cmd, useGpu=True)

    def createOutputAverage(self):
        processingFolder = self.getProcessingFolder()

        if not processingFolder:
            raise FileNotFoundError(f"No folder found matching the pattern: {MATCHING_PROCESSING_SPECIES_PATTERN}")

        # Output volume
        vol = AverageSubTomogram()
        sourceVolName = join(processingFolder, PROCESSING_SPECIES_AVERAGE)
        sourceHalf1 = join(processingFolder, PROCESSING_SPECIES_HALF1)
        sourceHalf2 = join(processingFolder, PROCESSING_SPECIES_HALF2)
        volName = self._getPath(PROCESSING_SPECIES_AVERAGE)
        half1 = self._getPath(PROCESSING_SPECIES_HALF1)
        half2 = self._getPath(PROCESSING_SPECIES_HALF2)

        copyFile(sourceVolName, volName)
        copyFile(sourceHalf1, half1)
        copyFile(sourceHalf2, half2)

        fixVolume(volName)  # Fix header for xmipp to consider it a volume instead of a stack
        fixVolume(half1)
        fixVolume(half2)
        vol.setFileName(volName)
        sr = self.angpix_resample.get()
        vol.setSamplingRate(sr)
        vol.setHalfMaps([half1, half2])
        self._defineOutputs(**{OUPUT_AVERAGE_SUBTOMOGRAM: vol})
        processingFolder = self.getProcessingFolder()
        speciesFile = join(processingFolder, PROCESSING_SPECIES_SPECIES)
        self.globalResolution.set(extractGlobalResolution(speciesFile))

    def createOutputMask(self):
        processingFolder = self.getProcessingFolder()
        if not processingFolder:
            raise FileNotFoundError(f"No folder found matching the pattern: {MATCHING_PROCESSING_SPECIES_PATTERN}")

        volMask = VolumeMask()
        sourceMaskFilePath = join(processingFolder, PROCESSING_SPECIES_MASK)
        maskFilePath = self._getPath(PROCESSING_SPECIES_MASK)
        copyFile(sourceMaskFilePath, maskFilePath)
        fixVolume(maskFilePath)
        volMask.setFileName(maskFilePath)
        sr = self.angpix_resample.get()
        volMask.setSamplingRate(sr)
        self._defineOutputs(**{OUTPUT_MASK_SUBTOMOGRAM: volMask})

    def getProcessingFolder(self):
        # Return the first matching folder (you can modify this logic if needed)
        matchingFolder = ''
        mFolder = self._getExtraPath(M_RESULT_FOLDER)
        speciesFolder = join(mFolder, SPECIES_FOLDER)
        if exists(speciesFolder):
            # Search for folder matching the pattern "processing_especies_*"
            folderPattern = join(speciesFolder, MATCHING_PROCESSING_SPECIES_PATTERN)
            matchingFolder = glob.glob(folderPattern)
            matchingFolder = matchingFolder[0]
        return matchingFolder

    def createOutputCTF(self, ts):
        tsId = ts.getTsId()
        processingFolder = abspath(self._getExtraPath(TILTSERIES_FOLDER))
        defocusFilePath = join(processingFolder, tsId + '.xml')

        if ts.isEnabled() and exists(defocusFilePath):
            outputSetOfCTFTomoSeries = self.getOutputSetOfCTFTomoSeries(OUTPUT_CTF_SERIES)

            newCTFTomoSeries = CTFTomoSeries(tsId=tsId)
            newCTFTomoSeries.copyInfo(ts)
            newCTFTomoSeries.setTiltSeries(ts)
            outputSetOfCTFTomoSeries.append(newCTFTomoSeries)

            ctfData, gridCtfData = parseCtfXMLFile(defocusFilePath)

            defaultDefocusDelta = float(ctfData['DefocusDelta']) * 1e4
            defaultDefocusAngle = float(ctfData['DefocusAngle'])

            # Grid nodes are indexed like the tilts of the XML
            tiltImages = self._getWarpTiltImages(ts)

            for index, ti in enumerate(tiltImages):
                if index not in gridCtfData["Nodes"]:
                    raise ValueError(
                        f'No Warp CTF defocus found for tilt index {index} '
                        f'in tilt-series {tsId}'
                    )

                defocus = gridCtfData["Nodes"][index]
                defocusDelta = gridCtfData["DeltaNodes"].get(
                    index, defaultDefocusDelta
                )
                defocusAngle = gridCtfData["AngleNodes"].get(
                    index, defaultDefocusAngle
                )

                defocusU = defocus + defocusDelta
                defocusV = defocus - defocusDelta

                itemIndex = index + 1

                newCTFTomo = CTFTomo()
                newCTFTomo.setAcquisitionOrder(ti.getAcquisitionOrder())
                newCTFTomo.setIndex(itemIndex)
                newCTFTomo.setDefocusU(defocusU)
                newCTFTomo.setDefocusV(defocusV)
                newCTFTomo.setDefocusAngle(defocusAngle)
                newCTFTomo.setResolution(0)
                newCTFTomo.setFitQuality(0)
                newCTFTomo.standardize()

                newCTFTomoSeries.append(newCTFTomo)

            outputSetOfCTFTomoSeries.update(newCTFTomoSeries)
            outputSetOfCTFTomoSeries.write()
            self._store(outputSetOfCTFTomoSeries)

    def getOutputSetOfCTFTomoSeries(self, outputSetName):
        outputSetOfCTFTomoSeries = getattr(self, outputSetName, None)

        if outputSetOfCTFTomoSeries:
            outputSetOfCTFTomoSeries.enableAppend()
        else:
            outputSetOfCTFTomoSeries = SetOfCTFTomoSeries.create(self._getPath(),
                                                                 template='CTFmodels%s.sqlite')
            tsSet = self._getInTsSetPointer()
            outputSetOfCTFTomoSeries.setSetOfTiltSeries(tsSet)
            outputSetOfCTFTomoSeries.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(**{outputSetName: outputSetOfCTFTomoSeries})
            self._defineCtfRelation(outputSetOfCTFTomoSeries, tsSet)

        return outputSetOfCTFTomoSeries

    def closeOutputStep(self):
        """Finalizes and closes the output sets. """
        self._closeOutputSet()

    def allowsDelete(self, obj):
        return True

    def _cleanExtraFiles(self):
        versionsFolder = join(self.getProcessingFolder(), 'versions')
        cleanPath(versionsFolder)

    def _markFailed(self, tsId: str, msg: str):
        """Drop a tilt-series from M's Warp project. create_source loads every
        .tomostar with an XML, so the series is removed from disk, together with any
        tilt image already written. That also keeps the failure across a resumed run
        (_getProjectTsIds)."""
        logger.error(redStr(f'tsId = {tsId} -> {msg}. Skipping... '))
        logger.error(traceback.format_exc())
        self.failedTsIds.append(tsId)
        cleanPath(join(self._getTsStarDir(), f'{tsId}{TOMOSTAR_EXT}'),
                  join(self._getTiltSeriesDir(), f'{tsId}.xml'),
                  *glob.glob(join(self._getAveragesDir(), f'{tsId}_*.mrc')),
                  *glob.glob(join(self._getTiltImagesDir(), f'{tsId}_*.mrc')))

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []
        if self.inputFromMProtocol.get():
            return errors
        exportProt = self.exportProt.get()
        if exportProt is None:
            return errors  # Already reported by Scipion's input validation
        if not exists(exportProt._getTiltSeriesSettingsFn()):
            errors.append('The Warp export particles run has no Warp settings file.')
        elif not self._getExportTsIds():
            errors.append('The Warp export particles run has no tilt-series prepared for Warp.')
        return errors

    def _summary(self):
        summary = []
        if hasattr(self, 'globalResolution') and self.globalResolution.get() is not None:
            summary.append("* Global Resolution: %s Å" % self.globalResolution.get())
        return summary
