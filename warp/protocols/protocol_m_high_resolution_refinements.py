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
import os
from collections import defaultdict

import emtools.metadata

from pwem.convert.headers import fixVolume
from pwem.objects import VolumeMask
from pyworkflow import BETA
import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from pyworkflow.object import Integer, Set, Float
from pyworkflow.protocol import GPU_LIST

from reliontomo.convert import convert50_tomo, readSetOfPseudoSubtomograms
from reliontomo.objects import createSetOfRelionPSubtomograms
from tomo.objects import CTFTomoSeries, CTFTomo, SetOfCTFTomoSeries,  AverageSubTomogram

from warp.constants import *
from warp.protocols.protocol_base import ProtWarpBase
from warp.utils import updateCtFXMLFile, parseCtfXMLFile, extractGlobalResolution, modifyStarFileMultiTable


class ProtWarpMHigResolutionRefinement(ProtWarpBase):
    """
    Performs high resolution subtomogram refinement for cryo-electron
    tomography datasets, with the goal of improving particle alignment,
    optical parameter estimation, and the final structural quality of
    three-dimensional reconstructions.
    More info:
        https://warpem.github.io/warp/user_guide/warptools/quick_start_warptools_tilt_series/#initial-3d-refinement-in-relion

    AI Generated:

    M High Resolution Refinement (ProtWarpMHigResolutionRefinement) - User Manual
        Overview

        The M High Resolution Refinement protocol is designed for
        high-precision refinement of subtomogram data in cryo-electron
        tomography. Its main objective is to improve the consistency
        between particle observations, tilt-series geometry, and optical
        models so that the resulting three-dimensional reconstruction
        reaches the highest possible biological interpretability.

        In practical cryo-ET workflows, this protocol is typically used
        after an initial subtomogram reconstruction has already been
        obtained. At this stage, the user already has particle
        coordinates, an initial reference, and an approximate structural
        model. The refinement process then improves particle poses,
        microscope parameters, and image distortions in a coordinated
        manner, allowing weak structural details to become more reliable.
        More info:
        https://warpem.github.io/warp/user_guide/warptools/quick_start_warptools_tilt_series/#initial-3d-refinement-in-relion

        Biological Context

        For biological users, the protocol becomes especially valuable
        when studying macromolecular complexes in crowded cellular
        environments, where imperfect alignment or inaccurate optical
        modeling can limit the attainable resolution. By refining these
        parameters jointly, the protocol helps recover more accurate
        molecular features while preserving the native contextual
        information provided by tomography.

        This refinement is commonly used for in situ structural biology,
        membrane protein analysis, viral assemblies, and other systems
        where particles are heterogeneous in orientation but expected to
        share a common structural core.

        Inputs and Starting Conditions

        The protocol requires a tilt-series dataset, particle
        definitions, contrast transfer information, and a reference
        subtomogram average. A corresponding mask is also important
        because it focuses refinement on the structurally meaningful
        region and reduces the influence of solvent, background density,
        or flexible peripheral regions.

        A particularly useful feature is the ability to continue from a
        previous refinement. This enables iterative improvement across
        several rounds of processing, which is often essential in
        high-resolution cryo-ET workflows. Users can therefore either
        start from a newly prepared dataset or continue from a previous
        refinement state without rebuilding the entire processing
        context.

        Species Definition and Biological Relevance

        The molecular diameter defines the approximate physical size of
        the target complex and strongly influences refinement behavior.
        A realistic value helps balance signal usage and noise
        suppression. If the diameter is too small, biologically relevant
        peripheral density may be excluded. If it is too large,
        irrelevant surrounding signal may degrade refinement quality.

        Symmetry can also be specified. This is biologically meaningful
        only when the molecular complex truly obeys a known symmetry.
        Applying incorrect symmetry may artificially sharpen maps while
        introducing structural artifacts. In most uncertain cases,
        asymmetric refinement remains the safer initial strategy.

        The protocol also supports temporal sampling, which is useful
        when particle trajectories or local motion are expected to vary
        over time. This can improve modeling of subtle dynamic behavior
        in challenging datasets.

        Refinement Strategy

        A central goal of this protocol is simultaneous optimization of
        particle alignment and imaging parameters. Particle pose
        refinement improves the relative orientation and position of
        each particle. For most biological applications, this is one of
        the most important contributors to improved map quality.

        Image warp refinement accounts for residual local distortions
        that may remain after initial preprocessing. This is especially
        useful in tilt-series data where local geometric consistency may
        vary across the field of view.

        Stage-angle refinement improves the geometric consistency of the
        tilt-series itself. This can be important when alignment errors
        remain after earlier tomographic preprocessing.

        Optical refinement can further improve reconstruction quality by
        refining defocus, phase shift, spherical aberration, and higher
        order aberration terms. These corrections become increasingly
        important as the target resolution improves. For biological
        users, such refinement is most meaningful when the underlying
        data quality is already good and when sufficient particle
        numbers are available.

        Resolution Control and Reference Conditioning

        The protocol allows optional resampling of the reference maps
        and masks to a different pixel size. This can be useful when
        harmonizing datasets that originate from different processing
        environments or when controlling computational cost during
        intermediate refinement stages.

        An optional low-pass filter may also be applied to the initial
        reference. This is often biologically sensible because it helps
        reduce model bias during early refinement rounds. Starting from
        an overly sharp reference can force incorrect high-frequency
        agreement that does not reflect the actual experimental signal.

        Outputs and Their Interpretation

        After completion, the protocol produces refined particle
        information, updated optical estimates, an improved subtomogram
        average, corresponding half maps, and an updated mask. These
        outputs form the basis for downstream interpretation, validation,
        visualization, or additional iterative refinement.

        The refined average represents the most important biological
        outcome. Improved continuity of secondary structure, clearer
        domain boundaries, and better local contrast are often the
        practical indicators that refinement has been successful.

        A global resolution estimate is also generated. This provides a
        useful overall indicator of refinement progress, but it should
        not be treated as the sole measure of biological quality. Local
        interpretability, preservation of expected features, and
        agreement with independent biological knowledge remain equally
        important.

        Practical Recommendations

        In routine cryo-ET practice, it is usually advisable to begin
        with conservative refinement settings and a reliable reference.
        Early rounds often focus on improving particle poses and defocus
        consistency. More advanced optical refinement is often better
        reserved for later stages, once alignment has stabilized.

        A carefully chosen mask usually has a major impact on success.
        For flexible complexes, focusing refinement on a rigid core
        often improves stability substantially. Overly broad masks may
        dilute signal, whereas overly tight masks may suppress relevant
        structural information.

        Iterative continuation from previous refinements is commonly the
        most effective strategy. Rather than expecting a single run to
        produce the final map, users generally obtain better results by
        progressively refining the dataset over multiple cycles.

        Final Perspective

        For cryo-electron tomography, high resolution refinement is not
        merely a technical optimization step. It is a biologically
        important stage that directly determines how faithfully the
        final reconstruction reflects molecular reality. Careful choice
        of reference, biologically sensible masking, and a gradual
        refinement strategy are the key factors that allow this
        protocol to extract the maximum structural information from
        subtomogram datasets.
    """

    _label = 'M high resolution refinement'
    _devStatus = BETA

    def __init__(self, **kwargs):
        ProtWarpBase.__init__(self, **kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection("Input")
        form.addParam("inputFromMProtocol", params.BooleanParam, default=False,
                      label="Input from M protocol ?",
                      help='Select a M high resolution refinement protocol or an input data')
        form.addParam("inputToMProtocol", params.PointerParam, default=None,
                      pointerClass="ProtWarpMHigResolutionRefinement",
                      label="M protocol",
                      help='M high resolution refinement protocol',
                      condition='inputFromMProtocol')
        form.addParam('inputSet',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      label='Input set of tilt-series',
                      help='Input set of tilt-series')
        form.addParam('inputSetOfCtfTomoSeries',
                      params.PointerParam,
                      condition='not inputFromMProtocol',
                      label="Input CTF estimation",
                      pointerClass='SetOfCTFTomoSeries',
                      help='Select the CTF estimation for the set '
                           'of tilt-series.')
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

        form.addHidden(GPU_LIST, params.StringParam, default='0',
                       label='Choose GPU IDs:', validators=[params.NonEmpty],
                       help='Space-separated list of GPU IDs to use for processing. '
                            'Default: all GPUs in the system For example: "0 1 5"')

    def _insertAllSteps(self):
        self.globalResolution = Float()
        inReParticles = self.getInputSetOfReParticles()
        inputTs = self.getInputSetTS()
        coordSet = inReParticles.getCoordinates3D()
        tsSr = inputTs.getSamplingRate()
        coordSr = coordSet.getSamplingRate()
        tomoDim = coordSet.getPrecedents().getDim()
        scaleFactor = coordSr / tsSr
        self.tomo_thickness = Integer(round(tomoDim[2] * scaleFactor))
        self.x_dimension = Integer(round(tomoDim[0] * scaleFactor))
        self.y_dimension = Integer(round(tomoDim[1] * scaleFactor))

        if not self.inputFromMProtocol.get():
            self._insertFunctionStep(self.prepareDataStep, needsGPU=True)
            self._insertFunctionStep(self.createPopulationStep, needsGPU=False)
            self._insertFunctionStep(self.createSourcesStep, needsGPU=False)
            self._insertFunctionStep(self.createSpeciesStep, needsGPU=True)
            self._insertFunctionStep(self.checkSetup, needsGPU=True)
        else:
            self._insertFunctionStep(self.prepareMDataStep, needsGPU=False)
        self._insertFunctionStep(self.refinementStep, needsGPU=True)
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)
        self._insertFunctionStep(self.closeOutputStep, needsGPU=False)

    def getInputSetTS(self):
        inputSet = self.inputSet.get()
        if self.inputFromMProtocol.get():
            mPrevProt = self.inputToMProtocol.get()
            while mPrevProt.inputFromMProtocol != None:
                mPrevProt = mPrevProt.inputFromMProtocol
            inputSet = mPrevProt.inputSet.get()
        return inputSet

    def getInputSetOfCtfTomoSeries(self):
        inputSetOfCtfTomoSeries = self.inputSetOfCtfTomoSeries.get()
        if self.inputFromMProtocol.get():
            mPrevProt = self.inputToMProtocol.get()
            inputSetOfCtfTomoSeries = mPrevProt.CTFTomoSeries
        return inputSetOfCtfTomoSeries

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

    def prepareMDataStep(self):
        mPrevProt = self.inputToMProtocol.get()
        mPrevProtExtraPath = mPrevProt._getExtraPath()
        pwutils.cleanPath(self._getExtraPath())
        pwutils.createLink(mPrevProtExtraPath, self.getPath('extra'))

    def prepareDataStep(self):
        inReParticles = self.getInputSetOfReParticles()
        outPath = self._getExtraPath()
        writer = convert50_tomo.Writer()
        writer.pseudoSubtomograms2Star(inReParticles, outPath, isWarp=True)

        inputTs = self.getInputSetTS()
        tsSr = inputTs.getSamplingRate()

        self.createTiltSeriesSetting(None)
        for ts in inputTs.iterItems(iterate=False):
            self.tsDataPrepare(ts)
            self.createImodFiles(ts)
        self.tsCtfEstimation(tsSr)
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
        extraPath = self._getExtraPath()
        averageSubTomoHalfMaps = self.averageSubtomogram.get().getHalfMaps(asList=True)
        populationPath = os.path.join(self._getExtraPath(M_RESULT_FOLDER))
        mask = self.refMask.get().getFileName()
        voxelSize = self.averageSubtomogram.get().getSamplingRate()

        argsDict = {
            "--population": os.path.join(populationPath, 'processing.population'),
            "--name": 'processing_species',
            "--diameter": self.diameter.get(),
            "--sym": self.symmetry.get(),
            "--temporal_samples": self.temporal_samples.get(),
            "--half1": averageSubTomoHalfMaps[0],
            "--half2": averageSubTomoHalfMaps[1],
            "--mask": mask,
            "--particles_relion": os.path.join(extraPath, 'inParticles.star'),
            "--lowpass": self.lowpass.get(),
            "--angpix": voxelSize

        }
        if self.angpix_resample.get() is not None:
            argsDict["--angpix_resample"] = self.angpix_resample.get()

        self.runProgram(argsDict, MTOOLS, CREATE_SPECIES)

    def refinementStep(self):
        populationPath = os.path.join(self._getExtraPath('m'))
        self.info(">>> Running M to check setup...")
        self.refinement(populationPath)

    def checkSetup(self):
        populationPath = os.path.join(self._getExtraPath('m'))
        argsDict = {
            "--population": os.path.join(populationPath, 'processing.population'),
            "--iter": 0
        }
        self.runProgram(argsDict, MCORE, None)

    def refinement(self, populationPath):
        msg = ">>>Refinement with 2D Image Warp: \n"
        argsDict = {
            "--population": os.path.join(populationPath, 'processing.population'),
        }

        cmd = ''

        if self.iter.get() is not None:
            argsDict["--iter"] = self.iter.get()

        if self.refine_imagewarp.get() is not None:
            msg += "* Refinement Image Warp \n"
            argsDict["--refine_imagewarp"] = self.refine_imagewarp.get()

        if self.refine_particles.get():
            msg += "* Refinement Particles Poses \n"
            cmd += '--refine_particles '

        if self.refine_mag.get():
            msg += "* Refine anisotropic magnification \n"
            cmd += '--refine_mag '

        if self.refine_stageangles.get():
            msg += "* Refine stage angles (tilt series only) \n"
            cmd += '--refine_stageangles '

        if self.ctf_defocus.get():
            msg += "* Refine defocus using a local search \n"
            cmd += '--ctf_defocus '

        if self.ctf_defocusexhaustive.get():
            msg += ("* Refine defocus using a more exhaustive grid search in the first sub-iteration; "
                    "only works in combination with ctf_defocus \n")
            cmd += '--ctf_defocusexhaustive '

        if self.ctf_phase.get():
            msg += "* Refine phase shift (phase plate data only) \n"
            cmd += '--ctf_phase '

        if self.ctf_cs.get():
            msg += "* Refine spherical aberration, which is also a proxy for pixel size \n"
            cmd += '--ctf_cs '

        if self.ctf_zernike2.get():
            msg += "* Refine Zernike polynomials of 2nd order (slow) \n"
            cmd += '--ctf_zernike2 '

        if self.ctf_zernike3.get():
            msg += "* Refine Zernike polynomials of 3rd order (beam tilt, trefoil – fast) \n"
            cmd += '--ctf_zernike3 '

        if self.ctf_zernike4.get():
            msg += "* Refine Zernike polynomials of 4th order (slow) \n"
            cmd += '--ctf_zernike4 '

        if self.ctf_zernike5.get():
            msg += "* Refine Zernike polynomials of 5th order (fast) \n"
            cmd += '--ctf_zernike5 '

        self.info(msg)
        self.runProgram(argsDict, MCORE, None, othersCmds=cmd)

    def tsCtfEstimation(self, tsSr):
        """CTF estimation"""

        self.info(">>> Generating ctf estimation...")
        settingFile = self._getExtraPath(TILTSERIE_SETTINGS)
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
        inputTs = self.getInputSetTS()
        for ts in inputTs.iterItems():
            ctfTomoSeries = self.getInputSetOfCtfTomoSeries().getItem('_tsId', ts.getTsId())
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

    def createOutputStep(self):
        self.info(">>> Creating outputs...")
        inputTs = self.getInputSetTS()
        for ts in inputTs.iterItems(iterate=False):
            self.createOutputCTF(ts)
        self.createOutputAverage()
        self.createOutputMask()
        self.createOutputParticles()

    def createOutputParticles(self):
        processingFolder = self.getProcessingFolder()

        if not processingFolder:
            raise FileNotFoundError(f"No folder found matching the pattern: {MATCHING_PROCESSING_SPECIES_PATTERN}")
            # Use the first matching folder (you can modify this logic if needed)

        particlesPath = os.path.join(processingFolder, PROCESSING_SPECIES_PARTICLES)
        sourceTable = emtools.metadata.StarFile(particlesPath).getTable('')
        outputFolder = self._getExtraPath(MATCHING_FOLDER)
        pwutils.makePath(outputFolder)

        groupedRows = defaultdict(list)
        for row in sourceTable:
            micrographName = row.wrpSourceName
            groupedRows[micrographName].append(row)

        for micrographName, rows in groupedRows.items():
            cleanName = os.path.splitext(os.path.basename(micrographName))[0]
            outputFilePath = os.path.join(outputFolder, f"{cleanName}.star")

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
                rowValues = list(row._asdict().values())
                rowValues[0] = rowValues[0] / (angpix * self.x_dimension.get())
                rowValues[1] = rowValues[1] / (angpix * self.y_dimension.get())
                rowValues[2] = rowValues[2] / (angpix * self.tomo_thickness.get())
                outputStar._writeRowValues(rowValues)

            outputStar.close()

        self.exportParticles()

        # Output Relion particles
        inReParticles = self.getInputSetOfReParticles()
        coordSet = inReParticles.getCoordinates3D()
        tsSet = self.inputSet.get()
        boxSize = self.AverageSubTomogram.getDim()[0]
        acq = tsSet.getAcquisition()
        relionFolder = self._getExtraPath(RELION_FOLDER)
        are2dStacks = True

        psubtomoSet = createSetOfRelionPSubtomograms(self._getPath(),
                                                     os.path.join(relionFolder, OPTIMISATION_SET_STAR),
                                                     coordSet,
                                                     template='pseudosubtomograms%s.sqlite',
                                                     tsSamplingRate=angpix,
                                                     relionBinning=angpix / tsSet.getSamplingRate(),
                                                     boxSize=boxSize,
                                                     are2dStacks=are2dStacks,
                                                     acquisition=acq)

        modifyStarFileMultiTable(os.path.join(relionFolder, MATCHING_PARTICLES_STAR),
                                 '_rlnImageName', lambda v: self.normalizeParticlesPath(v))
        modifyStarFileMultiTable(os.path.join(relionFolder, MATCHING_TOMOGRAMS_STAR),
                                 '_rlnTomoTiltSeriesName', lambda v: self.normalizeTomogramsPath(v))
        # Fill the set with the generated particles
        readSetOfPseudoSubtomograms(psubtomoSet)
        outDict = {OUTPUT_RELION_PARTICLES: psubtomoSet}
        self._defineOutputs(**outDict)

    def exportParticles(self):
        self.info(">>> Exporting particles...")
        settingFile = self._getExtraPath(TILTSERIE_SETTINGS)
        matchinFolder = self._getExtraPath(MATCHING_FOLDER)
        output = self._getExtraPath(RELION_FOLDER)
        pwutils.makePath(output)
        boxSixe = self.AverageSubTomogram.getDim()[0]
        argsDict = {
            "--settings": os.path.abspath(settingFile),
            "--input_directory": matchinFolder,
            "--input_pattern": "*.star",
            "--output_star": os.path.join(output, 'matching.star'),
            "--output_angpix": self.angpix_resample.get(),
            "--box": boxSixe,
            "--diameter": self.diameter.get(),
        }
        cmd = '--relative_output_paths --normalized_coords --2d'
        self.runProgram(argsDict, WARP_TOOLS, TS_EXPORT_PARTICLES, othersCmds=cmd)

    def createOutputAverage(self):
        processingFolder = self.getProcessingFolder()

        if not processingFolder:
            raise FileNotFoundError(f"No folder found matching the pattern: {MATCHING_PROCESSING_SPECIES_PATTERN}")

        # Output volume
        vol = AverageSubTomogram()
        sourceVolName = os.path.join(processingFolder, PROCESSING_SPECIES_AVERAGE)
        sourceHalf1 = os.path.join(processingFolder, PROCESSING_SPECIES_HALF1)
        sourceHalf2 = os.path.join(processingFolder, PROCESSING_SPECIES_HALF2)
        volName = self._getPath(PROCESSING_SPECIES_AVERAGE)
        half1 = self._getPath(PROCESSING_SPECIES_HALF1)
        half2 = self._getPath(PROCESSING_SPECIES_HALF2)

        pwutils.copyFile(sourceVolName, volName)
        pwutils.copyFile(sourceHalf1, half1)
        pwutils.copyFile(sourceHalf2, half2)

        fixVolume(volName)  # Fix header for xmipp to consider it a volume instead of a stack
        fixVolume(half1)
        fixVolume(half2)
        vol.setFileName(volName)
        sr = self.angpix_resample.get()
        vol.setSamplingRate(sr)
        vol.setHalfMaps([half1, half2])
        self._defineOutputs(**{OUPUT_AVERAGE_SUBTOMOGRAM: vol})
        processingFolder = self.getProcessingFolder()
        speciesFile = os.path.join(processingFolder, PROCESSING_SPECIES_SPECIES)
        self.globalResolution.set(extractGlobalResolution(speciesFile))

    def createOutputMask(self):
        processingFolder = self.getProcessingFolder()
        if not processingFolder:
            raise FileNotFoundError(f"No folder found matching the pattern: {MATCHING_PROCESSING_SPECIES_PATTERN}")

        volMask = VolumeMask()
        sourceMaskFilePath = os.path.join(processingFolder, PROCESSING_SPECIES_MASK)
        maskFilePath = self._getPath(PROCESSING_SPECIES_MASK)
        pwutils.copyFile(sourceMaskFilePath, maskFilePath)
        fixVolume(maskFilePath)
        volMask.setFileName(maskFilePath)
        sr = self.angpix_resample.get()
        volMask.setSamplingRate(sr)
        self._defineOutputs(**{OUTPUT_MASK_SUBTOMOGRAM: volMask})

    def getProcessingFolder(self):
        # Return the first matching folder (you can modify this logic if needed)
        matchingFolder = ''
        mFolder = self._getExtraPath(M_RESULT_FOLDER)
        speciesFolder = os.path.join(mFolder, SPECIES_FOLDER)
        if os.path.exists(speciesFolder):
            # Search for folder matching the pattern "processing_especies_*"
            folderPattern = os.path.join(speciesFolder, MATCHING_PROCESSING_SPECIES_PATTERN)
            matchingFolder = glob.glob(folderPattern)
            matchingFolder = matchingFolder[0]
        return matchingFolder

    def createOutputCTF(self, ts):
        tsId = ts.getTsId()
        processingFolder = os.path.abspath(self._getExtraPath(TILTSERIES_FOLDER))
        psdStack = os.path.join(processingFolder, POWERSPECTRUM_FOLDER, tsId + '.mrc')
        if ts.isEnabled():
            outputSetOfCTFTomoSeries = self.getOutputSetOfCTFTomoSeries(OUTPUT_CTF_SERIE)
            # CTF outputs
            newCTFTomoSeries = CTFTomoSeries(tsId=tsId)
            newCTFTomoSeries.copyInfo(ts)
            newCTFTomoSeries.setTiltSeries(ts)
            outputSetOfCTFTomoSeries.append(newCTFTomoSeries)
            defocusFilePath = os.path.join(processingFolder, ts.getTsId() + '.xml')
            ctfData, gridCtfData = parseCtfXMLFile(defocusFilePath)
            defocusDelta = float(ctfData['DefocusDelta']) * 1e4
            defocusAngle = float(ctfData['DefocusAngle'])

            index = 0
            for ti in ts.iterItems():
                if ti.isEnabled():
                    newCTFTomo = CTFTomo()
                    newCTFTomo.setAcquisitionOrder(ti.getAcquisitionOrder())
                    newCTFTomo.setIndex(index)
                    newCTFTomo.setObjId(index)
                    defocusU = 0
                    defocusV = 0
                    if index in gridCtfData["Nodes"]:
                        defocusU = gridCtfData["Nodes"][index] + defocusDelta
                        defocusV = gridCtfData["Nodes"][index] - defocusAngle
                    newCTFTomo.setDefocusU(defocusU)
                    newCTFTomo.setDefocusV(defocusV)
                    newCTFTomo.setDefocusAngle(defocusAngle)
                    newCTFTomo.setResolution(0)
                    newCTFTomo.setFitQuality(0)
                    newCTFTomo.standardize()
                    newCTFTomo.setPsdFile(f"{index}@" + psdStack)
                    newCTFTomoSeries.append(newCTFTomo)
                    index += 1

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
            tsSet = self.inputSet
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

    def _summary(self):
        summary = []
        if hasattr(self, 'globalResolution') and self.globalResolution.get() is not None:
            summary.append("* Global Resolution: %s Å" % self.globalResolution.get())
        return summary

    def _cleanExtraFiles(self):
        versionsFolder = os.path.join(self.getProcessingFolder(), 'versions')
        pwutils.cleanPath(versionsFolder)
