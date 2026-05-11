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
from enum import Enum

from emtable import Table

from pyworkflow import BETA
import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from pyworkflow.object import Integer
from pyworkflow.utils import Message
from reliontomo.convert import readSetOfPseudoSubtomograms
from reliontomo.objects import createSetOfRelionPSubtomograms, RelionSetOfPseudoSubtomograms
from tomo.constants import BOTTOM_LEFT_CORNER

from warp.constants import *
from warp.protocols.protocol_base import ProtWarpBase
from warp.utils import updateCtFXMLFile, getTransformInfoFromCoordOrSubtomo, modifyStarFileMultiTable, \
    modifyOptFileMultiTable


class outputObjects(Enum):
    relionParticles = RelionSetOfPseudoSubtomograms


class ProtWarpExportParticles(ProtWarpBase):
    """
    Exports particles from tilt-series data as either reconstructed 3D
    subtomograms or 2D particle image series, enabling downstream
    structural analysis from localized coordinates in tomographic
    experiments.
    More info:
       https://warpem.github.io/warp/user_guide/warptools/quick_start_warptools_tilt_series/#export-particles

    AI Generated:

    Export Particles (ProtWarpExportParticles) - User Manual
        Overview

        The Export Particles protocol prepares localized particles from
        cryo-electron tomography experiments for downstream analysis by
        converting previously identified particle coordinates into
        standardized particle datasets. Its main objective is to bridge
        the gap between tomographic localization and particle-based
        structural interpretation, allowing biologically meaningful
        regions detected inside reconstructed cellular or purified
        specimens to be extracted in a form suitable for refinement,
        classification, or visualization.

        In practical cryo-ET workflows, this protocol becomes relevant
        once particles have already been detected in tomograms and the
        user wishes to isolate them as individual analysis units. The
        exported particles can represent either volumetric subtomograms
        or aligned 2D particle series, depending on the intended
        downstream strategy. This flexibility makes the protocol useful
        both for subtomogram averaging pipelines and for workflows that
        rely on projection-based particle analysis.

        Inputs and Biological Context

        The protocol requires three biologically linked sources of
        information. First, it uses a set of three-dimensional particle
        coordinates that define the positions of the particles of
        interest inside the tomographic volume. Second, it requires the
        original tilt series used to reconstruct those tomograms, since
        these images contain the experimental signal from which particle
        information is ultimately derived. Third, it uses contrast
        transfer function estimations that provide optical correction
        parameters necessary for physically meaningful reconstruction.

        The biological quality of the exported particles depends heavily
        on the quality of these inputs. Accurate coordinate placement is
        essential because misplaced coordinates may isolate background
        density, neighboring complexes, or incomplete particles. In the
        same way, poorly aligned tilt series or inaccurate CTF
        estimation can degrade the interpretability of the final
        particle data.

        Reconstruction Strategy

        A central role of this protocol is to transform localized
        coordinates into particle-centered datasets with a controlled
        spatial sampling. The user defines the desired output pixel
        size, particle box dimensions, and an approximate particle
        diameter. These parameters together determine the physical scale
        and the amount of surrounding structural context preserved in
        the exported particle.

        From a biological perspective, choosing these values carefully
        is important. A box that is too small may truncate peripheral
        domains, flexible regions, or interaction partners. A box that
        is too large may include excessive solvent or neighboring
        densities that complicate later classification. Similarly, the
        output pixel size should balance computational efficiency
        against the structural detail needed for the intended analysis.

        Choosing Between 2D and 3D Export

        The protocol supports two conceptually different export modes.
        In the 2D mode, particles are represented as image series
        centered on each target location. This option is useful when
        users intend to preserve projection information or apply
        particle-based analyses closer to single-particle workflows.

        In the 3D mode, particles are exported as subtomograms. This is
        generally the preferred option when the biological goal is
        subtomogram averaging, structural classification, or local
        volumetric interpretation of macromolecular complexes in situ.

        The choice between these two modes should reflect the biological
        question being addressed. When studying native macromolecular
        organization inside cells, 3D export is often more directly
        informative. When emphasizing projection consistency or
        experimental image-space analysis, 2D export may be more
        appropriate.

        Coordinate Consistency and Geometric Interpretation

        An important biological feature of the protocol is the
        preservation of particle localization within the tomographic
        coordinate system. Each exported particle remains linked to its
        original spatial context through its position and orientation.
        This is particularly valuable in cryo-ET studies where spatial
        organization itself carries biological meaning, such as
        membrane-associated assemblies, intracellular filament systems,
        or molecular complexes arranged within crowded environments.

        Because of this, reliable coordinate definitions and consistent
        tomogram geometry are critical. Errors in scaling, mismatched
        sampling rates, or inaccurate tomogram references may produce
        particles that are mathematically valid but biologically
        misleading.

        Output Products and Their Interpretation

        After completion, the protocol generates a particle dataset
        compatible with downstream subtomogram-oriented processing
        environments. The resulting output preserves the connection
        between each exported particle, its originating coordinate, and
        its parent tilt-series information.

        For biological users, this means the exported particles are not
        merely cropped image regions. They represent structured
        particle-centered observations that can now be classified,
        aligned, averaged, or interpreted as localized molecular
        instances inside the original specimen.

        Practical Recommendations

        In most practical cryo-ET projects, it is advisable to choose
        an output sampling close to the intended refinement scale and a
        box size that comfortably encloses the expected particle
        diameter together with a moderate surrounding margin. This
        generally provides a useful balance between computational cost
        and structural completeness.

        Before large-scale export, users often benefit from verifying
        that the coordinate set is biologically sensible and that the
        selected particles correspond to recognizable molecular
        features. Small validation runs can prevent costly downstream
        processing of poorly centered or biologically irrelevant
        particle sets.

        Final Perspective

        For cryo-electron tomography users, particle export is not only
        a preparatory conversion step but an important biological
        transition from spatial localization to particle-centered
        structural analysis. When coordinates, sampling, and particle
        dimensions are chosen carefully, the resulting exported dataset
        becomes a reliable starting point for extracting meaningful
        molecular information from complex three-dimensional biological
        environments.
    """

    _label = 'Export particles'
    _devStatus = BETA

    def __init__(self, **kwargs):
        ProtWarpBase.__init__(self, **kwargs)

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

        form.addParam('writeStacks', params.EnumParam,
                      label='Export type',
                      default=0,
                      choices=['2D', '3D'],
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='If set to 2D, this program will write output particles as 2d image series centered on '
                           'the particle (particle series). If set to 3D, this program will write output '
                           'particles as 3d images (subtomograms)')

    def _insertAllSteps(self):
        self._insertFunctionStep(self.prepareDataStep, needsGPU=True)
        self._insertFunctionStep(self.exportParticlesStep, needsGPU=True)
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)
        self._insertFunctionStep(self.cleanIntermediateResults, needsGPU=False)

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
        self.tsCtfEstimation()
        self.updateCTFValues()
        self.tsImportAligments()

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

            angles, _ = getTransformInfoFromCoordOrSubtomo(coord, sRate)
            particlesTable.addRow(
                coord.getX(BOTTOM_LEFT_CORNER) / tomoDim[0],
                coord.getY(BOTTOM_LEFT_CORNER) / tomoDim[1],
                coord.getZ(BOTTOM_LEFT_CORNER) / tomoDim[2],
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
        settingFile = self._getExtraPath(TILTSERIE_SETTINGS)
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
        cmd = '--relative_output_paths --normalized_coords'
        if self.writeStacks.get() == 0:
            cmd += ' --2d'
        else:
            cmd += ' --3d'
        self.runProgram(argsDict, WARP_TOOLS, TS_EXPORT_PARTICLES, othersCmds=cmd)

    def createOutputStep(self):
        coordSet = self.coordinates.get()
        tsSet = self.inputSet.get()
        tsSRate = tsSet.getSamplingRate()
        boxSize = self.box.get()
        acq = tsSet.getAcquisition()
        relionFolder = self._getExtraPath(RELION_FOLDER)
        are2dStacks = self.writeStacks.get() == 0
        modifyOptFileMultiTable(os.path.join(relionFolder, OPTIMISATION_SET_STAR),
                                 '_rlnTomoParticlesFile', lambda v: self.normalizeOptimizationPath(v))
        modifyOptFileMultiTable(os.path.join(relionFolder, OPTIMISATION_SET_STAR),
                                 '_rlnTomoTomogramsFile', lambda v: self.normalizeOptimizationPath(v))
        psubtomoSet = createSetOfRelionPSubtomograms(self._getPath(),
                                                     os.path.join(relionFolder, OPTIMISATION_SET_STAR),
                                                     coordSet,
                                                     template='pseudosubtomograms%s.sqlite',
                                                     tsSamplingRate=tsSRate,
                                                     relionBinning=self.output_angpix.get()/tsSet.getSamplingRate(),
                                                     boxSize=boxSize,
                                                     are2dStacks=are2dStacks,
                                                     acquisition=acq)

        modifyStarFileMultiTable(os.path.join(relionFolder, MATCHING_PARTICLES_STAR),
                                      '_rlnImageName', lambda v: self.normalizeParticlesPath(v))
        modifyStarFileMultiTable(os.path.join(relionFolder, MATCHING_TOMOGRAMS_STAR),
                                      '_rlnTomoTiltSeriesName', lambda v: self.normalizeTomogramsPath(v))
        # Fill the set with the generated particles
        readSetOfPseudoSubtomograms(psubtomoSet)
        outDict = {outputObjects.relionParticles.name: psubtomoSet}
        self._defineOutputs(**outDict)
        self._defineSourceRelation(self.coordinates, psubtomoSet)
        self._defineSourceRelation(self.inputSetOfCtfTomoSeries, psubtomoSet)
        self._defineSourceRelation(self.inputSet, psubtomoSet)

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
        self.info(">>> Starting import alignments...")
        processingFolder = os.path.abspath(self._getExtraPath(TILTSERIES_FOLDER))
        tiltstackFolder = os.path.join(processingFolder, 'tiltstack')
        angpix = self.inputSet.get().getSamplingRate()
        settingFile = self._getExtraPath(TILTSERIE_SETTINGS)
        argsDict = {
            "--settings": os.path.abspath(settingFile),
            '--alignments': os.path.abspath(tiltstackFolder),
            "--alignment_angpix": angpix,
        }
        cmd = ' '.join(['%s %s' % (k, v) for k, v in argsDict.items()])
        self.runJob(self.getPlugin().getProgram(WARP_TOOLS, TS_IMPORT_ALIGNMENTS), cmd, executable='/bin/bash')

    def createImodFiles(self, ts):
        processingFolder = os.path.abspath(self._getExtraPath(TILTSERIES_FOLDER))
        tiltstackFolder = os.path.join(processingFolder, 'tiltstack', ts.getTsId())
        pwutils.makePath(tiltstackFolder)
        ts.writeImodFiles(tiltstackFolder, delimiter=' ', factor=1)

    def cleanIntermediateResults(self):
        self.info(">>> Cleaning intermediate results...")
        imagesFolder = self._getExtraPath(TILTIMAGES_FOLDER)
        pwutils.cleanPath(imagesFolder)




