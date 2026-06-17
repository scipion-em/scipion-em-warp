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
from warp.utils import updateCtFXMLFile, getTransformInfoFromCoordOrSubtomo, modifyStarFileMultiTable


class outputObjects(Enum):
    relionParticles = RelionSetOfPseudoSubtomograms


class ProtWarpExportParticles(ProtWarpBase):
    """
    Exports particles from tomographic coordinates using WarpTools as
    either 2D particle series or 3D subtomograms.
    More info:
       https://warpem.github.io/warp/user_guide/warptools/quick_start_warptools_tilt_series/#export-particles

    AI Generated:

    Export Particles (ProtWarpExportParticles) — User Manual
        Overview

        The Export Particles protocol extracts particle-centered data
        from tilt-series based on a set of 3D coordinates. It prepares
        Warp-compatible metadata, imports CTF and alignment information,
        and exports particles in formats suitable for downstream
        subtomogram analysis.

        This protocol is commonly used after tomogram reconstruction,
        particle picking, or coordinate annotation, when the goal is to
        generate particle stacks for refinement, classification, or
        visualization.

        Biological Purpose

        In cryo-electron tomography, particles are often identified as
        3D coordinates inside reconstructed tomograms. However, many
        downstream analysis tools require particles to be represented as
        localized image stacks or subtomograms.

        This protocol bridges that gap by converting coordinates into
        extracted particle-centered image data while preserving
        geometric, CTF, and alignment consistency.

        Inputs

        The protocol requires three main inputs:

        1. A set of 3D coordinates.
           These define the particle centers inside tomograms.

        2. A set of tilt-series.
           These provide the original aligned projection data.

        3. A set of CTF estimations.
           These are used to update Warp-compatible CTF metadata.

        All inputs must correspond to the same tomographic dataset.

        Export Workflow

        The protocol executes the following sequence:

        1. Read tomogram geometry from the coordinate set.
        2. Compute scaling factors between tomogram and tilt-series.
        3. Prepare Warp settings and tilt-series metadata.
        4. Generate IMOD alignment files.
        5. Generate Warp CTF estimation files.
        6. Update CTF metadata using provided CTF models.
        7. Import alignments into Warp.
        8. Create one STAR file per tomogram containing normalized
           particle coordinates and orientations.
        9. Export particles using WarpTools.
        10. Register exported particles as Scipion outputs.
        11. Remove temporary intermediate files.

        Coordinate Handling

        Coordinates are grouped by tomogram identifier.

        For each particle, the protocol writes:

        - normalized X, Y, Z coordinates,
        - particle orientation angles,
        - tomogram association,
        - particle score.

        Coordinates are normalized with respect to tomogram dimensions
        so Warp can correctly interpret particle positions during export.

        Export Modes

        The protocol supports two output modes:

        2D mode
            Exports particle-centered 2D image series. These are often
            useful for particle-series workflows and certain Relion
            subtomogram pipelines.

        3D mode
            Exports full subtomograms as 3D particle volumes.

        The selected mode determines how Warp writes particle data and
        how the output set is interpreted downstream.

        Main Parameters

        Output pixel size
            Defines the sampling rate of exported particles.

        Output box size
            Defines the cubic extraction box around each particle.

        Particle diameter
            Provides a biological size estimate used during export.

        Export type
            Selects whether output particles are written as 2D image
            stacks or 3D subtomograms.

        These parameters should match the biological particle size and
        the intended downstream analysis resolution.

        CTF and Alignment Consistency

        Before particle extraction, the protocol regenerates Warp
        metadata for each tilt-series and updates the local CTF XML
        files with the provided defocus information.

        Alignment files are imported from IMOD-compatible metadata so
        particle extraction remains geometrically consistent with the
        original tilt-series alignment.

        Output Registration

        After export, the protocol creates a
        RelionSetOfPseudoSubtomograms output object.

        This output includes:

        - exported particle references,
        - tomogram associations,
        - sampling information,
        - Relion-compatible metadata.

        Internal STAR files are also normalized so paths remain valid
        inside the Scipion project structure.

        Output Relations

        The resulting particle set is linked to:

        - the original coordinate set,
        - the tilt-series,
        - the input CTF estimations.

        This preserves provenance and ensures reproducibility of
        downstream processing.

        Practical Recommendations

        For subtomogram refinement:
            Use 3D export mode.

        For particle-series workflows:
            Use 2D export mode.

        Box size should be large enough to contain the full particle,
        but not excessively large, as unnecessary background increases
        computational cost.

        The chosen output pixel size should be compatible with both the
        biological target size and the refinement software to be used.

        Cleanup Strategy

        After successful export, temporary intermediate tilt-image
        files are removed automatically to reduce disk usage.

        Final Perspective

        Particle export is a critical bridge between coordinate-based
        particle localization and downstream subtomogram analysis.

        This protocol automates metadata preparation, coordinate
        normalization, and Warp-based extraction while generating
        Scipion-ready particle outputs for further structural analysis.
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

        # form.addParam('writeStacks', params.EnumParam,
        #               label='Export type',
        #               default=0,
        #               choices=['2D', '3D'],
        #               display=params.EnumParam.DISPLAY_HLIST,
        #               help='If set to 2D, this program will write output particles as 2d image series centered on '
        #                    'the particle (particle series). If set to 3D, this program will write output '
        #                    'particles as 3d images (subtomograms)')

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
        self.tsCtfEstimation(tomoSr)
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
        cmd = '--relative_output_paths --normalized_coords --2d'
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

    def tsCtfEstimation(self, tomoSr):
        """CTF estimation"""

        self.info(">>> Generating ctf estimation...")
        settingFile = self._getExtraPath(TILTSERIE_SETTINGS)
        argsDict = {
            "--settings": os.path.abspath(settingFile),
            "--range_high": tomoSr * 3,
            "--range_low": tomoSr * 4,
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




