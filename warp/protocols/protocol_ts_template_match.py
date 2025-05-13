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

import emtable
import numpy as np


from pwem.emlib.image.image_readers import ImageStack, ImageReadersRegistry
from pyworkflow import BETA
import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from pyworkflow.object import Integer, Set
import tomo.objects as tomoObj
from tomo.constants import BOTTOM_LEFT_CORNER
from tomo.protocols import ProtTomoPicking

from warp import Plugin, RLN_COORDINATE_X, RLN_COORDINATE_Y, RLN_COORDINATE_Z, RLN_ANGLE_ROT, RLN_ANGLE_TILT, \
    RLN_ANGLE_PSI, RLN_AUTOPICK_FIGURE_OF_MERTI
from warp.constants import (TILTSERIE_SETTINGS, TILTSERIES_FOLDER, TS_CTF,
                            OUTPUT_CTF_SERIE, TS_RECONSTRUCTION, MRC_EXT, OUTPUT_TOMOGRAMS_NAME,
                            RECONSTRUCTION_FOLDER, RECONSTRUCTION_ODD_FOLDER, RECONSTRUCTION_EVEN_FOLDER,
                            TS_TEMPLATE_MATCH, CREATE_SETTINGS, SETTINGS_FOLDER, TOMOSTAR_FOLDER, TILTIMAGES_FOLDER,
                            AVERAGE_FOLDER, TS_THRESHOLD_PICKS, MATCHING_FOLDER, TS_IMPORT_ALIGNMENTS)
from warp.protocols.protocol_base import ProtWarpBase
from warp.utils import tomoStarGenerate, genTransformMatrix, updateCtFXMLFile


class ProtWarpTSTemplateMatch(ProtWarpBase, ProtTomoPicking):
    """
    Match previously reconstructed tomograms against a 3D template, producing a list of the highest-scoring matches
    More info:
        https://warpem.github.io/warp/user_guide/warptools/quick_start_warptools_tilt_series/#particle-picking
        https://warpem.github.io/warp/user_guide/warptools/quick_start_warptools_tilt_series/#particle-picking
    """

    _label = 'tomo picking'
    _devStatus = BETA

    _possibleOutputs = {'output3DCoordinates': tomoObj.SetOfCoordinates3D}

    def __init__(self, **kwargs):
        ProtTomoPicking.__init__(self, **kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        ProtTomoPicking._defineParams(self, form)
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

        form.addSection(label="Template matching")

        form.addParam('templateVolume', params.PointerParam, pointerClass='Volume',
                      important=True,
                      label="Volume Template",
                      help='Map that will serve as a template')

        form.addParam('template_diameter', params.IntParam, default=None,
                      important=True,
                      label='Template diameter (Å)', help='Template diameter in Angstrom')

        form.addParam('template_flip', params.BooleanParam, default=True,
                      label='Mirror the template along the X axis?',
                      help="Mirror the template along the X axis to flip the handedness; '_flipx' will be added to the template's name")

        form.addParam('symmetry', params.StringParam, default='C1',
                      label='Symmetry',
                      help='Default: C1. Symmetry of the template, e.g. C1, D7, O')

        form.addParam('subdivisions', params.IntParam, default=3,
                      label='Number of subdivisions defining the angular search step',
                      help="Default: 3. Number of subdivisions defining the angular search step: 2 = 15° step, "
                           "3 = 7.5°, 4 = 3.75° and so on")

        form.addParam('tilt_range', params.IntParam, default=None,
                      allowsNull=True,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Tilt range',
                      help="Limit the range of angles between the reference's Z axis and the tomogram's XY plane to "
                           "plus/minus this value, in °; useful for matching filaments lying flat in the XY plane")

        form.addParam('batch_angles', params.IntParam, default=32,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Batch angles',
                      help="How many orientations to evaluate at once; memory consumption scales linearly with this; "
                           "higher than 32 probably won't lead to speed-ups")

        form.addParam('peak_distance', params.FloatParam, default=None,
                      allowsNull=True,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Peak distance (Å)', help='Minimum distance (in Angstrom) between peaks; '
                                                      'leave empty to use template diameter')

        form.addParam('npeaks', params.IntParam, default=2000,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Maximum number of peak',
                      help="Maximum number of peak positions to save")

        form.addParam('dont_normalize', params.BooleanParam, default=True,
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Don't normalize?",
                      help="Don't set score distribution to median = 0, stddev = 1")

        form.addParam('whiten', params.BooleanParam, default=True,
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Perform spectral whitening?",
                      help="Perform spectral whitening to give higher-resolution information more weight; "
                           "this can help when the alignments are already good and you need more selective matching")

        form.addParam('lowpass', params.FloatParam, default=1.0,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Gaussian low-pass filter',
                      help="Gaussian low-pass filter to be applied to template and tomogram, in "
                           "fractions of Nyquist; 1.0 = no low-pass, <1.0 = low-pass")

        form.addParam('lowpass_sigma', params.FloatParam, default=0.1,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Sigma (i.e. fall-off) of the Gaussian low-pass filter',
                      help="Sigma (i.e. fall-off) of the Gaussian low-pass filter, in fractions of Nyquist; larger "
                           "value = slower fall-off")

        form.addParam('max_missing_tilts', params.IntParam, default=2,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Maximum dismiss positions',
                      help="Dismiss positions not covered by at least this many tilts; set to -1 to disable position culling")

        form.addParam('reuse_results', params.BooleanParam, default=True,
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Reuse correlation volumes?",
                      help="Reuse correlation volumes from a previous run if available, only extract peak positions")

        form.addParam('check_hand', params.IntParam, default=0,
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Check hand",
                      help="Also try a flipped version of the template on this many tomograms to see what geometric hand they have")

        form.addParam('subvolume_size', params.IntParam, default=192,
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Subvolume size",
                      help="Matching is performed locally using sub-volumes of this size in pixele")

        form.addHidden(params.GPU_LIST, params.StringParam, default='0',
                       label='Choose GPU IDs:', validators=[params.NonEmpty],
                       help="Space-separated list of GPU IDs to use for processing. Default: all GPUs in the system."
                            " Warp can use multiple GPUs - in that case"
                            " set to i.e. *0 1 2*.")

        # form.addParam('apply_score', params.BooleanParam, default=False,
        #               label="Apply a score threshold to particles picked?",
        #               help="Apply a score threshold to particles picked through template-matching from tilt")
        #
        # form.addParam('minimum', params.IntParam, default=3,
        #               allowsNull=True,
        #               condition='apply_score',
        #               label="Minimum threshold",
        #               help="Remove all particles below this threshold")
        # form.addParam('maximum', params.IntParam, default=None,
        #               condition='apply_score',
        #               label="Maximun threshold",
        #               help="Remove all particles above this threshold")

    def _insertAllSteps(self):
        inputTomograms = self.inputTomograms.get()
        inputTs = self.inputSet.get()

        tsSr = inputTs.getSamplingRate()
        tomoSr = inputTomograms.getSamplingRate()
        tomo_dim = inputTomograms.getDim()

        scaleFactor = tomoSr / tsSr
        self.tomo_thickness = Integer(round(tomo_dim[2] * scaleFactor))
        self.x_dimension = Integer(round(tomo_dim[0] * scaleFactor))
        self.y_dimension = Integer(round(tomo_dim[1] * scaleFactor))

        self._insertFunctionStep(self.tomogramsPrepare, needsGPU=False)

        for tomogram in inputTomograms.iterItems(iterate=False):
            tsId = tomogram.getTsId()
            ts = inputTs.getItem('_tsId', tsId)
            self._insertFunctionStep(self.templateMatchStep, ts, needsGPU=True)
            self._insertFunctionStep(self.createOutputStep, ts, needsGPU=False)
            self._insertFunctionStep(self.cleanIntermediateResults, needsGPU=False)

            # if self.apply_score.get():
            #     self._insertFunctionStep(self.applyScoreStep, ts, needsGPU=True)

        self._insertFunctionStep(self._closeOutputSet, needsGPU=False)

    def tomogramsPrepare(self):
        inputTomograms = self.inputTomograms.get()
        processingFolder = os.path.abspath(self._getExtraPath(TILTSERIES_FOLDER))
        tomogramFolder = os.path.join(processingFolder, RECONSTRUCTION_FOLDER)
        pwutils.makePath(tomogramFolder)
        for tomogram in inputTomograms:
            fileName = tomogram.getFileName()
            destFolder = os.path.join(tomogramFolder, os.path.basename(fileName))
            os.symlink(fileName, destFolder)

    def tsCtfEstimation(self, ts):
        """CTF estimation"""

        self.info(">>> Generating ctf estimation file for %s..." % ts.getTsId())
        # inputTSAdquisition = ts.getAcquisition()
        settingFile = self._getExtraPath(SETTINGS_FOLDER, ts.getTsId() + '_' + TILTSERIE_SETTINGS)
        argsDict = {
            "--settings": os.path.abspath(settingFile),
            # "--window": 512,
            # "--range_low": 30,
            "--range_high": 7,
            # # "--range_high": self.inputSet.get().getSamplingRate() * 2 + 0.1,
            # "--defocus_min": 0.5,
            "--defocus_max": 8,
            # "--voltage": int(inputTSAdquisition.getVoltage()),
            # "--cs": inputTSAdquisition.getSphericalAberration(),
            # "--amplitude": inputTSAdquisition.getAmplitudeContrast(),
        }
        try:
            self.runProgram(argsDict, TS_CTF)
        except Exception:
            self.info(">>> Error generating ctf estimation file...")
        ctfTomoSeries = self.inputSetOfCtfTomoSeries.get().getItem('_tsId', ts.getTsId())
        processingFolder = os.path.abspath(self._getExtraPath(TILTSERIES_FOLDER))
        defocusFilePath = os.path.join(processingFolder, ts.getTsId() + '.xml')
        updateCtFXMLFile(defocusFilePath, ctfTomoSeries)

    def tsImportAligments(self, ts):
        processingFolder = os.path.abspath(self._getExtraPath(TILTSERIES_FOLDER))
        tiltstackFolder = os.path.join(processingFolder, 'tiltstack', ts.getTsId())
        pwutils.makePath(tiltstackFolder)
        ts.writeImodFiles(tiltstackFolder, delimiter=' ')
        self.info(">>> Starting import aligments...")
        angpix = ts.getSamplingRate()
        settingFile = self._getExtraPath(SETTINGS_FOLDER, ts.getTsId() + '_' + TILTSERIE_SETTINGS)
        argsDict = {
            "--settings": os.path.abspath(settingFile),
            '--alignments': os.path.abspath(tiltstackFolder),
            "--alignment_angpix": angpix,
        }
        cmd = ' '.join(['%s %s' % (k, v) for k, v in argsDict.items()])
        self.runJob(self.getPlugin().getProgram(TS_IMPORT_ALIGNMENTS), cmd, executable='/bin/bash')

    def templateMatchStep(self, ts):
        """Particle Picking"""
        self.createTiltSeriesSetting(ts)
        self.tsDataPrepare(ts)
        self.tsCtfEstimation(ts)
        self.tsImportAligments(ts)
        self.info(">>> Starting particle picking...")
        angpix = self.inputTomograms.get().getSamplingRate()
        settingFile = self._getExtraPath(SETTINGS_FOLDER, ts.getTsId() + '_' + TILTSERIE_SETTINGS)
        argsDict = {
            "--settings": os.path.abspath(settingFile),
            "--tomo_angpix": angpix,
            "--subdivisions": self.subdivisions.get(),
            "--template_path": self.templateVolume.get().getFileName(),
            "--template_diameter": self.template_diameter.get(),
            "--symmetry": self.symmetry.get(),
            "--check_hand": self.check_hand.get()
        }

        cmd = ''
        if self.whiten.get():
            cmd += " --whiten"

        self.runProgram(argsDict, TS_TEMPLATE_MATCH, othersCmds=cmd)

    def applyScoreStep(self, ts):
        """Apply a score threshold to particles picked through template-matching from tilt"""
        tsId = ts.getTsId()
        self.info(">>> Starting to apply a score threshold to particles picked to %s..." % tsId)
        settingFile = self._getExtraPath(SETTINGS_FOLDER, tsId + '_' + TILTSERIE_SETTINGS)
        argsDict = {
            "--settings": os.path.abspath(settingFile),
            "--in_suffix": tsId,
            "--out_suffix": 'clean',
            "--minimum": self.minimum.get()
        }
        self.runProgram(argsDict, TS_THRESHOLD_PICKS)

    def getOutputSetOfCoordinates3D(self, outputSetName):
        suffix = self._getOutputSuffix(tomoObj.SetOfCoordinates3D)
        inputTomos = self.inputTomograms.get()
        outputSetOfCoordinates3D = getattr(self, outputSetName, None)
        if outputSetOfCoordinates3D:
            outputSetOfCoordinates3D.enableAppend()
        else:
            outputSetOfCoordinates3D = self._createSetOfCoordinates3D(self.inputTomograms, suffix)
            outputSetOfCoordinates3D.setName("tomoCoord")
            outputSetOfCoordinates3D.setSamplingRate(inputTomos.getSamplingRate())
            outputSetOfCoordinates3D.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(**{outputSetName: outputSetOfCoordinates3D})
            self._defineSourceRelation(self.inputTomograms, outputSetOfCoordinates3D)

        return outputSetOfCoordinates3D

    def createOutputStep(self, ts):
        tsId = ts.getTsId()
        setOfTomograms = self.inputTomograms.get()
        tomogram = setOfTomograms.getItem('_tsId', tsId)
        tomoDim = tomogram.getDim()
        sr = tomogram.getSamplingRate()
        tomoFileName = os.path.basename(tomogram.getFileName())
        tomoFileBaseName = os.path.splitext(tomoFileName)[0]
        outputPath = self._getExtraPath(TILTSERIES_FOLDER, MATCHING_FOLDER)
        coordsFile = [f for f in os.listdir(outputPath) if f.startswith(tomoFileBaseName) and f.endswith(".star") and "_flipx" not in f][0]

        setOfCoord3D = self.getOutputSetOfCoordinates3D('output3DCoordinates')
        origin = BOTTOM_LEFT_CORNER

        coord = tomoObj.Coordinate3D()
        coordsFile = os.path.join(outputPath, coordsFile)
        table = emtable.Table(fileName=coordsFile, guessType=False)

        for row in table.iterRows(coordsFile,  guessType=False):
            # Clean up objId to add as a new coordinate
            coord.setObjId(None)
            coord.setVolume(tomogram.clone())
            x = float(row.get(RLN_COORDINATE_X)) * tomoDim[0]
            y = float(row.get(RLN_COORDINATE_Y)) * tomoDim[1]
            z = float(row.get(RLN_COORDINATE_Z)) * tomoDim[2]

            coord.setPosition(x, y, z, origin)
            M = genTransformMatrix(float(row.get(RLN_ANGLE_ROT)),
                                   float(row.get(RLN_ANGLE_TILT)),
                                   float(row.get(RLN_ANGLE_PSI)))
            coord.setMatrix(M)
            coord.setTomoId(tsId)
            coord.setScore(float((row.get(RLN_AUTOPICK_FIGURE_OF_MERTI))))
            setOfCoord3D.append(coord)
        setOfCoord3D.write()
        setOfCoord3D.setBoxSize(int(self.template_diameter.get()/sr))
        setOfCoord3D.initTomos()
        self._store(setOfCoord3D)

    def allowsDelete(self, obj):
        return True

    def _validate(self):
        errors = []
        # if self.apply_score.get() and self.maximum.get() is None and self.minimum.get() is None:
        #     errors.append('The minimum or maximum threshold value must be greater than 0')
        return errors

    def cleanIntermediateResults(self):
        self.info(">>> Cleaning intermediate results...")
        imagesFolder = self._getExtraPath(TILTIMAGES_FOLDER)
        pwutils.cleanPath(imagesFolder)


