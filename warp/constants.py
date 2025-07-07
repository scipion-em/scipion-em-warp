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


def getWarpEnvName(version):
    return "warp-%s" % version


V2_0_0 = "2.0.0"

VERSIONS = [V2_0_0]
WARP_DEFAULT_VER_NUM = V2_0_0

DEFAULT_ENV_NAME = getWarpEnvName(WARP_DEFAULT_VER_NUM)
DEFAULT_ACTIVATION_CMD = 'conda activate ' + DEFAULT_ENV_NAME
WARP_ENV_ACTIVATION = 'WARP_ENV_ACTIVATION'
WARP_LOADER = 'WARP_LOADER'
WARP_FORCE_MRC_FLOAT32 = 'WARP_FORCE_MRC_FLOAT32'


OUTPUT_CTF_SERIE = "CTFTomoSeries"
OUPUT_AVERAGE_SUBTOMOGRAM = "AverageSubTomogram"
OUTPUT_RELION_PARTICLES = "RelionParticles"
OUTPUT_TOMOGRAMS_NAME = "Tomograms"
OUTPUT_TILTSERIES = 'TiltSeries'
OUTPUT_HANDEDNESS = 'Handedness_OK'
EXT_MRC_EVEN_NAME = "even.mrc"
EXT_MRC_ODD_NAME = "odd.mrc"
EVEN = 'even'
ODD = 'odd'
MRC_EXT = 'mrc'


# --------- [WARPPTOOLS PROGRAMS] ---------------
WARP_TOOLS = 'WarpTools'
CREATE_SETTINGS = 'create_settings'
FS_MOTION = 'fs_motion'
TS_CTF = 'ts_ctf'
FS_CTF = 'fs_ctf'
FS_MOTION = 'fs_motion'
FS_MOTION_AND_CTF = 'fs_motion_and_ctf'
TS_RECONSTRUCTION = 'ts_reconstruct'
TS_DEFOCUS_HAND = 'ts_defocus_hand'
TS_TEMPLATE_MATCH = 'ts_template_match'
TS_THRESHOLD_PICKS = 'threshold_picks'
TS_IMPORT_ALIGNMENTS = 'ts_import_alignments'
TS_EXPORT_PARTICLES = 'ts_export_particles'

WARP_TOOLS_GPU_ALGORITHMS = [FS_MOTION, TS_CTF, FS_MOTION, FS_MOTION_AND_CTF, TS_RECONSTRUCTION, TS_TEMPLATE_MATCH,
                             TS_THRESHOLD_PICKS, TS_EXPORT_PARTICLES]

# ---------[MTOOLS PROGRAMS]-------------------
MTOOLS = 'MTools'
CREATE_POPULATION = 'create_population'
CREATE_SOURCES = 'create_source'
CREATE_SPECIES = 'create_species'
MCORE = 'MCore'


# -------------[LABELS]----------------
TOMOSTAR_FOLDER = 'tomostar'
TILTIMAGES_FOLDER = 'tiltimages'
FRAMES_FOLDER = 'frames'
AVERAGE_FOLDER = 'average'
SETTINGS_FOLDER = 'settings'
TILTSERIES_FOLDER = 'warp_tiltseries'
TILTSERIE_SETTINGS = "warp_tiltseries.settings"
POWERSPECTRUM_FOLDER = 'powerspectrum'
FRAMESERIES_FOLDER = 'warp_frameseries'
FRAMESERIES_SETTINGS = "warp_frameseries.settings"
RECONSTRUCTION_FOLDER = 'reconstruction'
RECONSTRUCTION_ODD_FOLDER = 'odd'
RECONSTRUCTION_EVEN_FOLDER = 'even'
MATCHING_FOLDER = 'matching'
RELION_FOLDER = 'relion'
M_RESULT_FOLDER = 'm'
SPECIES_FOLDER = 'species'
OPTIMISATION_SET_STAR = 'matching_optimisation_set.star'
MATCHING_PARTICLES_STAR = 'matching.star'
IN_PARTICLES_STAR = 'inParticles.star'
MATCHING_TOMOGRAMS_STAR = 'matching_tomograms.star'
MATCHING_PROCESSING_SPECIES_PATTERN = 'processing_species_*'
PROCESSING_SPECIES_AVERAGE = 'processing_species_filt.mrc'
PROCESSING_SPECIES_HALF1 = 'processing_species_half1.mrc'
PROCESSING_SPECIES_HALF2 = 'processing_species_half2.mrc'
PROCESSING_SPECIES_PARTICLES = 'processing_species_particles.star'

# RELION LABELS
RLN_COORDINATE_X = 'rlnCoordinateX'
RLN_COORDINATE_Y = 'rlnCoordinateY'
RLN_COORDINATE_Z = 'rlnCoordinateZ'
RLN_ANGLE_ROT = 'rlnAngleRot'
RLN_ANGLE_TILT = 'rlnAngleTilt'
RLN_ANGLE_PSI = 'rlnAnglePsi'
RLN_AUTOPICK_FIGURE_OF_MERIT = 'rlnAutopickFigureOfMerit'
RLN_MICROGRAPH_NAME = 'rlnMicrographName'

# WARP_LABELS
WRP_COORDINATE_X = 'wrpCoordinateX1'
WRP_COORDINATE_Y = 'wrpCoordinateY1'
WRP_COORDINATE_Z = 'wrpCoordinateZ1'
WRP_ANGLE_ROT = 'wrpAngleRot1'
WRP_ANGLE_TILT = 'wrpAngleTilt1'
WRP_ANGLE_PSI = 'wrpAnglePsi1'
WRP_RAMDOM_SUBSET = 'wrpRandomSubset'
WRP_SOURCE_NAME = 'wrpSourceName'
WRP_SOURCE_HASH = 'wrpSourceHash'

tomoStarFields = [
    RLN_COORDINATE_X,
    RLN_COORDINATE_Y,
    RLN_COORDINATE_Z,
    RLN_ANGLE_ROT,
    RLN_ANGLE_TILT,
    RLN_ANGLE_PSI,
    RLN_MICROGRAPH_NAME,
    RLN_AUTOPICK_FIGURE_OF_MERIT,
]