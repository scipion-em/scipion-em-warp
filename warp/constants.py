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
OUTPUT_TOMOGRAMS_NAME = "Tomograms"
OUTPUT_TILTSERIES = 'TiltSeries'
OUTPUT_HANDEDNESS = 'Warp handedness'
EXT_MRC_EVEN_NAME = "even.mrc"
EXT_MRC_ODD_NAME = "odd.mrc"
EVEN = 'even'
ODD = 'odd'
MRC_EXT = 'mrc'


# --------- [PROGRAMS] ---------------
CREATE_SETTINGS = 'create_settings'
FS_MOTION = 'fs_motion'
TS_CTF = 'ts_ctf'
FS_CTF = 'fs_ctf'
FS_MOTION = 'fs_motion'
FS_MOTION_AND_CTF = 'fs_motion_and_ctf'
TS_RECONSTRUCTION = 'ts_reconstruct'
TS_DEFOCUS_HAND = 'ts_defocus_hand'

# -------------[LABELS]----------------
TOMOSTAR_FOLDER = 'tomostar'
TILTIMAGES_FOLDER = 'tiltimages'
FRAMES_FOLDER = 'frames'
AVERAGE_FOLDER = 'average'
TILTSERIES_FOLDER = 'warp_tiltseries'
TILTSERIE_SETTINGS = "warp_tiltseries.settings"
FRAMESERIES_FOLDER = 'warp_frameseries'
FRAMESERIES_SETTINGS = "warp_frameseries.settings"
RECONSTRUCTION_FOLDER = 'reconstruction'
RECONSTRUCTION_ODD_FOLDER = 'odd'
RECONSTRUCTION_EVEN_FOLDER = 'even'
