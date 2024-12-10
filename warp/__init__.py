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

import os

import pwem
import pyworkflow.utils as pwutils

from warp.constants import *


__version__ = '3.4.3'
_references = ['Nickell2005', 'Tegunov2019']
_logo = "warp_logo.png"


class Plugin(pwem.Plugin):
    _url = "https://github.com/scipion-em/scipion-em-warp"

    @classmethod
    def _defineVariables(cls):
        cls._defineVar(WARP_ENV_ACTIVATION, DEFAULT_ACTIVATION_CMD, description="Command to activate warp environment.")
        cls._defineVar(WARP_LOADER, None, description="Command to load dependencies needed by the warp environment (e.g: moduleload ...)")
        cls._defineVar(WARP_FORCE_MRC_FLOAT32, "1", description="Force writing 32 bit MRCs. A 1 will produce bigger files but compatible with the rest or the methods.")

    @classmethod
    def getEnviron(cls):
        """ Setup the environment variables needed to launch Warp. """
        environ = pwutils.Environ(os.environ)
        if 'PYTHONPATH' in environ:
            # this is required for python virtual env to work
            del environ['PYTHONPATH']

        environ[WARP_FORCE_MRC_FLOAT32] = cls.getVar(WARP_FORCE_MRC_FLOAT32)
        return environ

    @classmethod
    def getDependencies(cls):
        """ Return a list of dependencies. Include conda if
        activation command was not found. """
        neededProgs = []
        if cls.getVar(WARP_LOADER) is not None:
            return neededProgs
        else:
            condaActivationCmd = cls.getCondaActivationCmd()
            if not condaActivationCmd:
                neededProgs.append('conda')
            return neededProgs

    @classmethod
    def defineBinaries(cls, env):
        for ver in VERSIONS:
            cls.addWarpPackage(env, ver,
                               default=ver == WARP_DEFAULT_VER_NUM)

    @classmethod
    def addWarpPackage(cls, env, version, default=False):
        ENV_NAME = getWarpEnvName(version)
        FLAG = f"warp_{version}_installed"

        # try to get CONDA activation command
        installCmds = [
            cls.getCondaActivationCmd(),
            f'conda create -y -n {ENV_NAME} warp=2.0.0dev31 -c warpem -c nvidia/label/cuda-11.8.0 -c pytorch -c conda-forge &&',
            f'touch {FLAG}'  # Flag installation finished
        ]

        envPath = os.environ.get('PATH', "")
        # keep path since conda likely in there
        installEnvVars = {'PATH': envPath} if envPath else None
        warpCmds = [(" ".join(installCmds), FLAG)]

        env.addPackage('warp', version=version,
                       tar='void.tgz',
                       commands=warpCmds,
                       neededProgs=cls.getDependencies(),
                       default=default,
                       vars=installEnvVars)

    @classmethod
    def getActivationCmd(cls):
        """ Return the activation command. """
        if cls.getVar(WARP_LOADER) is None:
            return f"{cls.getCondaActivationCmd()} {cls.getVar(WARP_ENV_ACTIVATION)}"
        else:
            return cls.getVar(WARP_LOADER)

    @classmethod
    def getProgram(cls, program):
        """ Create Warp command line. """
        return f"{cls.getActivationCmd()} && WarpTools {program}"
