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

import numpy as np
import math
import scipy.fft


""" This code is adapted from https://github.com/dtegunov/tom_deconv
    and https://github.com/Heng-Z/IsoNet

    Original code is from:

    Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
    Journal of Structural Biology, 149 (2005), 227-234.

    Copyright (c) 2004-2007
    TOM toolbox for Electron Tomography
    Max-Planck-Institute of Biochemistry
    Dept. Molecular Structural Biology
    82152 Martinsried, Germany
    http://www.biochem.mpg.de/tom
"""


def tom_ctf1d(length, pixelsize, voltage, cs, defocus,
              amplitude, phaseshift, bfactor):

    ny = 1 / pixelsize
    lambda1 = 12.2643247 / np.sqrt(voltage * (1.0 + voltage * 0.978466e-6)) * 1e-10
    lambda2 = lambda1 * 2

    points = np.arange(0, length)
    points = points/(2 * length) * ny
    k2 = points ** 2
    term1 = lambda1**3 * cs * k2**2

    w = np.pi/2 * (term1 + lambda2 * defocus * k2) - phaseshift

    acurve = np.cos(w) * amplitude
    pcurve = -np.sqrt(1 - amplitude**2) * np.sin(w)
    bfactor = np.exp(-bfactor * k2 * 0.25)
    ctf = (pcurve + acurve) * bfactor

    return ctf


def tom_deconv(vol, angpix, voltage, cs, defocus, snrfalloff=1.1, deconvstrength=1,
               highpassnyquist=0.02, phaseflipped=False, phaseshift=0, ncpu=1):
    """
    :param vol: tomogram volume (or 2D image)
    :param angpix: angstrom per pixel
    :param voltage: voltage in kilovolts
    :param cs: cs in millimeters
    :param defocus: defocus in micrometers, positive = underfocus
    :param snrfalloff: how fast does SNR fall off, i.e. higher values will downweight high frequencies; values like 1.0 or 1.2 seem reasonable
    :param deconvstrength: how much will the signal be deconvoluted overall, i.e. a global scale for SNR; exponential scale: 1.0 is SNR = 1000 at zero frequency, 0.67 is SNR = 100, and so on
    :param highpassnyquist: fraction of Nyquist frequency to be cut off on the lower end (since it will be boosted the most)
    :param phaseflipped: whether the data are already phase-flipped
    :param phaseshift: CTF phase shift in degrees (e.g. from a phase plate)
    :param ncpu: number of CPUs for FFT

    Example:
    deconv = tom_deconv(mytomo, 3.42, 300, 2.7, 6, 1.1, 1, 0.02, False, 0);

    """

    data = np.arange(0, 1+1/2047, 1/2047)
    highpass = np.minimum(1, data / highpassnyquist) * np.pi
    highpass = 1 - np.cos(highpass)

    snr = np.exp(-data * snrfalloff * 100 / angpix) * (10**(3 * deconvstrength)) * highpass
    ctf = tom_ctf1d(2048, angpix * 1e-10, voltage * 1000, cs * 1e-3,
                    -defocus * 1e-6, 0.07, phaseshift / 180 * np.pi, 0)
    if phaseflipped:
        ctf = np.abs(ctf)

    wiener = ctf / (ctf * ctf + np.divide(1, np.maximum(1e-15, snr)))

    s1 = -math.floor(vol.shape[0] / 2)
    f1 = s1 + vol.shape[0] - 1
    m1 = np.arange(s1, f1 + 1)
    s2 = -math.floor(vol.shape[1] / 2)
    f2 = s2 + vol.shape[1] - 1
    m2 = np.arange(s2, f2 + 1)

    if vol.ndim == 3:
        s3 = -math.floor(vol.shape[2] / 2)
        f3 = s3 + vol.shape[2] - 1
        m3 = np.arange(s3, f3 + 1)
        x, y, z = np.meshgrid(m1, m2, m3, indexing='ij')
        x = np.divide(x, abs(s1))
        y = np.divide(y, abs(s2))
        z = np.divide(z, max(1, abs(s3)))
        r = np.sqrt(x**2 + y**2 + z**2)
        del x, y, z
    else:
        x, y = np.meshgrid(m1, m2, indexing='ij')
        x = np.divide(x, abs(s1))
        y = np.divide(y, abs(s2))
        r = np.sqrt(x ** 2 + y ** 2)
        del x, y

    r = np.minimum(1, r)
    r = np.fft.ifftshift(r)
    ramp = np.interp(r, data, wiener).astype(np.float32)
    vol = vol.astype(np.float32)

    deconv = np.real(scipy.fft.ifftn(scipy.fft.fftn(vol, overwrite_x=True, workers=ncpu) * ramp,
                                     overwrite_x=True, workers=ncpu))

    return deconv
