#
# Copyright (c) 2011, 2012
# Claudio Kopper <claudio.kopper@icecube.wisc.edu>
# and the IceCube Collaboration <http://www.icecube.wisc.edu>
# 
# Permission to use, copy, modify, and/or distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
# 
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
# SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
# OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
# CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
# 
# 
# $Id: GetIceCubeFlasherSpectrum.py 108199 2013-07-12 21:33:08Z nwhitehorn $
# 
# @file GetIceCubeFlasherSpectrum.py
# @version $Revision: 108199 $
# @date $Date: 2013-07-12 17:33:08 -0400 (Fri, 12 Jul 2013) $
# @author Claudio Kopper
#

from icecube import icetray, dataclasses

from icecube.clsim import I3CLSimFunctionFromTable
from icecube.clsim import I3CLSimFunctionDeltaPeak
from icecube.clsim import I3CLSimFlasherPulse

from I3Tray import I3Units

import numpy, math
from os.path import expandvars

def GetIceCubeFlasherSpectrumData(spectrumType):
    if spectrumType == I3CLSimFlasherPulse.FlasherPulseType.LED340nm:
        data = numpy.loadtxt(expandvars("$I3_BUILD/clsim/resources/flasher_data/flasher_led_340nm_emission_spectrum_cw_measured_20mA_pulseCurrent.txt"), unpack=True)
        #data = numpy.loadtxt(expandvars("$I3_BUILD/clsim/resources/flasher_data/flasher_led_340nm_emission_spectrum_cw_measured_200mA_pulseCurrent.txt"), unpack=True)
        data[0] *= I3Units.nanometer # apply the correct units
        data[1] /= 24.306508         # (20mA pulse)  pre-calculated normalization constant (not really necessary for the generation algorithm)
        #data[1] /= 22.323254        # (200ma pulse) pre-calculated normalization constant (not really necessary for the generation algorithm)
    elif spectrumType == I3CLSimFlasherPulse.FlasherPulseType.LED370nm:
        data = numpy.loadtxt(expandvars("$I3_BUILD/clsim/resources/flasher_data/flasher_led_370nm_emission_spectrum_cw_measured.txt"), unpack=True)
        data[0] *= I3Units.nanometer # apply the correct units
        data[1] /= 15.7001863        # pre-calculated normalization constant (not really necessary for the generation algorithm)
    elif spectrumType == I3CLSimFlasherPulse.FlasherPulseType.LED405nm:
        data = numpy.loadtxt(expandvars("$I3_BUILD/clsim/resources/flasher_data/flasher_led_405nm_emission_spectrum_datasheet.txt"), unpack=True)
        data[0] *= I3Units.nanometer # apply the correct units
        data[1] /= 8541585.10324     # pre-calculated normalization constant (not really necessary for the generation algorithm)
    elif spectrumType == I3CLSimFlasherPulse.FlasherPulseType.LED450nm:
        data = numpy.loadtxt(expandvars("$I3_BUILD/clsim/resources/flasher_data/flasher_led_450nm_emission_spectrum_datasheet.txt"), unpack=True)
        data[0] *= I3Units.nanometer # apply the correct units
        data[1] /= 21.9792812618     # pre-calculated normalization constant (not really necessary for the generation algorithm)
    elif spectrumType == I3CLSimFlasherPulse.FlasherPulseType.LED505nm:
        data = numpy.loadtxt(expandvars("$I3_BUILD/clsim/resources/flasher_data/flasher_led_505nm_emission_spectrum_cw_measured.txt"), unpack=True)
        data[0] *= I3Units.nanometer # apply the correct units
        data[1] /= 38.1881           # pre-calculated normalization constant (not really necessary for the generation algorithm)
    else:
        raise RuntimeError("invalid spectrumType")
        
    ## auto-normalize (not active in order not to depend on scipy)
    #import scipy, scipy.integrate
    #integral = scipy.integrate.trapz(y=data[1], x=data[0]/I3Units.nanometer)
    #data[1] /= integral
    #print integral
    
    return data

def GetIceCubeFlasherSpectrum(spectrumType = I3CLSimFlasherPulse.FlasherPulseType.LED405nm):
    if spectrumType in [I3CLSimFlasherPulse.FlasherPulseType.SC1, I3CLSimFlasherPulse.FlasherPulseType.SC2]:
        # special handling for Standard Candles:
        # these currently have single-wavelength spectra
        spectrum = I3CLSimFunctionDeltaPeak(337.*I3Units.nanometer)
        return spectrum
    else:
        data = GetIceCubeFlasherSpectrumData(spectrumType)
        spectrum = I3CLSimFunctionFromTable(data[0], data[1])
        return spectrum
