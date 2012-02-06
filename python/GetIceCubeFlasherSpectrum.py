from icecube import icetray, dataclasses

from . import I3CLSimWlenDependentValueFromTable
from . import I3CLSimFlasherPulse

from I3Tray import I3Units

import numpy, math
from os.path import expandvars

def GetIceCubeFlasherSpectrumData(spectrumType):
    if spectrumType == I3CLSimFlasherPulse.FlasherPulseType.LED340nm:
        data = numpy.loadtxt(expandvars("$I3_SRC/clsim/resources/flasher_data/flasher_led_340nm_emission_spectrum_cw_measured_20mA_pulseCurrent.txt"), unpack=True)
        #data = numpy.loadtxt(expandvars("$I3_SRC/clsim/resources/flasher_data/flasher_led_340nm_emission_spectrum_cw_measured_200mA_pulseCurrent.txt"), unpack=True)
        data[0] *= I3Units.nanometer # apply the correct units
        data[1] /= 24.306508         # (20mA pulse)  pre-calculated normalization constant (not really necessary for the generation algorithm)
        #data[1] /= 22.323254        # (200ma pulse) pre-calculated normalization constant (not really necessary for the generation algorithm)
    elif spectrumType == I3CLSimFlasherPulse.FlasherPulseType.LED370nm:
        data = numpy.loadtxt(expandvars("$I3_SRC/clsim/resources/flasher_data/flasher_led_370nm_emission_spectrum_cw_measured.txt"), unpack=True)
        data[0] *= I3Units.nanometer # apply the correct units
        data[1] /= 15.7001863        # pre-calculated normalization constant (not really necessary for the generation algorithm)
    elif spectrumType == I3CLSimFlasherPulse.FlasherPulseType.LED405nm:
        data = numpy.loadtxt(expandvars("$I3_SRC/clsim/resources/flasher_data/flasher_led_405nm_emission_spectrum_datasheet.txt"), unpack=True)
        data[0] *= I3Units.nanometer # apply the correct units
        data[1] /= 8541585.10324     # pre-calculated normalization constant (not really necessary for the generation algorithm)
    elif spectrumType == I3CLSimFlasherPulse.FlasherPulseType.LED450nm:
        data = numpy.loadtxt(expandvars("$I3_SRC/clsim/resources/flasher_data/flasher_led_450nm_emission_spectrum_datasheet.txt"), unpack=True)
        data[0] *= I3Units.nanometer # apply the correct units
        data[1] /= 21.9792812618     # pre-calculated normalization constant (not really necessary for the generation algorithm)
    elif spectrumType == I3CLSimFlasherPulse.FlasherPulseType.LED505nm:
        data = numpy.loadtxt(expandvars("$I3_SRC/clsim/resources/flasher_data/flasher_led_505nm_emission_spectrum_cw_measured.txt"), unpack=True)
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
    data = GetIceCubeFlasherSpectrumData(spectrumType)
    spectrum = I3CLSimWlenDependentValueFromTable(data[0], data[1])
    return spectrum
