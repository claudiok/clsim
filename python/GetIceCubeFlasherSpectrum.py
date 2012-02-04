from icecube import icetray, dataclasses
from icecube.clsim import I3CLSimRandomValueInterpolatedDistribution

from I3Tray import I3Units

import numpy, math
from os.path import expandvars

def GetIceCubeFlasherSpectrumData(colorIndex = 0):
    if colorIndex == 0:
        data = numpy.loadtxt(expandvars("$I3_SRC/clsim/resources/flasher_data/flasher_led_405nm_emission_spectrum_datasheet.txt"), unpack=True)
        data[0] *= I3Units.nanometer # apply the correct units
        data[1] /= 8541585.10324     # pre-calculated normalization constant (not really necessary for the generation algorithm)
    else:
        raise RuntimeError("invalid colorIndex")
    return data

def GetIceCubeFlasherSpectrumGenerator(colorIndex = 0):
    data = GetIceCubeFlasherSpectrumData(colorIndex)
    spectrumGen = I3CLSimRandomValueInterpolatedDistribution(data[0], data[1])
    return spectrumGen
