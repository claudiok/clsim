from icecube import icetray, dataclasses
from icecube.clsim import I3CLSimLightSourceParameterization
from icecube.clsim import I3CLSimFlasherPulse, I3CLSimLightSourceToStepConverterFlasher
from icecube.clsim import GetIceCubeFlasherSpectrum

def GetFlasherParameterizationList(spectrumTable):
    spectrumTypes = [I3CLSimFlasherPulse.FlasherPulseType.LED340nm,
                     I3CLSimFlasherPulse.FlasherPulseType.LED370nm,
                     I3CLSimFlasherPulse.FlasherPulseType.LED405nm,
                     I3CLSimFlasherPulse.FlasherPulseType.LED450nm,
                     I3CLSimFlasherPulse.FlasherPulseType.LED505nm]
    
    parameterizations = []
    for flasherSpectrumType in spectrumTypes:
        theSpectrum = GetIceCubeFlasherSpectrum(spectrumType=flasherSpectrumType)
        theConverter = I3CLSimLightSourceToStepConverterFlasher(flasherSpectrumNoBias=theSpectrum, spectrumTable=spectrumTable)
        parameterization = I3CLSimLightSourceParameterization(converter=theConverter, forFlasherPulseType=flasherSpectrumType)
        parameterizations.append(parameterization)
        
    return parameterizations
