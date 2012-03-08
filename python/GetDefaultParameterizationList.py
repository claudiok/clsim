import string
from os.path import expandvars, exists, isdir, isfile

from icecube import icetray, dataclasses
from icecube.clsim import I3CLSimLightSourceParameterization

def GetDefaultParameterizationList(theConverter, muonOnly=False):
    fromEnergy=0.
    toEnergy=float('Inf')

    muons    = [dataclasses.I3Particle.MuMinus,
                dataclasses.I3Particle.MuPlus]

    cascades = [dataclasses.I3Particle.Neutron,
                dataclasses.I3Particle.Hadrons,
                dataclasses.I3Particle.Pi0,
                dataclasses.I3Particle.PiPlus,
                dataclasses.I3Particle.PiMinus,
                dataclasses.I3Particle.K0_Long,
                dataclasses.I3Particle.KPlus,
                dataclasses.I3Particle.KMinus,
                dataclasses.I3Particle.PPlus,
                dataclasses.I3Particle.PMinus,
                dataclasses.I3Particle.K0_Short,
                dataclasses.I3Particle.EMinus,
                dataclasses.I3Particle.EPlus,
                dataclasses.I3Particle.Gamma,
                dataclasses.I3Particle.Brems,
                dataclasses.I3Particle.DeltaE,
                dataclasses.I3Particle.PairProd,
                dataclasses.I3Particle.NuclInt]

    parameterizationsMuon = []
    for type in muons:
        converter = \
          I3CLSimLightSourceParameterization(
            converter=theConverter,
            forParticleType=type,
            fromEnergy=fromEnergy,
            toEnergy=toEnergy, 
            needsLength=True)
        parameterizationsMuon.append(converter)

    if muonOnly:
        return parameterizationsMuon

    parameterizationsOther = []
    for type in cascades:
        converter = \
          I3CLSimLightSourceParameterization(
            converter=theConverter,
            forParticleType=type,
            fromEnergy=fromEnergy,
            toEnergy=toEnergy, 
            needsLength=False)
        parameterizationsOther.append(converter)

    return parameterizationsMuon+parameterizationsOther
