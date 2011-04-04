#!/usr/bin/env python

from optparse import OptionParser

usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)
parser.add_option("-o", "--outfile",
                  dest="OUTFILE", help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_option("-i", "--infile",
                  dest="INFILE", help="Read input from INFILE (.i3{.gz} format)")
parser.add_option("-s", "--seed",type="int",default=12345,
                                  dest="SEED", help="Initial seed for the random number generator")
#parser.add_option("-i", "--indexspectrum",type="float",default=2.,
#                              dest="INDEX", help="The spectral INDEX used for the generated neutrinos. The simulated spectrum will be E^-INDEX.")
parser.add_option("-r", "--runnumber", type="int", default=1,
                  dest="RUNNUMBER", help="The run number for this simulation")
#parser.add_option("-n", "--numevents", type="int", default=100000,
#                                  dest="NUMEVENTS", help="The number of events per run")

# parse cmd line args, bail out if anything is not understood
(options,args) = parser.parse_args()
if len(args) != 0:
        crap = "Got undefined options:"
        for a in args:
                crap += a
                crap += " "
        parser.error(crap)
#
if not options.OUTFILE:
    print "No OUTFILE specified!"
    parser.print_help()
    exit(-1)

if not options.INFILE:
    print "No INFILE specified!"
    parser.print_help()
    exit(-1)


from I3Tray import *
from os.path import expandvars
import os
import sys

from icecube import icetray, dataclasses, dataio, phys_services
from icecube import clsim
#load("libsim-services")

tray = I3Tray()

# a random number generator
randomService = phys_services.I3SPRNGRandomService(
    seed = options.SEED,
    nstreams = 10000,
    streamnum = options.RUNNUMBER)

# ice properties (SPICE-Mie model)
mediumProperties = clsim.MakeIceCubeMediumProperties()

# parameterizations for fast simulation (bypassing Geant4)
# converters first:
cascadeConverter = clsim.I3CLSimParticleToStepConverterCascadeParameterization(randomService=randomService)
cascadeConverter.SetMediumProperties(mediumProperties)
cascadeConverter.Initialize()


# now set up a list of converters with particle types and valid energy ranges
parameterizations = \
[
 #clsim.I3CLSimParticleParameterization(converter=cascadeConverter,
 #                                      forParticleType=dataclasses.I3Particle.MuMinus,
 #                                      fromEnergy=0.5*I3Units.GeV,
 #                                      toEnergy=1000.*I3Units.GeV),
 #clsim.I3CLSimParticleParameterization(converter=cascadeConverter,
 #                                      forParticleType=dataclasses.I3Particle.MuPlus,
 #                                      fromEnergy=0.5*I3Units.GeV,
 #                                      toEnergy=1000.*I3Units.GeV),

 # do we need some special handling for neutrons?
 clsim.I3CLSimParticleParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.Neutron,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.GeV),

 clsim.I3CLSimParticleParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.Hadrons,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.GeV),
 clsim.I3CLSimParticleParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.Pi0,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.GeV),
 clsim.I3CLSimParticleParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.PiPlus,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.GeV),
 clsim.I3CLSimParticleParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.PiMinus,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.GeV),
 clsim.I3CLSimParticleParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.K0_Long,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.GeV),
 clsim.I3CLSimParticleParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.KPlus,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.GeV),
 clsim.I3CLSimParticleParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.KMinus,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.GeV),
 clsim.I3CLSimParticleParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.PPlus,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.GeV),
 clsim.I3CLSimParticleParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.PMinus,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.GeV),
 clsim.I3CLSimParticleParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.K0_Short,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.GeV),
 clsim.I3CLSimParticleParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.NuclInt,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.GeV),

 clsim.I3CLSimParticleParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.EMinus,
                                       fromEnergy=0.5*I3Units.GeV,
                                       toEnergy=1000.*I3Units.GeV),
 clsim.I3CLSimParticleParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.EPlus,
                                       fromEnergy=0.5*I3Units.GeV,
                                       toEnergy=1000.*I3Units.GeV),
 clsim.I3CLSimParticleParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.Gamma,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.GeV),
 clsim.I3CLSimParticleParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.Brems,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.GeV),
 clsim.I3CLSimParticleParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.DeltaE,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.GeV),
 clsim.I3CLSimParticleParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.PairProd,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.GeV)
]

tray.AddModule("I3Reader","reader",
               Filename=options.INFILE)

tray.AddModule("I3CLSimModule", "clsim",
               RandomService=randomService,
               MediumProperties=mediumProperties,
               ParameterizationList=parameterizations,
               MaxNumParallelEvents=100,
               IgnoreMuons=False)

#tray.AddModule("Dump","dumper")

tray.AddModule("I3Writer","writer",
    Filename = options.OUTFILE)

tray.AddModule("TrashCan", "the can")

tray.Execute()
tray.Finish()
