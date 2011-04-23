#!/usr/bin/env python

from optparse import OptionParser

usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)
parser.add_option("-o", "--outfile", default="test_muons_photons.i3",
                  dest="OUTFILE", help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_option("-i", "--infile", default="test_muons.i3",
                  dest="INFILE", help="Read input from INFILE (.i3{.gz} format)")
parser.add_option("-s", "--seed",type="int",default=12345,
                                  dest="SEED", help="Initial seed for the random number generator")
parser.add_option("-r", "--runnumber", type="int", default=1,
                  dest="RUNNUMBER", help="The run number for this simulation")
parser.add_option("-p", "--max-parallel-events", type="int", default=100,
                  dest="MAXPARALLELEVENTS", help="maximum number of events(==frames) that will be processed in parallel")
parser.add_option("-m", "--oversize-factor", type="float", default=1.,
                  dest="OVERSIZEFACTOR", help="scale the OM radius by this factor, but scale back the OM acceptance accordingly")
parser.add_option("--apply-mmc", action="store_true", default=False,
                  dest="APPLYMMC", help="apply MMC to the I3MCTree before passing it to CLSim")
parser.add_option("--mmc-with-recc", action="store_true", default=False,
                  dest="MMCWITHRECC", help="add the -recc and -cont options to MMC to make it output muon slices with energy losses taken into account")
parser.add_option("--chop-muons", action="store_const", default=-1, const=1,
                  dest="CHOPMUONS", help="Tries to estimate the muon energy between each pair of cascades along its track")
parser.add_option("--no-chop-muons", action="store_const", default=-1, const=0,
                  dest="CHOPMUONS", help="Tries to estimate the muon energy between each pair of cascades along its track")

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

if options.MMCWITHRECC and (not options.APPLYMMC):
    print "using the --mmc-with-recc without --apply-mmc will have no effect"

if options.CHOPMUONS==-1:
    if options.MMCWITHRECC and options.APPLYMMC:
        options.CHOPMUONS=False
        print "auto-configured --chop-muons=False"
    else:
        options.CHOPMUONS=True
        print "auto-configured --chop-muons=True"
else:
    options.CHOPMUONS = (options.CHOPMUONS==1)
    
    if options.APPLYMMC:
        if options.CHOPMUONS and options.MMCWITHRECC:
            print "you cannot use the --chop-muons and the --mmc-with-recc option together"
            exit(-1)
        elif (not options.CHOPMUONS) and (not options.MMCWITHRECC):
            print "you should consider using either the --chop-muons or --mmc-with-recc options"


if options.APPLYMMC:
    if  options.MMCWITHRECC:
        print "applying MMC (with -recc -cont)"
    else:
        print "applying MMC (without -recc -cont)"
if options.CHOPMUONS:
    print "chopping muons"
else:
    print "not chopping muons"


from I3Tray import *
from os.path import expandvars
import os
import sys

from icecube import icetray, dataclasses, dataio, phys_services
from icecube import clsim

if options.APPLYMMC:
    load("libc2j-icetray")
    load("libmmc-icetray")
    MMCseed=options.SEED

radiusOverSizeFactor=options.OVERSIZEFACTOR

if radiusOverSizeFactor != 1.:
    print "using a OM radius oversize factor of {0}".format(radiusOverSizeFactor)

tray = I3Tray()

if options.APPLYMMC:
    tray.AddService("I3JavaVMFactory","javavm",
                    options = [expandvars("-Djava.class.path=$I3_BUILD/lib/mmc.jar"), "-server", "-Xms64m", "-Xmx512m"])

# a random number generator
randomService = phys_services.I3SPRNGRandomService(
    seed = options.SEED,
    nstreams = 10000,
    streamnum = options.RUNNUMBER)

# ice properties (SPICE-Mie model)
mediumProperties = clsim.MakeIceCubeMediumProperties()

domAcceptance = clsim.GetIceCubeDOMAcceptance(domRadius = 0.16510*I3Units.m*radiusOverSizeFactor)
domAngularSensitivity = clsim.GetIceCubeDOMAngularSensitivity(holeIce=True)

# parameterizations for fast simulation (bypassing Geant4)
# converters first:
cascadeConverter = clsim.I3CLSimParticleToStepConverterCascadeParameterization(randomService=randomService, photonsPerStep=200)

# now set up a list of converters with particle types and valid energy ranges
parameterizationsMuon = [
 clsim.I3CLSimParticleParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.MuMinus,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.PeV,
                                       needsLength=True),
 clsim.I3CLSimParticleParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.MuPlus,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.PeV,
                                       needsLength=True)
]

parameterizationsOther = [
 ## do we need some special handling for neutrons?
 #clsim.I3CLSimParticleParameterization(converter=cascadeConverter,
 #                                      forParticleType=dataclasses.I3Particle.Neutron,
 #                                      fromEnergy=0.0*I3Units.GeV,
 #                                      toEnergy=1000.*I3Units.GeV),

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
                                       forParticleType=dataclasses.I3Particle.EMinus,
                                       fromEnergy=0.5*I3Units.GeV,
                                       toEnergy=1000.*I3Units.GeV),
 clsim.I3CLSimParticleParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.EPlus,
                                       fromEnergy=0.5*I3Units.GeV,
                                       toEnergy=1000.*I3Units.GeV),
 clsim.I3CLSimParticleParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.Gamma,
                                       fromEnergy=0.5*I3Units.GeV,
                                       toEnergy=1000.*I3Units.GeV),

 clsim.I3CLSimParticleParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.Brems,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.PeV),
 clsim.I3CLSimParticleParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.DeltaE,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.PeV),
 clsim.I3CLSimParticleParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.PairProd,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.PeV),
 clsim.I3CLSimParticleParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.NuclInt,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.PeV),
                                        
]

tray.AddModule("I3Reader","reader",
               Filename=options.INFILE)

if options.APPLYMMC:
    mmcOpts = "-seed=%i -radius=900 -length=1600" % (MMCseed)
    if options.MMCWITHRECC:
        mmcOpts = "-cont -recc " + mmcOpts
    
    tray.AddModule("I3PropagatorMMC","propagate",
                   PrimaryTreeName = "I3MCTree",
                   mode=-1,
                   opts=mmcOpts,
                   ShiftParticles = False,
                   )

if options.CHOPMUONS:
    tray.AddModule("I3MuonSlicer", "chopMuons",
                   InputMCTreeName="I3MCTree",
                   MMCTrackListName="MMCTrackList",
                   OutputMCTreeName="I3MCTree_sliced")
    clSimMCTreeName = "I3MCTree_sliced"
else:
    clSimMCTreeName = "I3MCTree"

tray.AddModule("I3CLSimModule", "clsim",
               MCTreeName=clSimMCTreeName,
               DOMRadius = 0.16510*I3Units.m*radiusOverSizeFactor, # 13" diameter
               RandomService=randomService,
               MediumProperties=mediumProperties,
               IgnoreNonIceCubeOMNumbers=True, # ignore AMANDA and IceTop OMKeys (do NOT use for any other detector!)
               GenerateCherenkovPhotonsWithoutDispersion=False,
               WavelengthGenerationBias=domAcceptance,
               ParameterizationList=parameterizationsMuon+parameterizationsOther,
               #ParameterizationList=parameterizationsMuon,
               MaxNumParallelEvents=options.MAXPARALLELEVENTS,
               
               #OpenCLPlatformName="NVIDIA CUDA",
               #OpenCLDeviceName="GeForce GTX 580",
               #OpenCLUseNativeMath=True,
               #OpenCLApproximateNumberOfWorkItems=512000,
               
               #OpenCLPlatformName="ATI Stream",
               #OpenCLDeviceName="Intel(R) Core(TM) i5 CPU         760  @ 2.80GHz",
               #OpenCLUseNativeMath=False,
               #OpenCLApproximateNumberOfWorkItems=51200,

               #OpenCLPlatformName="Apple",
               #OpenCLDeviceName="GeForce 9600M GT",
               #OpenCLUseNativeMath=True,
               #OpenCLApproximateNumberOfWorkItems=1024,

               #OpenCLPlatformName="Apple",
               #OpenCLDeviceName="GeForce 9400M",
               #OpenCLUseNativeMath=True,
               #OpenCLApproximateNumberOfWorkItems=512,
               
               OpenCLPlatformName="Apple",
               OpenCLDeviceName="Intel(R) Core(TM)2 Duo CPU     T9600  @ 2.80GHz",
               OpenCLUseNativeMath=False,
               OpenCLApproximateNumberOfWorkItems=51200,
               
               )

tray.AddModule("I3PhotonToMCHitConverter", "make_hits",
               RandomService = randomService,
               MCTreeName=clSimMCTreeName,
               InputPhotonSeriesMapName = "PropagatedPhotons",
               OutputMCHitSeriesMapName = "MCHitSeriesMap_clsim",
               DOMRadiusWithoutOversize=0.16510*I3Units.m,
               DOMOversizeFactor=radiusOverSizeFactor,
               WavelengthAcceptance = domAcceptance,
               AngularAcceptance = domAngularSensitivity)

tray.AddModule("I3Writer","writer",
    Filename = options.OUTFILE)

tray.AddModule("TrashCan", "the can")

tray.Execute()
tray.Finish()
