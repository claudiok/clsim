#!/usr/bin/env python

from __future__ import print_function

from optparse import OptionParser
import os

usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)
parser.add_option("-o", "--outfile", default="test_muons_photons.i3",
                  dest="OUTFILE", help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_option("-i", "--infile", default="test_muons.i3",
                  dest="INFILE", help="Read input from INFILE (.i3{.gz} format)")
parser.add_option("-s", "--seed",type="int",
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
parser.add_option("--remove-photon-data", action="store_true", default=False,
                  dest="REMOVEPHOTONDATA", help="Remove I3Photons before writing the output file (only keep hits)")

# parse cmd line args, bail out if anything is not understood
(options,args) = parser.parse_args()
if len(args) != 0:
        crap = "Got undefined options:"
        for a in args:
                crap += a
                crap += " "
        parser.error(crap)

########################
if options.INFILE:
        filename = options.INFILE
        if os.access(filename,os.R_OK) == False:
                raise IOError("cannot find input file!")
        infile = filename
        print('using input file %s' % infile)
else:
        print("No input file!")
        parser.print_help()
        exit(-1)

infileRoot, infileExt = os.path.splitext(infile)
if infileExt == ".gz":
    infileRoot2, infileExt2 = os.path.splitext(infileRoot)
    if infileExt2 == ".i3":
        infileRoot=infileRoot2
        infileExt = ".i3.gz"

if infileExt != ".i3" and infileExt != ".i3.gz":
        raise Exception("you have to specify either a .i3 or an .i3.gz file!")

########################
outdir=""
outfile=None
if options.OUTFILE:
        outfile = options.OUTFILE
        # did the user specify a directory? then use that and auto-generate
        if os.path.isdir(outfile):
            outdir = outfile
            outfile = None
        else:
            outdir, outfile = os.path.split(outfile)

# add a trailing slash to the output directory name if not already there
if outdir and outdir!="":
    if outdir[-1] != "/":
        outdir += "/"

if not outfile:
        # automatically generate the output filename
        infileRootDir, infileRootFile = os.path.split(infileRoot)
        outfile = infileRootFile + "_clsim"
        outfile = outfile + infileExt
print("output dir is %s" % outdir)
print("output file is %s" % outdir + outfile)

########################



if options.MMCWITHRECC and (not options.APPLYMMC):
    print("using the --mmc-with-recc without --apply-mmc will have no effect")

if options.CHOPMUONS==-1:
    if options.MMCWITHRECC and options.APPLYMMC:
        options.CHOPMUONS=False
        print("auto-configured --chop-muons=False")
    else:
        options.CHOPMUONS=True
        print("auto-configured --chop-muons=True")
else:
    options.CHOPMUONS = (options.CHOPMUONS==1)
    
    if options.APPLYMMC:
        if options.CHOPMUONS and options.MMCWITHRECC:
            print("you cannot use the --chop-muons and the --mmc-with-recc option together")
            exit(-1)
        elif (not options.CHOPMUONS) and (not options.MMCWITHRECC):
            print("you should consider using either the --chop-muons or --mmc-with-recc options")


if options.APPLYMMC:
    if  options.MMCWITHRECC:
        print("applying MMC (with -recc -cont)")
    else:
        print("applying MMC (without -recc -cont)")
if options.CHOPMUONS:
    print("chopping muons")
else:
    print("not chopping muons")


if options.REMOVEPHOTONDATA:
    print("not storing I3Photons")
else:
    print("storing I3Photons")


if options.SEED is None:
    # this is not the best possible option, but at least it's not a fixed seed
    theSeed = hash(options.INFILE)
    if theSeed < 0: theSeed = -theSeed
    theSeed = theSeed % 100000
    print("using auto-seed generated from input filename:", theSeed)
else:
    theSeed = options.SEED
    

from I3Tray import *
from os.path import expandvars
import os
import sys

from icecube import icetray, dataclasses, dataio, phys_services
from icecube import clsim

if options.APPLYMMC:
    load("libc2j-icetray")
    load("libmmc-icetray")
    MMCseed=theSeed

radiusOverSizeFactor=options.OVERSIZEFACTOR

if radiusOverSizeFactor != 1.:
    print("using a OM radius oversize factor of {0}".format(radiusOverSizeFactor))

tray = I3Tray()

if options.APPLYMMC:
    tray.AddService("I3JavaVMFactory","javavm",
                    options = [expandvars("-Djava.class.path=$I3_BUILD/lib/mmc.jar"), "-server", "-Xms64m", "-Xmx512m"])

# a random number generator
randomService = phys_services.I3SPRNGRandomService(
    seed = theSeed,
    nstreams = 10000,
    streamnum = options.RUNNUMBER)

# water properties (partic-0.0075)
mediumProperties = clsim.MakeAntaresMediumProperties()

domAcceptance = clsim.GetAntaresOMAcceptance(domRadius = (17./2.) * 0.0254*I3Units.m*radiusOverSizeFactor)
domAngularSensitivity = clsim.GetAntaresOMAngularSensitivity(name='Spring09')

# parameterizations for fast simulation (bypassing Geant4)
# converters first:
cascadeConverter = clsim.I3CLSimLightSourceToStepConverterPPC(photonsPerStep=200, highPhotonsPerStep=2000, useHighPhotonsPerStepStartingFromNumPhotons=1e8)

# now set up a list of converters with particle types and valid energy ranges
parameterizationsMuon = [
 clsim.I3CLSimLightSourceParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.MuMinus,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.PeV,
                                       needsLength=True),
 clsim.I3CLSimLightSourceParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.MuPlus,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.PeV,
                                       needsLength=True)
]

parameterizationsOther = [
 ## do we need some special handling for neutrons?
 clsim.I3CLSimLightSourceParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.Neutron,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.PeV),

 clsim.I3CLSimLightSourceParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.Pi0,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.PeV),
 clsim.I3CLSimLightSourceParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.PiPlus,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.PeV),
 clsim.I3CLSimLightSourceParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.PiMinus,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.PeV),
 clsim.I3CLSimLightSourceParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.K0_Long,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.PeV),
 clsim.I3CLSimLightSourceParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.KPlus,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.PeV),
 clsim.I3CLSimLightSourceParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.KMinus,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.PeV),
 clsim.I3CLSimLightSourceParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.PPlus,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.PeV),
 clsim.I3CLSimLightSourceParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.PMinus,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.PeV),
 clsim.I3CLSimLightSourceParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.K0_Short,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.PeV),

 clsim.I3CLSimLightSourceParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.EMinus,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.PeV),
 clsim.I3CLSimLightSourceParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.EPlus,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.PeV),
 clsim.I3CLSimLightSourceParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.Gamma,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.PeV),

 clsim.I3CLSimLightSourceParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.Brems,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.PeV),
 clsim.I3CLSimLightSourceParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.DeltaE,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.PeV),
 clsim.I3CLSimLightSourceParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.PairProd,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.PeV),
 clsim.I3CLSimLightSourceParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.NuclInt,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.PeV),
 clsim.I3CLSimLightSourceParameterization(converter=cascadeConverter,
                                       forParticleType=dataclasses.I3Particle.Hadrons,
                                       fromEnergy=0.0*I3Units.GeV,
                                       toEnergy=1000.*I3Units.PeV),

 # catch all. nothing should go to Geant4
 clsim.I3CLSimLightSourceParameterization(converter=cascadeConverter,
                                       forParticleType=clsim.I3CLSimLightSourceParameterization.AllParticles_t(),
                                       fromEnergy=0.,
                                       toEnergy=float('inf'))
                                        
]

tray.AddModule("I3Reader","reader",
               Filename=infile)

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

openCLDevices = [device for device in clsim.I3CLSimOpenCLDevice.GetAllDevices() if device.gpu]
for device in openCLDevices:
    if string.count(device.device, 'Tesla') > 0 or string.count(device.device, 'GTX') > 0:
        device.useNativeMath=True
        device.approximateNumberOfWorkItems=1024000
    else:
        device.useNativeMath=False
        device.approximateNumberOfWorkItems=10240

tray.AddModule("I3CLSimModule", "clsim",
               MCTreeName=clSimMCTreeName,
               DOMRadius = (17./2.) * 0.0254*I3Units.m*radiusOverSizeFactor, # 13" diameter
               RandomService=randomService,
               MediumProperties=mediumProperties,
               IgnoreNonIceCubeOMNumbers=False, # ignore AMANDA and IceTop OMKeys (do NOT use for any other detector!)
               GenerateCherenkovPhotonsWithoutDispersion=False,
               WavelengthGenerationBias=domAcceptance,
               ParameterizationList=parameterizationsMuon+parameterizationsOther,
               #ParameterizationList=parameterizationsMuon,
               MaxNumParallelEvents=options.MAXPARALLELEVENTS,
               
               OpenCLDeviceList=openCLDevices
               )

#pmtPhotonSimulator = clsim.I3CLSimPMTPhotonSimulatorIceCube(randomService=randomService, jitter=2.*I3Units.ns)
tray.AddModule("I3PhotonToMCHitConverter", "make_hits",
               RandomService = randomService,
               MCTreeName = clSimMCTreeName,
               InputPhotonSeriesMapName = "PropagatedPhotons",
               OutputMCHitSeriesMapName = "MCHitSeriesMap",
               DOMRadiusWithoutOversize=(17./2.) * 0.0254*I3Units.m,
               DOMOversizeFactor = radiusOverSizeFactor,
               WavelengthAcceptance = domAcceptance,
               AngularAcceptance = domAngularSensitivity,
               PMTPhotonSimulator = None, # no jitter or after-pulse simulation here for Antares
               IgnoreDOMsWithoutDetectorStatusEntry = False,
               OnlyWarnAboutInvalidPhotonPositions=True)

if options.REMOVEPHOTONDATA:
    tray.AddModule("Delete", "delete_photons",
        Keys = ["PropagatedPhotons"])

tray.AddModule("I3Writer","writer",
    Filename = outdir+outfile)



tray.Execute()

