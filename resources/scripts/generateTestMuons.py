#!/usr/bin/env python

from optparse import OptionParser

usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)
parser.add_option("-o", "--outfile",default="test_muons.i3",
                  dest="OUTFILE", help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_option("-s", "--seed",type="int",default=12345,
                  dest="SEED", help="Initial seed for the random number generator")
parser.add_option("-r", "--runnumber", type="int", default=1,
                  dest="RUNNUMBER", help="The run number for this simulation")
parser.add_option("-n", "--numevents", type="int", default=100,
                  dest="NUMEVENTS", help="The number of events per run")
parser.add_option("--apply-mmc", action="store_true", default=False,
                  dest="APPLYMMC", help="apply MMC to the I3MCTree after generating them")
parser.add_option("--mmc-with-recc", action="store_true", default=False,
                  dest="MMCWITHRECC", help="add the -recc and -cont options to MMC to make it output muon slices with energy losses taken into account")

# parse cmd line args, bail out if anything is not understood
(options,args) = parser.parse_args()
if len(args) != 0:
        crap = "Got undefined options:"
        for a in args:
                crap += a
                crap += " "
        parser.error(crap)

if options.MMCWITHRECC and (not options.APPLYMMC):
    print "using the --mmc-with-recc without --apply-mmc will have no effect"

from I3Tray import *
from os.path import expandvars
import os
import sys

from icecube import icetray, dataclasses, dataio, phys_services, sim_services

from icecube import simple_injector
from icecube.simple_injector.simple_muon import I3SimpleMuon

if options.APPLYMMC:
    load("libc2j-icetray")
    load("libmmc-icetray")
    MMCseed=options.SEED

tray = I3Tray()

if options.APPLYMMC:
    tray.AddService("I3JavaVMFactory","javavm",
                    options = [expandvars("-Djava.class.path=$I3_BUILD/lib/mmc.jar"), "-server", "-Xms64m", "-Xmx512m"])

# a random number generator
randomService = phys_services.I3SPRNGRandomService(
    seed = options.SEED,
    nstreams = 10000,
    streamnum = options.RUNNUMBER)

tray.AddService("I3ReaderServiceFactory", "gcd_reader",
    Filename = "GeoCalibDetectorStatus_IC86.55040_official.i3.gz",
    OmitGeometry=False,
    OmitCalibration=False,
    OmitStatus=False,
    OmitEvent=True)

tray.AddService("I3MCTimeGeneratorServiceFactory", "events",
    Year=2009,
    DAQTime=179500000000000000,
    RunNumber=1,
    EventID=1,
    IncrementEventID=True)

tray.AddModule("I3Muxer", "muxme")

tray.AddModule(I3SimpleMuon, "injectMuon",
               I3RandomService = randomService,
               Type = dataclasses.I3Particle.ParticleType.MuMinus,
               NEvents = options.NUMEVENTS,
               EnergyMin = 1.*I3Units.TeV,
               EnergyMax = .1*I3Units.PeV,
               DiskRadius = 800.*I3Units.m,
               SphereRadius = 800.*I3Units.m
               )

if options.APPLYMMC:
    mmcOpts = "-seed=%i -radius=900 -length=1600" % (MMCseed)
    if options.MMCWITHRECC:
        mmcOpts = "-cont -recc " + mmcOpts
    
    tray.AddModule("I3PropagatorMMC","propagate",
                   PrimaryTreeName = "I3MCTree",
                   mode=-1,
                   opts=mmcOpts,
                   )
               

tray.AddModule("I3Writer","writer",
    Filename = options.OUTFILE)

tray.AddModule("TrashCan", "the can")

tray.Execute(options.NUMEVENTS+3)
tray.Finish()
