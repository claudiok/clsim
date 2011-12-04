#!/usr/bin/env python
from os.path import expandvars
from optparse import OptionParser

usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)
parser.add_option("-o", "--outfile",default="PINGU_proton_decay.i3",
                  dest="OUTFILE", help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_option("-g", "--gcd",default=expandvars("$I3_PORTS/test-data/sim/GeoCalibDetectorStatus_IC86.55380_corrected.i3.gz"),
                  dest="GCDFILE", help="Read geometry from GCDFILE (.i3{.gz} format)")
parser.add_option("-s", "--seed",type="int",default=12345,
                  dest="SEED", help="Initial seed for the random number generator")
parser.add_option("-r", "--runnumber", type="int", default=1,
                  dest="RUNNUMBER", help="The run number for this simulation")
parser.add_option("-n", "--numevents", type="int", default=10000,
                  dest="NUMEVENTS", help="The number of events per run")

# parse cmd line args, bail out if anything is not understood
(options,args) = parser.parse_args()
if len(args) != 0:
        crap = "Got undefined options:"
        for a in args:
                crap += a
                crap += " "
        parser.error(crap)

###############################

from I3Tray import *

from icecube import icetray, dataclasses, dataio, phys_services, sim_services
from icecube import clsim

from proton_decay_generator import I3SimpleProtonDecayGenerator

import math, numpy

# a random number generator
randomService = phys_services.I3SPRNGRandomService(
    seed = options.SEED,
    nstreams = 10000,
    streamnum = options.RUNNUMBER)

tray = I3Tray()

# make a stream of events
tray.AddModule("I3InfiniteSource","streams",
               Prefix=options.GCDFILE,
               Stream=icetray.I3Frame.DAQ)
tray.AddModule("I3MCEventHeaderGenerator","gen_header",
               Year=2009,
               DAQTime=158100000000000000,
               RunNumber=options.RUNNUMBER,
               EventID=1,
               IncrementEventID=True)
tray.AddModule("Delete", "cleanup", Keys=["MCTimeIncEventID"])

# decay some protons
tray.AddModule(I3SimpleProtonDecayGenerator, "I3SimpleProtonDecayGenerator",
               posZRange=(-500.*I3Units.m,-100.*I3Units.m),
               radius=60.*I3Units.m,
               centerX=45.771827*I3Units.m,
               centerY=-34.411053*I3Units.m,
               randomService=randomService,
               Streams=[icetray.I3Frame.DAQ])

# propagate particles, generate C'kov photons, track them and create hits in IceCube DOMs
tray.AddSegment(clsim.I3CLSimMakeHits, "makeCLSimHits",
                ParallelEvents = 2000,     # this assumes low energies, you will probably run out of memory if your energies are too high
                RandomService = randomService,
                UseGPUs=False,
                UseCPUs=True, 
                UseGeant4=True,         # use Geant4 (i.e. no parameterizations)
                MMCTrackListName=None,  # no MMC!
                IceModelLocation=expandvars("$I3_SRC/clsim/resources/ice/spice_mie"))

# write the results to disk
tray.AddModule("I3Writer", "writer",
    filename = options.OUTFILE)

tray.AddModule("TrashCan","trash")

tray.Execute(options.NUMEVENTS+3)
tray.Finish()
