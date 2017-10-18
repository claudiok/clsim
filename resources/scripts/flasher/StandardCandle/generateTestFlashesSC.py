#!/usr/bin/env python

from optparse import OptionParser
from os.path import expandvars

usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)
parser.add_option("-o", "--outfile",default="test_flashesSC.i3",
                  dest="OUTFILE", help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_option("-s", "--seed",type="int",default=12344,
                  dest="SEED", help="Initial seed for the random number generator")
parser.add_option("-g", "--gcd",default=expandvars("$I3_TESTDATA/sim/GeoCalibDetectorStatus_IC86.55380_corrected.i3.gz"),
                  dest="GCDFILE", help="Read geometry from GCDFILE (.i3{.gz} format)")
parser.add_option("-r", "--runnumber", type="int", default=1,
                  dest="RUNNUMBER", help="The run number for this simulation")
parser.add_option("-n", "--numevents", type="int", default=1,
                  dest="NUMEVENTS", help="The number of events per run")

# parse cmd line args, bail out if anything is not understood
(options,args) = parser.parse_args()
if len(args) != 0:
        crap = "Got undefined options:"
        for a in args:
                crap += a
                crap += " "
        parser.error(crap)

from I3Tray import *
import os
import sys

from icecube import icetray, dataclasses, dataio, phys_services, clsim, sim_services

import math
import numpy


tray = I3Tray()

# a random number generator
randomService = phys_services.I3SPRNGRandomService(
    seed = options.SEED,
    nstreams = 10000,
    streamnum = options.RUNNUMBER)

tray.AddModule("I3InfiniteSource","streams",
               Prefix=options.GCDFILE,
               Stream=icetray.I3Frame.DAQ)

tray.AddModule("I3MCEventHeaderGenerator","gen_header",
               Year=2012,
               DAQTime=7968509615844458,
               RunNumber=1,
               EventID=1,
               IncrementEventID=True)

tray.AddModule(clsim.StandardCandleFlasherPulseSeriesGenerator, "StandardCandleFlasherPulseSeriesGenerator",
               FlasherPulseSeriesName = "SCFlashes",
               PhotonsPerPulse = 1.7e10,    # @  0.5% nominal output (SC2) [see http://wiki.icecube.wisc.edu/index.php/Standard_Candle#Data_runs]
               # PhotonsPerPulse = 3.5e11,  # @  1.0% nominal output (SC2)
               # PhotonsPerPulse = 7.3e11,  # @  3.0% nominal output (SC2)
               # PhotonsPerPulse = 2.1e12,  # @ 10.0% nominal output (SC2)
               # PhotonsPerPulse = 6.5e12,  # @ 30.5% nominal output (SC2)
               # PhotonsPerPulse = 9.3e12,  # @ 51.0% nominal output (SC2)
               # PhotonsPerPulse = 2.5e13,  # @  100% nominal output (SC2)
               # PhotonsPerPulse = 2.3e11,  # @  5.0% nominal output (SC1)
               # PhotonsPerPulse = 4.4e11,  # @ 10.0% nominal output (SC1)
               # PhotonsPerPulse = 1.29e12, # @ 30.0% nominal output (SC1)
               # PhotonsPerPulse = 2.2e12,  # @ 50.0% nominal output (SC1)
               # PhotonsPerPulse = 4.0e12,  # @  100% nominal output (SC1)
               FlashTime = 0.*I3Units.ns, # this time is arbitrary (the trigger modules should time-shift this anyway)
               CandleNumber = 2)          # simulate SC2


tray.AddModule("I3Writer","writer",
    Filename = options.OUTFILE)



tray.Execute(options.NUMEVENTS+3)




