#!/usr/bin/env python

from optparse import OptionParser
from os.path import expandvars

usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)
parser.add_option("-o", "--outfile",default="test_flashes.i3",
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

from icecube import icetray, dataclasses, dataio, phys_services, sim_services, clsim

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

tray.AddModule(clsim.FakeFlasherInfoGenerator, "FakeFlasherInfoGenerator",
               FlashingDOM = icetray.OMKey(57,10),
               #FlashingDOM = icetray.OMKey(36,22), # a cDOM
               FlasherTime = 0.*I3Units.ns,
               # FlasherMask = 0b111111000000, # only the 6 horizontal LEDs,  0-5 are tilted, 6-11 are horizontal [cDOM LEDs are all horizantal]
               FlasherMask = 0b000000111111, # only the 6 tilted LEDs
               #FlasherMask = 0b10101, # 505nm LEDs only (on a cDOM) 

               # FlasherBrightness = 127, # full brightness
               # FlasherWidth = 127,      # full width
               FlasherBrightness = 15,
               FlasherWidth = 31,

               )


tray.AddModule("I3Writer","writer",
    Filename = options.OUTFILE)



tray.Execute(options.NUMEVENTS+3)




