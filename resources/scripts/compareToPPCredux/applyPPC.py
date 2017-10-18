#!/usr/bin/env python

from __future__ import print_function
from optparse import OptionParser

usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)
parser.add_option("-o", "--outfile",
                  dest="OUTFILE", help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_option("-i", "--infile",
                  dest="INFILE", help="Read input from INFILE (.i3{.gz} format)")
parser.add_option("-s", "--seed",type="int",default=12345,
                                  dest="SEED", help="Initial seed for the random number generator")
parser.add_option("-r", "--runnumber", type="int", default=1,
                  dest="RUNNUMBER", help="The run number for this simulation")

parser.add_option("--icemodel", default="test_ice_models/lea",
                  dest="ICEMODEL", help="A clsim ice model file/directory")

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
    print("No OUTFILE specified!")
    parser.print_help()
    exit(-1)

if not options.INFILE:
    print("No INFILE specified!")
    parser.print_help()
    exit(-1)


from I3Tray import *
from os.path import expandvars
import os
import sys

from icecube import icetray, dataclasses, dataio, phys_services

os.putenv("PPCTABLESDIR", options.ICEMODEL)


#load("libcudart")
load("libxppc")
load("libppc")
load("libppc-eff")

tray = I3Tray()

tray.AddService("I3GSLRandomServiceFactory", "random",
    Seed = options.SEED,
)

tray.AddModule("I3Reader","reader",
               Filename=options.INFILE)

#tray.AddModule("Dump", "dump")

tray.AddModule("i3ppc", "ppc",
    gpu = 0,
    bad = [],
    #JPALpulses=False,
)

# WARNING: this will adjust *all* the MCHitSeriesMaps in the frame.
tray.AddModule("AdjEff", "ppc-eff",
    eff = 0.9,
    )

tray.AddModule("I3Writer","writer",
    Filename = options.OUTFILE)



tray.Execute()

