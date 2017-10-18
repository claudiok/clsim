#!/usr/bin/env python

from __future__ import print_function
from optparse import OptionParser
import os
import string

usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)
parser.add_option("-o", "--outfile", default=None,
                  dest="OUTFILE", help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_option("-i", "--infile", default="test_flashesSC.i3",
                  dest="INFILE", help="Read input from INFILE (.i3{.gz} format)")
parser.add_option("-s", "--seed",type="int",default=12346,
                  dest="SEED", help="Initial seed for the random number generator")
parser.add_option("-r", "--runnumber", type="int", default=1,
                  dest="RUNNUMBER", help="The run number for this simulation")
parser.add_option("-p", "--max-parallel-events", type="int", default=10,
                  dest="MAXPARALLELEVENTS", help="maximum number of events(==frames) that will be processed in parallel")
parser.add_option("--keep-photon-data", action="store_false", default=True,
                  dest="REMOVEPHOTONDATA", help="Keep I3Photons before writing the output file (in addition to I3MCPEs")

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
                raise RuntimeError("cannot find input file!")
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

if options.REMOVEPHOTONDATA:
    print("not storing I3Photons")
else:
    print("storing I3Photons")


from I3Tray import *
from os.path import expandvars
import os
import sys

from icecube import icetray, dataclasses, dataio, phys_services
from icecube import clsim

# a random number generator
randomService = phys_services.I3SPRNGRandomService(
    seed = options.SEED,
    nstreams = 10000,
    streamnum = options.RUNNUMBER)

tray = I3Tray()


tray.AddModule("I3Reader","reader",
               Filename=infile)


if options.REMOVEPHOTONDATA:
    photonSeriesName = None
else:
    photonSeriesName = "PropagatedPhotons"

tray.AddSegment(clsim.I3CLSimMakeHits, "makeCLSimHits",
    PhotonSeriesName = photonSeriesName,
    ParallelEvents = options.MAXPARALLELEVENTS,
    RandomService = randomService,
    UseGPUs=False,
    UseCPUs=True,
    IceModelLocation=expandvars("$I3_SRC/clsim/resources/ice/spice_mie"),
    FlasherPulseSeriesName="SCFlashes")

tray.AddModule("I3Writer","writer",
    Filename = outdir+outfile)

tray.AddModule("TrashCan", "the can")

tray.Execute()

