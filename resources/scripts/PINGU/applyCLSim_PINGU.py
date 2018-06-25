#!/usr/bin/env python

from __future__ import print_function

from optparse import OptionParser
import os
import string

usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)
parser.add_option("-o", "--outfile", default="test_muons_photons.i3",
                  dest="OUTFILE", help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_option("-i", "--infile", default="test_muons.i3",
                  dest="INFILE", help="Read input from INFILE (.i3{.gz} format)")
parser.add_option("-g", "--gcdfile", default=None,
                  dest="GCDFILE", help="Read an optional GCDFILE (.i3{.gz} format) before starting with the actual data")
parser.add_option("-s", "--seed",type="int",default=12345,
                                  dest="SEED", help="Initial seed for the random number generator")
parser.add_option("-r", "--runnumber", type="int", default=1,
                  dest="RUNNUMBER", help="The run number for this simulation")
parser.add_option("-p", "--max-parallel-events", type="int", default=100,
                  dest="MAXPARALLELEVENTS", help="maximum number of events(==frames) that will be processed in parallel")
parser.add_option("--apply-mmc", action="store_true", default=False,
                  dest="APPLYMMC", help="apply MMC to the I3MCTree before passing it to CLSim")
parser.add_option("--remove-photon-data", action="store_true", default=False,
                  dest="REMOVEPHOTONDATA", help="Remove I3Photons before writing the output file (only keep hits)")
parser.add_option("--qify", action="store_true", default=False,
                  dest="QIFY", help="convert a P-frame only file to Q-frame only before using it")
parser.add_option("--clean-input", action="store_true", default=False,
                  dest="CLEANINPUT", help="only read frame items relevant for simulation from the input file")

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
        raise RuntimeError("you have to specify either a .i3 or an .i3.gz file!")

if options.GCDFILE is not None:
    print("using GCD file", options.GCDFILE) 
    infiles = [options.GCDFILE, infile]
else:
    print("using no extra GCD file")
    infiles = [infile]

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


# tray to adapt to condor stupidity
if "_CONDOR_SLOT" in os.environ:
    if "CUDA_VISIBLE_DEVICES" in os.environ:
        print("running in CONDOR, but CUDA_VISIBLE_DEVICES is already set. no further configuration necessary.")
    else:
        condorSlotNumber = int(os.environ["_CONDOR_SLOT"])
        print("script seems to be running in condor (slot %u). auto-configuring CUDA_VISIBLE_DEVICES!" % condorSlotNumber)
        os.environ["CUDA_VISIBLE_DEVICES"] = str(condorSlotNumber-1)


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

if options.APPLYMMC:
    load("libc2j-icetray")
    load("libmmc-icetray")
    MMCseed=options.SEED

tray = I3Tray()


tray.AddModule("I3Reader","reader",
               FilenameList=infiles)

if options.QIFY:
    tray.AddModule("QConverter","qconvert")

if options.CLEANINPUT:
    def PythonKeep(frame, Keys):
        for key in frame.keys():
            if key not in Keys:
                del frame[key]
    tray.AddModule(PythonKeep, "keeper", 
                   Keys=["I3EventHeader", "I3MCWeightDict", "I3MCTree", "MMCTrackList"],
                   Streams=[icetray.I3Frame.DAQ])

if options.APPLYMMC:
    mmcOpts = "-seed=%i -radius=900 -length=1600" % (MMCseed)
    
    tray.AddModule("I3PropagatorMMC","propagate",
                   PrimaryTreeName = "I3MCTree",
                   mode=-1,
                   opts=mmcOpts,
                   ShiftParticles = False,
                   )



photonSeriesName = "PropagatedPhotons"


tray.AddSegment(clsim.I3CLSimMakePhotons, "makeCLSimPhotons",
                PhotonSeriesName = photonSeriesName,
                OutputMCTreeName="I3MCTree_sliced",
                ParallelEvents = options.MAXPARALLELEVENTS,
                RandomService = randomService,
                UseGPUs=True,
                UseCPUs=False,
                IceModelLocation=expandvars("$I3_BUILD/ice-models/resources/models/spice_mie"),
                ExtraArgumentsToI3CLSimModule={"EnableDoubleBuffering":True,            # will sligthly speed up things (and still work on Teslas)
                                               "UseHardcodedDeepCoreSubdetector":True}, # needed for PINGU in order to squeeze all geometry
                DOMOversizeFactor=1.,                                                   # PINGU needs this (no oversizing, DOMs on strings are too dense)
                )

@icetray.traysegment
def clsimMakePhotonsAndWrite(tray, name, 
                             Filename,
                             RandomService,
                             PhotonSeriesName,
                             DOMOversizeFactor=5.,
                             RemovePhotonData=False
                             ):
    tray.AddSegment(clsim.I3CLSimMakeHitsFromPhotons, name+"_MakeHitsFromPhotons",
                    MCTreeName="I3MCTree_sliced",
                    PhotonSeriesName=PhotonSeriesName,
                    RandomService=RandomService,
                    DOMOversizeFactor=DOMOversizeFactor
                    )
    if RemovePhotonData:
        tray.AddModule("Delete", name+"_DeletePhotons",
                       Keys = [photonSeriesName])
    tray.AddModule("I3Writer", name+"_writer",
                   Filename = outdir+outfile,
                   Streams=[icetray.I3Frame.DAQ])


tray.AddModule(clsim.AsyncTap, "clsimMakePhotonsAndWrite_Async", 
               BufferSize = options.MAXPARALLELEVENTS*3,
               Segment = clsimMakePhotonsAndWrite,
               Args = dict(Filename=outdir+outfile,
                           RandomService=randomService,
                           PhotonSeriesName=photonSeriesName,
                           RemovePhotonData=options.REMOVEPHOTONDATA,
                           DOMOversizeFactor=1.,       # PINGU needs this (no oversizing, DOMs on strings are too dense)
                           )
               )

#tray.AddSegment(clsimMakePhotonsAndWrite, "clsimMakePhotonsAndWrite",
#                Filename=outdir+outfile,
#                RandomService=randomService,
#                PhotonSeriesName=photonSeriesName,
#                RemovePhotonData=options.REMOVEPHOTONDATA,
#                DOMOversizeFactor=1.,       # PINGU needs this (no oversizing, DOMs on strings are too dense)
#                )



tray.Execute()

