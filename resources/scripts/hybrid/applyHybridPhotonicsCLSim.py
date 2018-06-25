#!/usr/bin/env python
"""

"""

from __future__ import print_function
from optparse import OptionParser
import os
import string

usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)
parser.add_option("-o", "--outfile", default=None,
                  dest="OUTFILE", help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_option("-i", "--infile", default=None,
                  dest="INFILE", help="Read input from INFILE (.i3{.gz} format)")
parser.add_option("-s", "--seed",type="int",default=12345,
                  dest="SEED", help="Initial seed for the random number generator")
parser.add_option("-r", "--runnumber", type="int", default=1,
                  dest="RUNNUMBER", help="The run number for this simulation")
parser.add_option("-p", "--max-parallel-events", type="int", default=100,
                  dest="MAXPARALLELEVENTS", help="maximum number of events(==frames) that will be processed in parallel")
parser.add_option("--apply-mmc", action="store_true", default=False,
                  dest="APPLYMMC", help="apply MMC to the I3MCTree before passing it to CLSim")

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

from I3Tray import *
from os.path import expandvars
import os
import sys


# pull in all modules
from icecube import icetray, dataclasses, dataio
from icecube import phys_services, sim_services
from icecube import photonics_service
from icecube import clsim


# set up a random number generator
randomService = phys_services.I3SPRNGRandomService(
    seed = options.SEED,
    nstreams = 10000,
    streamnum = options.RUNNUMBER)


# additional things that need to be loaded for MMC
if options.APPLYMMC:
    load("libc2j-icetray")
    load("libmmc-icetray")
    MMCseed=options.SEED



# the tray!
tray = I3Tray()



# read the input file
tray.AddModule("I3Reader","reader",
               Filename=infile)



# optionally apply MMC
if options.APPLYMMC:
    mmcOpts = "-seed=%i -radius=900 -length=1600" % (MMCseed)
    
    tray.AddModule("I3PropagatorMMC","propagate",
                   PrimaryTreeName = "I3MCTree",
                   mode=-1,
                   opts=mmcOpts,
                   ShiftParticles = False,
                   )



# split the MCTree into a cascade-only and a track-only version
tray.AddModule("I3MCTreeHybridSimulationSplitter", "splitMCTree",
    InputMCTreeName="I3MCTree",
    OutputMCTreeNameTracks="I3MCTreeTracks",
    OutputMCTreeNameCascades="I3MCTreeCascades")



# simulate cascades (with photonics-hit-maker)
cascade_service = photonics_service.I3PhotoSplineService(
    amplitudetable='/Users/claudio/Documents/Work/IceTray/PhotonTables/spline-tables/ems_spice1_z20_a10.abs.fits',
    timingtable='/Users/claudio/Documents/Work/IceTray/PhotonTables/spline-tables/ems_spice1_z20_a10.prob.fits',
    timingSigma=0.)

tray.AddModule("I3PhotonicsHitMaker", "hitsFromTheTable",
    CascadeService = cascade_service,
    TrackService = None, # tracks are handled by clsim
    UnshadowedFraction = 0.9,
    Input = "I3MCTreeCascades",
    Output = "I3MCPESeriesMapCascades",
    RandomService = randomService
    )



# simulate tracks (with clsim)
tray.AddSegment(clsim.I3CLSimMakeHits, "makeCLSimHits",
    PhotonSeriesName = None,
    MCTreeName = "I3MCTreeTracks",
    MCPESeriesName = "I3MCPESeriesMapTracks",
    MMCTrackListName = "MMCTrackList",
    ParallelEvents = options.MAXPARALLELEVENTS,
    RandomService = randomService,
    # DoNotParallelize=True, # you may need to turn this on for clusters that assume "1 job == 1 core"
    UseGeant4=False, # never use this with Geant4!
    UseGPUs=False,
    UseCPUs=True,
    IceModelLocation=expandvars("$I3_BUILD/clsim/resources/ice/spice_1"),
    DisableTilt=True, # cannot use tilted ice layers with tables (tables cannot describe tilt)
    )

tray.AddModule("Delete", "cleanup_clsim_sliced_MCTree",
    Keys = ["I3MCTreeTracks_sliced"])


# combine the resulting I3MCPESeriesMaps
tray.AddModule("I3CombineMCPE", "combine_hits",
    InputResponses = ["I3MCPESeriesMapTracks", "I3MCPESeriesMapCascades"],
    OutputResponse = "I3MCPESeriesMap")



# delete the original maps and the split I3MCTrees
tray.AddModule("Delete", "cleanup_hitseriesmaps",
    Keys = ["I3MCPESeriesMapTracks", "I3MCPESeriesMapCascades"])

tray.AddModule("Delete", "cleanup_MCTree",
    Keys=["I3MCTreeTracks", "I3MCTreeCascades"])



# write the results to a file
tray.AddModule("I3Writer","writer",
    Filename = outdir+outfile)



tray.Execute()

