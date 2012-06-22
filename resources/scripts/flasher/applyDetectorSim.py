#!/usr/bin/env python
"""
This is an EXAMPLE script of how to apply detector simulation.
There does not seem to be a collection of simulation tray segments
anywhere, so this uses a bunch of modules copied from simprod
scripts. This simulation should be valid for IC86, but don't trust
it blindly. Check the modules against simprod to make sure this
does the same thing!
"""

from optparse import OptionParser
import os
import string

usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)
parser.add_option("-o", "--outfile", default=None,
                  dest="OUTFILE", help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_option("-i", "--infile", default="test_flashes_clsim.i3",
                  dest="INFILE", help="Read input from INFILE (.i3{.gz} format)")
parser.add_option("-r", "--runnumber", type="int", default=1,
                  dest="RUNNUMBER", help="The run number for this simulation")
parser.add_option("-s", "--seed",type="int",default=12345,
                  dest="SEED", help="Initial seed for the random number generator")
parser.add_option("--keep-mchits", action="store_true", default=False,
                  dest="KEEPMCHITS", help="Keep I3MCHits before writing the output file")
parser.add_option("--use-domlauncher", action="store_true", default=False,
                  dest="USEDOMLAUNCHER", help="Use DOMLauncher instead of PMTsimulator/DOMsimulator")

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
        print 'using input file %s' % infile
else:
        print "No input file!"
        parser.print_help()
        exit(-1)

infileRoot, infileExt = os.path.splitext(infile)
if infileExt == ".gz":
    infileRoot2, infileExt2 = os.path.splitext(infileRoot)
    if infileExt2 == ".i3":
        infileRoot=infileRoot2
        infileExt = ".i3.gz"

if infileExt != ".i3" and infileExt != ".i3.gz":
        raise Exception, "you have to specify either a .i3 or an .i3.gz file!"

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
        outfile = infileRootFile + "_detsim"
        outfile = outfile + infileExt
print "output dir is %s" % outdir
print "output file is %s" % outdir + outfile

########################


from I3Tray import *
from os.path import expandvars
import os
import sys

from icecube import icetray, dataclasses, dataio, phys_services
from icecube import noise_generator, trigger_sim
from icecube.sim_services import bad_dom_list_static

if options.USEDOMLAUNCHER:
    from icecube import DOMLauncher
else:
    from icecube import pmt_simulator, DOMsimulator

tray = I3Tray()

tray.AddService("I3SPRNGRandomServiceFactory","random",
    seed = options.SEED,
    nstreams = 10000,
    streamnum = options.RUNNUMBER)

tray.AddModule("I3Reader","reader",
               Filename=infile)

# copy DrivingTime from I3EventHeader (I3GlobalTriggerSim needs it)
def makeDrivingTime(frame):
    if "DrivingTime" in frame: return
    header = frame["I3EventHeader"]
    frame["DrivingTime"] = header.start_time
tray.AddModule(makeDrivingTime, "makeDrivingTime", Streams=[icetray.I3Frame.DAQ])

# stolen from http://x2100.icecube.wisc.edu/svn/projects/simprod-scripts/trunk/python/simulation/detector.py [ IC86() ],
# bad DOM list changed to IC86, parameters inlined (ScaleFactor=1.0, FilterMode=True)
#
# don't trust this blindly for IC86 detector simulation, always compare to simprod scripts!
#

tray.AddModule("I3NoiseGeneratorModule","noiseic",
    ScaleFactor = 1.0,
    InIce = True,
    IceTop = False,
    EndWindow = 10.*I3Units.microsecond,
    StartWindow = 10.*I3Units.microsecond,
    IndividualRates = True, 
    InputHitSeriesMapName = "MCHitSeriesMap",
    DOMstoExclude = bad_dom_list_static.IC86_static_bad_dom_list())

if options.USEDOMLAUNCHER:
    tray.AddModule("PMTResponseSimulator","rosencrantz",
        Input="MCHitSeriesMap",
        Output="weightedMCHitSeriesMap")
    tray.AddModule("DOMLauncher", "guildenstern",
        Input="weightedMCHitSeriesMap",
        Output="InIceRawData")
else:
    tray.AddModule("I3PMTSimulator","pmt")
    tray.AddModule("I3DOMsimulator","domsimulator")

## The usual SMT8
tray.AddModule("SimpleMajorityTrigger","IISMT8",
     TriggerConfigID = 1006)

####### The simprod script has this commented out:
# 
# ## The new SMT3 w/ the reduced DomSet
# tray.AddModule("SimpleMajorityTrigger","DCSMT3",
#      TriggerConfigID = 1011)

## The usual 5/7 string trigger a.k.a. ClusterTrigger
tray.AddModule("ClusterTrigger","string",
     TriggerConfigID = 1007)

tray.AddModule("I3GlobalTriggerSim","globaltrigger",
    FilterMode = True)

tray.AddModule("I3Pruner","pruner",
    GlobalTriggerName = "I3TriggerHierarchy",
    DOMLaunchSeriesMapNames = ["InIceRawData"])

if options.USEDOMLAUNCHER:
    MCPMTResponseMapNames = []
    MCHitSeriesMapNames = ["MCHitSeriesMap", "weightedMCHitSeriesMap"]
else:
    MCPMTResponseMapNames = ["MCPMTResponseMap"]
    MCHitSeriesMapNames = ["MCHitSeriesMap"]

tray.AddModule("I3TimeShifter","timeshifter",
    I3DOMLaunchSeriesMapNames = ["InIceRawData"],
    I3MCPMTResponseMapNames = MCPMTResponseMapNames,
    FlasherInfoName = "I3FlasherInfo", # the flasher info object needs to be time-shifted
    I3MCTreeNames = [],
    I3MCHitSeriesMapNames = MCHitSeriesMapNames,
    ShiftUntriggeredEvents = False)


# clean up
tray.AddModule("Delete", "cleanup",
    Keys = ["MCTimeIncEventID",
            "MCPMTResponseMap",
            "weightedMCHitSeriesMap"])

if not options.KEEPMCHITS:
    tray.AddModule("Delete", "cleanup_I3MCHits",
        Keys = ["MCHitSeriesMap"])


tray.AddModule("I3Writer","writer",
    Filename = outdir+outfile)

tray.AddModule("TrashCan", "the can")

tray.Execute()
tray.Finish()
