#!/usr/bin/env python
"""
This is an EXAMPLE script of how to apply detector simulation.
There does not seem to be a collection of simulation tray segments
anywhere, so this uses a bunch of modules copied from simprod
scripts. This simulation should be valid for IC86, but don't trust
it blindly. Check the modules against simprod to make sure this
does the same thing!
"""

from __future__ import print_function
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
parser.add_option("--keep-mcpes", action="store_true", default=False,
                  dest="KEEPMCPES", help="Keep I3MCPEs before writing the output file")

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
        outfile = infileRootFile + "_detsim"
        outfile = outfile + infileExt
print("output dir is %s" % outdir)
print("output file is %s" % outdir + outfile)

########################


from I3Tray import *
from os.path import expandvars
import os
import sys

from icecube import icetray, dataclasses, dataio, phys_services, trigger_sim, vuvuzela
from icecube import DOMLauncher

tray = I3Tray()

tray.Add("I3SPRNGRandomServiceFactory",
    seed = options.SEED,
    nstreams = 10000,
    streamnum = options.RUNNUMBER)

tray.Add("I3Reader", Filename = infile)
                
# Add noise parameters
noise_file = expandvars("$I3_SRC/vuvuzela/resources/data/parameters.dat")
tray.Add("Inject", InputNoiseFile = noise_file)

tray.Add("Vuvuzela", 
	 InputHitSeriesMapName  = "",
	 OutputHitSeriesMapName = "I3MCPESeriesMap",
	 StartWindow            = 0,
	 EndWindow              = 25*I3Units.millisecond,
	 IceTop                 = False,
	 InIce                  = True,
	 RandomServiceName      = "I3RandomService",
	 UseIndividual          = True,
	 DisableLowDTCutoff     = True,
)

tray.Add("PMTResponseSimulator")
tray.Add("DOMLauncher")
tray.Add(trigger_sim.TriggerSim)
    
if not options.KEEPMCPES:
    tray.Add("Delete", Keys = ["I3MCPESeriesMap"])        

tray.Add("I3Writer", Filename = outdir+outfile)    

tray.Execute()

