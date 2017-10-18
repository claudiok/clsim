#!/usr/bin/env python

from __future__ import print_function
from optparse import OptionParser

usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)
parser.add_option("-o","--outfile",dest="OUTFILE",help="Write output to OUTFILE",default=None)
parser.add_option("-f","--format",dest="format",help="format to output [hdf5, root, or csv]",default='hdf5')
parser.add_option("-z","--compress",dest="compression",help="compression level",default=1,type=int)
parser.add_option("-n","--nevents",dest="nevents",help="number of events to read",default=None,type=int)

# parse cmd line args, bail out if anything is not understood
(options,args) = parser.parse_args()
if len(args) != 1:
        crap = "Got undefined options:"
        for a in args:
                crap += a
                crap += " "
        parser.error(crap)

infile = args[0]

if options.OUTFILE is None:
    outfile = infile + '.' + options.format
else:
    outfile = options.OUTFILE


from I3Tray import *
from os.path import expandvars
import os
import sys

from icecube import icetray,dataclasses,dataio,tableio,clsim

if options.format == 'hdf5':
    try:
        from icecube import hdfwriter
    except ImportError:
        raise ImportError("Couldn't find the HDF writer service")
    tabler = hdfwriter.I3HDFTableService(outfile,options.compression)
elif options.format == 'root':
    try:
        from icecube import rootwriter
    except ImportError:
        raise ImportError("Couldn't find the ROOT writer service")
    tabler = rootwriter.I3ROOTTableService(outfile,options.compression)
elif options.format == 'csv':
    tabler = tableio.I3CSVTableService(outfile[:-4] + '_csv')
else:
    raise ValueError("I don't have a writer service for format '%s'"%options.format)

tray = I3Tray()

tray.AddModule('I3Reader','reader',filename=infile)

tray.AddModule('I3NullSplitter', 'nullsplit')

count = 0
def counter(frame):
    global count
    if (count%100==0):
        print("%d frames"%count)
    count +=1

def extractMCTreeParticles(frame):
    mctree = frame["I3MCTree"]
    
    # get tracks in ice/water
    mctracksinice = mctree.get_filter(lambda p: p.location_type==p.InIce)
    mcmuontrack = None
    highestEnergy=-1.
    for track in mctracksinice:
        if track.energy > highestEnergy:
            highestEnergy = track.energy
            mcmuontrack = track

    if mcmuontrack is None:
        mcmuontrack = dataclasses.I3Particle()
    primary = dataclasses.get_most_energetic_primary(mctree)
    
    if primary is None:
        primary = dataclasses.I3Particle()
    
    frame["MCMostEnergeticInIce"] = mcmuontrack
    frame["MCMostEnergeticPrimary"] = primary
tray.AddModule(extractMCTreeParticles, "extractMCTreeParticles")

tray.AddModule(tableio.I3TableWriter,'writer',
    tableservice = tabler,
    keys=[
          dict(key="MCMostEnergeticInIce", name="MCMostEnergeticInIce"),
          dict(key="MCMostEnergeticPrimary", name="MCMostEnergeticPrimary"),
          dict(key="MCPESeriesMap", converter=clsim.converters.I3MCPESeriesMapConverterWithIDs(True), name="MCPESeriesMap"),
         ],
    SubEventStreams=['nullsplit',])


tray.AddModule(counter,'countvoncount')


if options.nevents is None:
    tray.Execute()
else:
    tray.Execute(options.nevents+4)

