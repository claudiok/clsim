#!/usr/bin/env python

from __future__ import print_function
from optparse import OptionParser

usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)
parser.add_option("-f","--format",dest="format",help="format to output [hdf5, root, or csv]",default='hdf5')
parser.add_option("-z","--compress",dest="compression",help="compression level",default=1,type=int)
parser.add_option("-n","--frames",dest="nframes",help="number of frames to process",default=None,type=int)

import os

options,args = parser.parse_args()

if len(args) != 1:
    parser.error("You must supply an input filename")
    
infile = args[0]
outfile = os.path.basename(args[0]) + '.' + options.format

from icecube import icetray,dataclasses,dataio,tableio,clsim
from I3Tray import I3Tray

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

# tray.AddModule(tableio.I3TableWriter,'writer1',
#     tableservice = tabler, keys=['LineFit', 'InIceRawData'], types=[dataclasses.I3Particle])

def extractMCTreeParticles(frame):
    mctree = frame["I3MCTree"]
    
    # get tracks in ice/water
    mctracksinice = mctree.in_ice
    mcmuontrack = None
    highestEnergy=-1.
    for track in mctracksinice:
        if (track.type == dataclasses.I3Particle.ParticleType.MuPlus) or (track.type == dataclasses.I3Particle.ParticleType.MuMinus):
            if track.energy > highestEnergy:
                highestEnergy = track.energy
                mcmuontrack = track

    if mcmuontrack is None:
        mcmuontrack = dataclasses.I3Particle()
    primary = mctree.most_energetic_primary
    
    if primary is None:
        primary = dataclasses.I3Particle()
    
    frame["MCMostEnergeticMuon"] = mcmuontrack
    frame["MCMostEnergeticPrimary"] = primary
tray.AddModule(extractMCTreeParticles, "extractMCTreeParticles")

tray.AddModule(tableio.I3TableWriter,'writer1',
    tableservice = tabler,
    keys=[
          dict(key="MCMostEnergeticMuon", name="MCMostEnergeticMuon"),
          dict(key="MCMostEnergeticPrimary", name="MCMostEnergeticPrimary"),
          dict(key="MCHitSeriesMap", converter=clsim.converters.I3MCHitSeriesMapConverterWithIDs(True), name="MCHitSeriesMap"),
          dict(key="MCPESeriesMap_clsim", converter=clsim.converters.I3MCPESeriesMapConverterWithIDs(True), name="MCPESeriesMap_clsim"),
          dict(key="PropagatedPhotons", name="PropagatedPhotons"),
         ],
    SubEventStreams=['nullsplit',])

#tray.AddModule(tableio.I3TableWriter,'writer2',
#    tableservice = tabler, types=[dataclasses.I3Particle])



tray.AddModule(counter,'count-count')


if options.nframes is not None:
    tray.Execute(options.nframes)
else:
    tray.Execute()

