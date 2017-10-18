#!/usr/bin/env python

from optparse import OptionParser
import sys, os

usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)
parser.add_option("-n", "--eventnum", type="int", default=1,
                  dest="EVENTNUM", help="The number of events to display")
parser.add_option("-s", "--skipevents", type="int", default=0,
                  dest="SKIPNUM", help="Number of events to skip")

# parse cmd line args, bail out if anything is not understood
(options,args) = parser.parse_args()

if len(args) == 0:
    parser.error("No input filenames specified")

firstFileName = os.path.basename(args[0])
fileBase, fileExt = os.path.splitext(firstFileName)
outFile = fileBase + ".pdf"

###############################

from event_plotter import eventPlotter
from event_plotter2 import eventPlotter2

from I3Tray import *
from icecube import icetray, dataclasses, dataio, interfaces, phys_services
from icecube import clsim

tray = I3Tray()

tray.AddModule("I3Reader", "reader",
               FilenameList = args)

skippedEvents=0
pickedEvents=0
def counter(frame, SkipNevents=0, NEventsToPick=0):
    global skippedEvents
    global pickedEvents
    
    if skippedEvents < SkipNevents:
        skippedEvents+=1
        return False
    
    pickedEvents+=1
    if NEventsToPick>0 and pickedEvents>NEventsToPick:
        return False
    return True

tray.AddModule(counter, "counter", Streams=[icetray.I3Frame.DAQ],
               SkipNevents=options.SKIPNUM,
               NEventsToPick=options.EVENTNUM)

tray.AddModule(eventPlotter2, "plotModule",
               filename=outFile,
               pi0Only=False,
               eplusOnly=False)



tray.Execute()

