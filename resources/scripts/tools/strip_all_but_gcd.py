#!/usr/bin/env python
#
#
from __future__ import print_function
from I3Tray import *

from os.path import expandvars
from glob import glob

argc = len(sys.argv)
if argc < 3:
    print('usage:', sys.argv[0], 'input_file output_file')
    sys.exit(1);

input_file = sys.argv[1]
output_file = sys.argv[2]

from icecube.icetray import I3Frame

import os
import sys

load("libdataclasses")
load("libdataio")

tray = I3Tray()

tray.AddModule("I3Reader","reader")(
    ("FileNameList", [input_file]),
    #("DontMergeStreamTypes", [I3Frame.Stream('!'), I3Frame.Stream('G'), I3Frame.Stream('C'), I3Frame.Stream('D')]),
    ("SkipKeys", ["I3MCTreeBeforeG4Sim"])
    )

tray.AddModule("Delete", "del",
    Keys = ["I3MCTreeBeforeG4Sim"])

#tray.AddModule("Dump","dump")

tray.AddModule("I3Writer", "writer")(
    ("Filename", output_file),
    ("Streams", ["Geometry", "Calibration", "DetectorStatus"]),
    )
 


tray.Execute(10)

