#!/usr/bin/env python
#
#
from I3Tray import *

from os.path import expandvars
from glob import glob

argc = len(sys.argv)
if argc < 3:
    print 'usage:', sys.argv[0], 'input_file output_file'
    sys.exit(1);

input_file = sys.argv[1]
output_file = sys.argv[2]


import os
import sys

load("libdataclasses")
load("libdataio")

tray = I3Tray()

tray.AddModule("I3Reader","reader")(
    ("FileNameList", [input_file])
    )

#tray.AddModule("Dump","dump")

tray.AddModule("I3Writer", "writer")(
    ("Filename", output_file),
    ("Streams", ["Physics"]),
    )
 
tray.AddModule("TrashCan", "the can");

tray.Execute()
tray.Finish()
