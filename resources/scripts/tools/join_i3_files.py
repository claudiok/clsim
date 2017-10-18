#!/usr/bin/env python
#
#
from __future__ import print_function
from I3Tray import *

from os.path import expandvars
from glob import glob

argc = len(sys.argv)
if argc < 3:
    print('usage:', sys.argv[0], 'output_file input_files')
    sys.exit(1);

output_file = sys.argv[1];

file_list = sys.argv;

# remove the first two entries from the file list (the executable name and the
# output file name)
file_list.pop(0)
file_list.pop(0)

print("==============")
print("Merging files: ")
print(file_list)
print("==============")

import os
import sys

load("libdataclasses")
load("libdataio")

tray = I3Tray()

tray.AddModule("I3Reader","reader")(
    ("FileNameList", file_list)
    )

#tray.AddModule("Dump","dump")

tray.AddModule("I3Writer", "writer")(
    ("Filename", output_file),
    ("Streams", ["Physics", "Geometry", "Calibration", "DetectorStatus"]),
    ("SkipKeys", ["I3SkipNEventFilter"])
    )
 


tray.Execute()

