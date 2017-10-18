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


import os
import sys

from icecube import icetray, dataclasses, dataio

tray = I3Tray()

tray.AddModule("I3Reader","reader",
    FileNameList = [input_file])

tray.AddModule("I3Writer", "writer",
    Filename = output_file,
    Streams = [icetray.I3Frame.Physics, icetray.I3Frame.DAQ])
 


tray.Execute()

