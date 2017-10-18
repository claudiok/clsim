#!/usr/bin/env python
#
#
from __future__ import print_function
from I3Tray import *

from os.path import expandvars
from glob import glob

import os

argc = len(sys.argv)
if argc < 3:
    print('usage:', sys.argv[0], 'physics_file_base output_file')
    sys.exit(1);

physics_file_base = sys.argv[1];
output_file = sys.argv[2];

#  use 'glob' to convert the string with the 'star' in it to a list of real filenames
file_list = glob(physics_file_base + "*.i3.gz")

#  puts them in alpha/numeric order
file_list.sort()

non_zero_file_list = [x for x in file_list if os.path.getsize(x)!=0]

print("working with") 
print(non_zero_file_list)

import os
import sys

load("libdataclasses")
load("libdataio")

tray = I3Tray()

tray.AddModule("I3Reader","reader")(
    ("FileNameList", non_zero_file_list)
    )

#tray.AddModule("Dump","dump")

tray.AddModule("I3Writer", "writer")(
    ("Filename", output_file),
    ("Streams", ["Physics"])
    )
 


tray.Execute()

