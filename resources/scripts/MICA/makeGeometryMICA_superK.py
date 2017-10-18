#!/usr/bin/env python
from __future__ import print_function
from I3Tray import *

from icecube import icetray, dataclasses, dataio, sim_services, phys_services

from icecube.geometry_maker.detectors import *
from icecube.geometry_maker.strings import SimpleString
from icecube.geometry_maker.floors import SimpleFloor
from icecube.geometry_maker.oms import mDOM, ANTARESOM
from icecube.geometry_maker.pmts import MultiPMT

from os.path import expandvars
import math, numpy

tray = I3Tray()

tray.AddModule("I3InfiniteSource","streams",
               Stream=icetray.I3Frame.DAQ,
               Prefix=expandvars("$I3_TESTDATA/sim/GeoCalibDetectorStatus_IC86.55380_corrected.i3.gz"))

tray.AddModule("I3MCEventHeaderGenerator","gen_header",
               Year=2009,
               DAQTime=158100000000000000,
               RunNumber=1,
               EventID=1,
               IncrementEventID=True)

tray.AddModule("I3GeometryDecomposer", "decompose")

ringRadius = numpy.array([57.*I3Units.m, 50.*I3Units.m])
numberOfStrings = 126
stringDistance = 2.*math.pi*ringRadius/float(numberOfStrings/len(ringRadius))
domsPerString=176
domSpacing=1.*I3Units.m
stringLength=float(domsPerString-1)*domSpacing
totalNumDOMs=domsPerString*numberOfStrings
depthAtZ0 = 1948.07*I3Units.m # depth @ IceCube z==0
bedrockDepth = 2810.*I3Units.m
bottomDOMDepth = depthAtZ0+503.*I3Units.m

print("                  radius: {0}m".format(ringRadius/I3Units.m))
print("         string distance: {0}m".format(stringDistance/I3Units.m))
print("       number of strings: {0}".format(numberOfStrings))
print("         DOMs per string: {0}".format(domsPerString))
print("    DOM vertical spacing: {0}m".format(domSpacing/I3Units.m))
print("instrumented string len.: {0}m".format(stringLength/I3Units.m))
print("    total number of DOMs: {0}".format(totalNumDOMs))

myPMT = MultiPMT()
myOM = mDOM(pmtDescription=myPMT)
myFloor = SimpleFloor(omDescription=myOM)
myString = SimpleString(numberOfFloors=domsPerString, 
                        instrumentedStringLength=stringLength, 
                        distToFirstFloor=0., # distance up from the bedrock to the bottommost DOM
                        floorDescription=myFloor)
myDetector = RingDetector(numberOfStrings=numberOfStrings,
                          ringRadius=ringRadius,
                          stringDescription=myString)

tray.AddModule("I3GeometryMaker","geo_maker",
    SubdetectorName = "MICA",
    DetectorService = myDetector,
    BedrockDepth = bedrockDepth,
    DepthAtZ0 = depthAtZ0,
    CenterDetector = False, # do not re-center the geometry
    DetectorShiftPosX = 0.*I3Units.m,
    DetectorShiftPosY = 0.*I3Units.m,
    DetectorShiftPosZ = depthAtZ0-bottomDOMDepth)

tray.AddModule("I3Writer", "writer",
    filename = "geometry_MICA.i3",
    Streams = [icetray.I3Frame.Geometry, icetray.I3Frame.Calibration, icetray.I3Frame.DetectorStatus])

tray.AddModule("Dump", "dump")

tray.AddModule("FrameCheck", "check",
    ensure_physics_has = ["I3Calibration", "I3Geometry", 
                          "I3DetectorStatus"])



tray.Execute(5)

