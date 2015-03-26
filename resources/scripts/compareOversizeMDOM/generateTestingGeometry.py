#!/usr/bin/env python

from optparse import OptionParser

usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-x", "--xpos", type="float", default=0.,
                  dest="XPOS", help="The x coordinate in meters")
parser.add_option("-y", "--ypos", type="float", default=0.,
                  dest="YPOS", help="The y coordinate in meters")
parser.add_option("-z", "--zpos", type="float", default=0.,
                  dest="ZPOS", help="The z coordinate in meters")

# parse cmd line args, bail out if anything is not understood
(options,args) = parser.parse_args()
if len(args) != 0:
        crap = "Got undefined options:"
        for a in args:
                crap += a
                crap += " "
        parser.error(crap)


outfile = "testingGCD.i3"

import numpy

from I3Tray import *
import os
import sys
import copy

from icecube import icetray, dataclasses, dataio


radius = 40.*I3Units.m

omPos = numpy.array(
    [[ 0.,  1., 0.],
     [ 1.,  1., 0.],
     [ 1.,  0., 0.],
     [ 1., -1., 0.],
     [ 0., -1., 0.],
     [-1., -1., 0.],
     [-1.,  0., 0.],
     [-1.,  1., 0.]]
    )
# normalize and scale
omPos = (omPos.T/numpy.sqrt(numpy.sum(omPos**2, 1))).T * radius

omPosLower = numpy.array(omPos)
omPosLower.T[2] = omPosLower.T[2] - radius
omPosUpper = numpy.array(omPos)
omPosUpper.T[2] = omPosUpper.T[2] + radius

omPositions = numpy.concatenate((omPosUpper, omPos, omPosLower), axis=0)


omKeys = [
    icetray.OMKey(1,1),
    icetray.OMKey(2,1),
    icetray.OMKey(3,1),
    icetray.OMKey(4,1),
    icetray.OMKey(5,1),
    icetray.OMKey(6,1),
    icetray.OMKey(7,1),
    icetray.OMKey(8,1),

    icetray.OMKey(1,2),
    icetray.OMKey(2,2),
    icetray.OMKey(3,2),
    icetray.OMKey(4,2),
    icetray.OMKey(5,2),
    icetray.OMKey(6,2),
    icetray.OMKey(7,2),
    icetray.OMKey(8,2),

    icetray.OMKey(1,3),
    icetray.OMKey(2,3),
    icetray.OMKey(3,3),
    icetray.OMKey(4,3),
    icetray.OMKey(5,3),
    icetray.OMKey(6,3),
    icetray.OMKey(7,3),
    icetray.OMKey(8,3),

]

def get_orientations():
    for i in range(4):
        d = dataclasses.I3Direction(0,0,1)
        d.rotate_y((90-57)*icetray.I3Units.degree)
        d.rotate_z((22.5 + 90*i)*icetray.I3Units.degree)
        yield d
    for i in range(4,12):
        d = dataclasses.I3Direction(0,0,1)
        d.rotate_y((90-16)*icetray.I3Units.degree)
        d.rotate_z((0 + 45*i)*icetray.I3Units.degree)
        yield d
    for i in range(12,20):
        d = dataclasses.I3Direction(0,0,1)
        d.rotate_y((90+16)*icetray.I3Units.degree)
        d.rotate_z((22.5 + 45*i)*icetray.I3Units.degree)
        yield d
    for i in range(20,24):
        d = dataclasses.I3Direction(0,0,1)
        d.rotate_y((90+57)*icetray.I3Units.degree)
        d.rotate_z((-22.5 + 90*i)*icetray.I3Units.degree)
        yield d

def generate_multipmts(module_geo):
	# inner radius of the 14" mDOM sphere
	r_housing = (356*icetray.I3Units.mm/2. - 14*icetray.I3Units.mm)
	# radius of 3" PMT window
	r_pmt = 3*2.54*icetray.I3Units.cm/2.
	# distance from mDOM center to center of PMT window (assuming that the
	# outside edge of the window lies directly on inner surface of the
	# pressure housing)
	r_window = numpy.sqrt(r_housing**2 - r_pmt**2) - 5*icetray.I3Units.mm
	
	for d in get_orientations():
		omgeo = copy.copy(module_geo)
		omgeo.position += d*r_window
		omgeo.orientation = dataclasses.I3Orientation(d)
		omgeo.area = dataclasses.I3Constants.pi*r_pmt**2
		# omgeo.omtype = omgeo.OMType.UnknownType
		yield omgeo

geometry = dataclasses.I3Geometry()
calibration = dataclasses.I3Calibration()
detectorStatus = dataclasses.I3DetectorStatus()

# fill the geometry map
omgeomap = geometry.omgeo
modgeomap = dataclasses.I3ModuleGeoMap()
subdetectors = dataclasses.I3MapModuleKeyString()
domcalmap = calibration.dom_cal
domstatusmap = detectorStatus.dom_status

for i, pos in enumerate(omPositions):
    shiftedPos = pos
    shiftedPos[0] += options.XPOS*I3Units.m
    shiftedPos[1] += options.YPOS*I3Units.m
    shiftedPos[2] += options.ZPOS*I3Units.m

    omkey = omKeys[i]

    newomgeo = dataclasses.I3OMGeo()
    newomgeo.omtype = dataclasses.I3OMGeo.OMType.IceCube
    newomgeo.orientation = dataclasses.I3Orientation(dataclasses.I3Direction(0.,0.,-1.))
    newomgeo.position = dataclasses.I3Position(shiftedPos[0], shiftedPos[1], shiftedPos[2])
    
    for k, pmtgeo in enumerate(generate_multipmts(newomgeo)):
        pmtkey = icetray.OMKey(omkey.string, omkey.om, k)
        omgeomap[pmtkey] = pmtgeo
    
    mkey = dataclasses.ModuleKey(omkey.string, omkey.om)
    modgeo = dataclasses.I3ModuleGeo()
    modgeo.module_type = modgeo.mDOM
    modgeo.pos = newomgeo.position
    modgeo.radius = 356*I3Units.mm/2.
    modgeomap[mkey] = modgeo
    subdetectors[mkey] = "Gen2"

    newdomcal = dataclasses.I3DOMCalibration()
    newdomcal.relative_dom_eff = 1.0
    domcalmap[omkey] = newdomcal


    newdomstatus = dataclasses.I3DOMStatus()
    newdomstatus.pmt_hv = 1345.*I3Units.V # some arbitrary setting: >0 and not NaN
    domstatusmap[omkey] = newdomstatus



Gframe = icetray.I3Frame(icetray.I3Frame.Geometry)
Cframe = icetray.I3Frame(icetray.I3Frame.Calibration)
Dframe = icetray.I3Frame(icetray.I3Frame.DetectorStatus)

Gframe["I3Geometry"] = geometry
Gframe["I3ModuleGeoMap"] = modgeomap
Gframe["I3OMGeoMap"] = omgeomap
Gframe["Subdetectors"] = subdetectors
Cframe["I3Calibration"] = calibration
Dframe["I3DetectorStatus"] = detectorStatus

i3file = dataio.I3File(outfile, 'w')
i3file.push(Gframe)
i3file.push(Cframe)
i3file.push(Dframe)
i3file.close()



