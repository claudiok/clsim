#!/usr/bin/env python

from optparse import OptionParser
from os.path import expandvars

usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)
parser.add_option("-n", "--numevents", type="int", default=1,
                  dest="NUMEVENTS", help="The number of events per run")
parser.add_option("-s", "--seed",type="int",default=12345,
                  dest="SEED", help="Initial seed for the random number generator")
parser.add_option("-r", "--runnumber", type="int", default=1,
                  dest="RUNNUMBER", help="The run number for this simulation")
parser.add_option("-x", "--xmlfile", default="benchmark.xml",
                  dest="XMLFILE", help="Write statistics to XMLFILE")
parser.add_option("-p", "--max-parallel-events", type="int", default=100,
                  dest="MAXPARALLELEVENTS", help="maximum number of events(==frames) that will be processed in parallel")
parser.add_option("--icemodel", default=expandvars("$I3_SRC/clsim/resources/ice/spice_lea"),
                  dest="ICEMODEL", help="A clsim ice model file/directory (ice models *will* affect performance metrics, always compare using the same model!)")
parser.add_option("--use-cpu",  action="store_true", default=False,
                  dest="USECPU", help="simulate using CPU instead of GPU")

parser.add_option("-d", "--device", type="int", default=None,
                  dest="DEVICE", help="device number")

# parse cmd line args, bail out if anything is not understood
(options,args) = parser.parse_args()
if len(args) != 0:
        crap = "Got undefined options:"
        for a in args:
                crap += a
                crap += " "
        parser.error(crap)

if options.DEVICE is not None:
    print " ** DEVICE selected using the \"--device\" command line option. Only do this if you know what you are doing!"
    print " ** You should be using the CUDA_VISIBLE_DEVICES and/or GPU_DEVICE_ORDINAL environment variables instead."

from I3Tray import *
import os
import sys
import math
import numpy

from icecube import icetray, dataclasses, dataio, phys_services, sim_services, clsim

# icetray.I3Logger.global_logger.set_level(icetray.I3LogLevel.LOG_INFO)
icetray.I3Logger.global_logger.set_level(icetray.I3LogLevel.LOG_WARN)

radius = 120.*I3Units.m

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


class generateEvent(icetray.I3Module):
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddParameter("I3RandomService", "the service", None)
        self.AddParameter("Type", "", dataclasses.I3Particle.ParticleType.EMinus)
        self.AddParameter("Energy", "", 10.*I3Units.TeV)
        self.AddParameter("NEvents", "", 1)
        self.AddParameter("XCoord", "", 0.)
        self.AddParameter("YCoord", "", 0.)
        self.AddParameter("ZCoord", "", 0.)

        self.AddOutBox("OutBox")        

    def Configure(self):
        self.rs = self.GetParameter("I3RandomService")
        self.particleType = self.GetParameter("Type")
        self.energy = self.GetParameter("Energy")
        self.nEvents = self.GetParameter("NEvents")
        self.xCoord = self.GetParameter("XCoord")
        self.yCoord = self.GetParameter("YCoord")
        self.zCoord = self.GetParameter("ZCoord")
        
        self.eventCounter = 0

    def DAQ(self, frame):
        daughter = dataclasses.I3Particle()
        daughter.type = self.particleType
        daughter.energy = self.energy
        daughter.pos = dataclasses.I3Position(self.xCoord,self.yCoord,self.zCoord)
        daughter.dir = dataclasses.I3Direction(0.,0.,-1.)
        daughter.time = 0.
        daughter.location_type = dataclasses.I3Particle.LocationType.InIce

        primary = dataclasses.I3Particle()
        primary.type = dataclasses.I3Particle.ParticleType.NuE
        primary.energy = self.energy
        primary.pos = dataclasses.I3Position(self.xCoord,self.yCoord,self.zCoord)
        primary.dir = dataclasses.I3Direction(0.,0.,-1.)
        primary.time = 0.
        primary.location_type = dataclasses.I3Particle.LocationType.Anywhere

        mctree = dataclasses.I3MCTree()
        mctree.add_primary(primary)
        mctree.append_child(primary,daughter)

        frame["I3MCTree"] = mctree

        self.PushFrame(frame)

        self.eventCounter += 1
        if self.eventCounter==self.nEvents:
            self.RequestSuspension()


class injectFakeGCD(icetray.I3Module):
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddParameter("OMKeys",      "", [])
        self.AddParameter("OMPositions", "", [])
        self.AddParameter("XCoord", "", 0.)
        self.AddParameter("YCoord", "", 0.)
        self.AddParameter("ZCoord", "", 0.)

        self.AddOutBox("OutBox")        

    def Configure(self):
        self.omkeys = self.GetParameter("OMKeys")
        self.ompositions = self.GetParameter("OMPositions")
        self.xCoord = self.GetParameter("XCoord")
        self.yCoord = self.GetParameter("YCoord")
        self.zCoord = self.GetParameter("ZCoord")

        self.has_been_injected = False

    def DAQ(self, frame):
        # only inject it once
        if self.has_been_injected:
            self.PushFrame(frame)
            return
        self.has_been_injected = True

        geometry = dataclasses.I3Geometry()
        calibration = dataclasses.I3Calibration()
        detectorStatus = dataclasses.I3DetectorStatus()

        # fill the geometry map
        omgeomap = geometry.omgeo
        domcalmap = calibration.dom_cal
        domstatusmap = detectorStatus.dom_status

        for i, pos in enumerate(omPositions):
            shiftedPos = pos
            shiftedPos[0] += self.xCoord*I3Units.m
            shiftedPos[1] += self.yCoord*I3Units.m
            shiftedPos[2] += self.zCoord*I3Units.m

            omkey = omKeys[i]

            newomgeo = dataclasses.I3OMGeo()
            newomgeo.omtype = dataclasses.I3OMGeo.OMType.IceCube
            newomgeo.orientation = dataclasses.I3Orientation(dataclasses.I3Direction(0.,0.,-1.))
            newomgeo.position = dataclasses.I3Position(shiftedPos[0], shiftedPos[1], shiftedPos[2])
            omgeomap[omkey] = newomgeo


            newdomcal = dataclasses.I3DOMCalibration()
            newdomcal.relative_dom_eff = 1.0
            domcalmap[omkey] = newdomcal


            newdomstatus = dataclasses.I3DOMStatus()
            newdomstatus.pmt_hv = 1345.*I3Units.V # some arbitrary setting: >0 and not NaN
            domstatusmap[omkey] = newdomstatus


        # make GCD frames and fill them with objects
        Gframe = icetray.I3Frame(icetray.I3Frame.Geometry)
        Cframe = icetray.I3Frame(icetray.I3Frame.Calibration)
        Dframe = icetray.I3Frame(icetray.I3Frame.DetectorStatus)

        Gframe["I3Geometry"] = geometry
        Cframe["I3Calibration"] = calibration
        Dframe["I3DetectorStatus"] = detectorStatus

        # push the new GCD frames
        self.PushFrame(Gframe)
        self.PushFrame(Cframe)
        self.PushFrame(Dframe)

        # push the original Q-frame
        self.PushFrame(frame)


tray = I3Tray()

tray.AddService("I3XMLSummaryServiceFactory","summary",
    OutputFileName = options.XMLFILE)

# a random number generator
randomService = phys_services.I3SPRNGRandomService(
    seed = options.SEED,
    nstreams = 10000,
    streamnum = options.RUNNUMBER)

# use a real GCD file for a real-world test
tray.AddModule("I3InfiniteSource","streams",
    Prefix = expandvars("$I3_PORTS/test-data/sim/GeoCalibDetectorStatus_IC86.55697_corrected_V2.i3.gz"),
    Stream=icetray.I3Frame.DAQ)

## this would be an alternative to a real GCD file. comment out "Prefix" above and uncomment this
# tray.AddModule(injectFakeGCD,"gcd",
#     OMKeys = omKeys,
#     OMPositions = omPositions,
#     # XCoord = xCoord,
#     # YCoord = yCoord,
#     # ZCoord = zCoord,
#     )

tray.AddModule("I3MCEventHeaderGenerator","gen_header",
    Year=2009,
    DAQTime=158100000000000000,
    RunNumber=1,
    EventID=1,
    IncrementEventID=True)

tray.AddModule(generateEvent, "generateEvent",
    I3RandomService = randomService,
    NEvents = options.NUMEVENTS,
    Energy = 40.*I3Units.TeV,
    # Energy = 1000.*I3Units.TeV,
    # XCoord = xCoord,
    # YCoord = yCoord,
    # ZCoord = zCoord,
    )

MCTreeName="I3MCTree"
MMCTrackListName=None
photonSeriesName = None

tray.AddSegment(clsim.I3CLSimMakeHits, "makeCLSimHits",
    PhotonSeriesName = photonSeriesName,
    MCTreeName = MCTreeName,
    MMCTrackListName = MMCTrackListName,
    ParallelEvents = options.MAXPARALLELEVENTS,
    RandomService = randomService,
    MCPESeriesName = "MCPESeriesMap",
    UnshadowedFraction = 0.95,
    UseGPUs=not options.USECPU,
    UseCPUs=options.USECPU,
    UseOnlyDeviceNumber=options.DEVICE,
    IceModelLocation=options.ICEMODEL,
    ExtraArgumentsToI3CLSimModule={"EnableDoubleBuffering":True}
    )

tray.AddModule("TrashCan", "the can")

tray.Execute()
tray.Finish()

del tray

########### this is optional and just parses the generated XML file

# retrieve the XML file and parse it

import xml.etree.ElementTree as ET
tree = ET.parse(options.XMLFILE)
root = tree.getroot()
ns_per_photon = [float(item.find('second').text) for item in root.find('I3XMLSummaryService').find('map').findall('item') if item.find('first').text=="I3CLSimModule_makeCLSimHits_makePhotons_clsim_AverageDeviceTimePerPhoton"][0]
ns_per_photon_with_util = [float(item.find('second').text) for item in root.find('I3XMLSummaryService').find('map').findall('item') if item.find('first').text=="I3CLSimModule_makeCLSimHits_makePhotons_clsim_AverageHostTimePerPhoton"][0]
device_util = [float(item.find('second').text) for item in root.find('I3XMLSummaryService').find('map').findall('item') if item.find('first').text=="I3CLSimModule_makeCLSimHits_makePhotons_clsim_DeviceUtilization"][0]

print " "
print "# these numbers are performance figures for the GPU:"
print "time per photon (GPU):", ns_per_photon, "ns"
print "photons per second (GPU):", 1e9/ns_per_photon, "photons per second"

print " "
print "# these numbers include the host utilization and are probably not meaningful for --numevents=1 (the default). You need more events to even out the startup/setup time."
print "time per photon (actual, including under-utilization):", ns_per_photon_with_util, "ns"
print "photons per second (actual, including under-utilization):", 1e9/ns_per_photon_with_util, "photons per second"

print "device utilization:", device_util*100., "%"

