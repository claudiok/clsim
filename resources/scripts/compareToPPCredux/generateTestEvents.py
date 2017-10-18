#!/usr/bin/env python

from optparse import OptionParser
from os.path import expandvars

usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)
parser.add_option("-o", "--outfile",default="test_events.i3",
                  dest="OUTFILE", help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_option("-s", "--seed",type="int",default=12345,
                  dest="SEED", help="Initial seed for the random number generator")
parser.add_option("-g", "--gcd",default="testingGCD.i3",
                  dest="GCDFILE", help="Read geometry from GCDFILE (.i3{.gz} format)")
parser.add_option("-r", "--runnumber", type="int", default=1,
                  dest="RUNNUMBER", help="The run number for this simulation")
parser.add_option("-n", "--numevents", type="int", default=100,
                  dest="NUMEVENTS", help="The number of events per run")

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

from I3Tray import *
import os
import sys

from icecube import icetray, dataclasses, dataio, phys_services

xCoord = options.XPOS*I3Units.m
yCoord = options.YPOS*I3Units.m
zCoord = options.ZPOS*I3Units.m

load("libsim-services")

import math
import numpy

class generateEvent(icetray.I3Module):
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddParameter("I3RandomService", "the service", None)
        self.AddParameter("Type", "", dataclasses.I3Particle.ParticleType.EMinus)
        self.AddParameter("Energy", "", 100*I3Units.TeV)
        self.AddParameter("NEvents", "", 100)
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


tray = I3Tray()

# a random number generator
randomService = phys_services.I3GSLRandomService(options.SEED)

tray.AddModule("I3InfiniteSource","streams",
               Prefix=options.GCDFILE,
               Stream=icetray.I3Frame.DAQ)

tray.AddModule("I3MCEventHeaderGenerator","gen_header",
               Year=2009,
               DAQTime=158100000000000000,
               RunNumber=1,
               EventID=1,
               IncrementEventID=True)

tray.AddModule(generateEvent, "generateEvent",
               I3RandomService = randomService,
               NEvents = options.NUMEVENTS,
               Energy = 1000.*I3Units.TeV,
               XCoord = xCoord,
               YCoord = yCoord,
               ZCoord = zCoord,
               )

tray.AddModule("I3Writer","writer",
    Filename = options.OUTFILE)



tray.Execute()




