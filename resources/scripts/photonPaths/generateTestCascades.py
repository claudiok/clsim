#!/usr/bin/env python

from optparse import OptionParser
from os.path import expandvars

usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)
parser.add_option("-o", "--outfile",default="test_cascades.i3",
                  dest="OUTFILE", help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_option("-s", "--seed",type="int",default=12345,
                  dest="SEED", help="Initial seed for the random number generator")
parser.add_option("-g", "--gcd",default=expandvars("$I3_TESTDATA/sim/GeoCalibDetectorStatus_IC86.55380_corrected.i3.gz"),
                  dest="GCDFILE", help="Read geometry from GCDFILE (.i3{.gz} format)")
parser.add_option("-r", "--runnumber", type="int", default=1,
                  dest="RUNNUMBER", help="The run number for this simulation")
parser.add_option("-n", "--numevents", type="int", default=1,
                  dest="NUMEVENTS", help="The number of events per run")
parser.add_option("--no-apply-mmc", action="store_false", default=True,
                  dest="APPLYMMC", help="do not apply MMC to the I3MCTree after generating the muons")

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

from icecube import icetray, dataclasses, dataio, phys_services, sim_services

from PropagateMuons import PropagateMuons

#
import math
import numpy

class mySimpleCascade(icetray.I3Module):
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddParameter("I3RandomService", "the service", None)
        self.AddParameter("Type", "", dataclasses.I3Particle.ParticleType.MuMinus)
        self.AddParameter("EnergyMin", "", 10*I3Units.TeV)
        self.AddParameter("EnergyMax", "", 10*I3Units.TeV)
        self.AddParameter("ZenithMin", "", 0*I3Units.degree)
        self.AddParameter("ZenithMax", "", 180*I3Units.degree)
        self.AddParameter("AzimuthMin", "", 0*I3Units.degree)
        self.AddParameter("AzimuthMax", "", 360*I3Units.degree)
        self.AddParameter("CylinderRadius", "", 500 * I3Units.m)
        self.AddParameter("CylinderHeight", "", 1000 * I3Units.m)
        self.AddParameter("NEvents", "", 100)

        self.AddOutBox("OutBox")        

    def Configure(self):
        self.rs = self.GetParameter("I3RandomService")
        self.particleType = self.GetParameter("Type")
        self.energyMin = self.GetParameter("EnergyMin")
        self.energyMax = self.GetParameter("EnergyMax")
        self.zenithMin = self.GetParameter("ZenithMin")
        self.zenithMax = self.GetParameter("ZenithMax")
        self.azimuthMin = self.GetParameter("AzimuthMin")
        self.azimuthMax = self.GetParameter("AzimuthMax")
        self.cylinderRadius = self.GetParameter("CylinderRadius")
        self.cylinderHeight = self.GetParameter("CylinderHeight")
        self.nEvents = self.GetParameter("NEvents")

    def DAQ(self, frame):

        azi = self.rs.uniform(self.azimuthMin,self.azimuthMax)

        cos_zen_low = math.cos(self.zenithMin / I3Units.radian)
        cos_zen_high = math.cos(self.zenithMax / I3Units.radian )
        zen = math.acos(self.rs.uniform(cos_zen_low,cos_zen_high))

        while True:
          posX = self.rs.uniform(-self.cylinderRadius,self.cylinderRadius)
          posY = self.rs.uniform(-self.cylinderRadius,self.cylinderRadius)
          if posX**2+posY**2 < self.cylinderRadius**2: break
        posZ = self.rs.uniform(-self.cylinderHeight/2.,self.cylinderHeight/2.)

        # set the particle's energy
        energy = self.rs.uniform(self.energyMin,self.energyMax) * I3Units.GeV        

        daughter = dataclasses.I3Particle()
        daughter.type = self.particleType
        daughter.energy = energy
        daughter.pos = dataclasses.I3Position(posX, posY, posZ)
        daughter.dir = dataclasses.I3Direction(zen,azi)
        daughter.time = 0.
        daughter.location_type = dataclasses.I3Particle.LocationType.InIce

        primary = dataclasses.I3Particle()
        primary.type = dataclasses.I3Particle.ParticleType.NuMu
        primary.energy = energy
        primary.pos = dataclasses.I3Position(posX, posY, posZ)
        primary.dir = dataclasses.I3Direction(zen,azi)
        primary.time = 0.
        primary.location_type = dataclasses.I3Particle.LocationType.Anywhere

        mctree = dataclasses.I3MCTree()
        mctree.add_primary(primary)
        mctree.append_child(primary,daughter)

        frame["I3MCTree_preMuonProp"] = mctree

        self.PushFrame(frame)




tray = I3Tray()

# a random number generator
randomService = phys_services.I3SPRNGRandomService(
    seed = options.SEED,
    nstreams = 10000,
    streamnum = options.RUNNUMBER)

tray.AddModule("I3InfiniteSource","streams",
               Prefix=options.GCDFILE,
               Stream=icetray.I3Frame.DAQ)

tray.AddModule("I3MCEventHeaderGenerator","gen_header",
               Year=2009,
               DAQTime=158100000000000000,
               RunNumber=1,
               EventID=1,
               IncrementEventID=True)


tray.AddModule(mySimpleCascade, "injectCascade",
               I3RandomService = randomService,
               Type = dataclasses.I3Particle.ParticleType.EMinus,
               NEvents = options.NUMEVENTS,
               EnergyMin = 1*I3Units.TeV,
               EnergyMax = 1*I3Units.TeV,
               CylinderRadius = 400.*I3Units.m,
               CylinderHeight = 800.*I3Units.m,
               )


if options.APPLYMMC:
    tray.AddSegment(PropagateMuons, "PropagateMuons",
            RandomService = randomService)
               

tray.AddModule("I3Writer","writer",
    Filename = options.OUTFILE)

tray.AddModule("TrashCan", "the can")

tray.Execute(options.NUMEVENTS+3)




