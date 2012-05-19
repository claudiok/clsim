#!/usr/bin/env python
"""
This generates photons originating from a low-energy cascade
for photon table-making.
"""

from optparse import OptionParser
import os
import string

usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)
parser.add_option("-o", "--outfile", default="photon_table.fits",
                  dest="OUTFILE", help="Write output table to OUTFILE (.fits)")
parser.add_option("-s", "--seed",type="int",default=12345,
                  dest="SEED", help="Initial seed for the random number generator")
parser.add_option("-r", "--runnumber", type="int", default=1,
                  dest="RUNNUMBER", help="The run number for this simulation (used for the RNG seed)")

# parse cmd line args, bail out if anything is not understood
(options,args) = parser.parse_args()
if len(args) != 0:
        crap = "Got undefined options:"
        for a in args:
                crap += a
                crap += " "
        parser.error(crap)

########################
if not options.OUTFILE:
        print "No output file!"
        parser.print_help()
        exit(-1)


########################


from I3Tray import *
from os.path import expandvars
import os
import sys

from icecube import icetray, dataclasses, dataio, phys_services, sim_services
from icecube import clsim

import tabulator

# a random number generator
randomService = phys_services.I3SPRNGRandomService(
    seed = options.SEED,
    nstreams = 10000,
    streamnum = options.RUNNUMBER)


tray = I3Tray()

# this is how you can dump some of the simulation timings&statistics to an XML file:
tray.AddService("I3XMLSummaryServiceFactory","summary",
    OutputFileName = "makePhotonsForTabulation.xml")


class MakeParticle(icetray.I3Module):
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddParameter("I3RandomService", "the service", None)
        self.AddParameter("Type", "", dataclasses.I3Particle.ParticleType.EMinus)
        self.AddParameter("Energy", "", 1.*I3Units.GeV)
        self.AddParameter("Zenith", "", 90.*I3Units.degree)
        self.AddParameter("Azimuth", "", 0.*I3Units.degree)
        self.AddParameter("ZCoordinate", "", 0. * I3Units.m)
        self.AddParameter("NEvents", "", 100)

        self.AddOutBox("OutBox")        

    def Configure(self):
        self.rs = self.GetParameter("I3RandomService")
        self.particleType = self.GetParameter("Type")
        self.energy = self.GetParameter("Energy")
        self.zenith = self.GetParameter("Zenith")
        self.azimuth = self.GetParameter("Azimuth")
        self.zCoordinate = self.GetParameter("ZCoordinate")
        self.nEvents = self.GetParameter("NEvents")
        
        self.emittedEvents=0

    def DAQ(self, frame):
        daughter = dataclasses.I3Particle()
        daughter.type = self.particleType
        daughter.energy = self.energy
        daughter.pos = dataclasses.I3Position(0., 0., self.zCoordinate)
        daughter.dir = dataclasses.I3Direction(self.zenith,self.azimuth)
        daughter.time = 0.
        daughter.location_type = dataclasses.I3Particle.LocationType.InIce

        primary = dataclasses.I3Particle()
        primary.type = dataclasses.I3Particle.ParticleType.NuE
        primary.energy = self.energy
        primary.pos = dataclasses.I3Position(0., 0., self.zCoordinate)
        primary.dir = dataclasses.I3Direction(self.zenith,self.azimuth)
        primary.time = 0.
        primary.location_type = dataclasses.I3Particle.LocationType.Anywhere

        mctree = dataclasses.I3MCTree()
        mctree.add_primary(primary)
        mctree.append_child(primary,daughter)
    
        # clsim likes I3MCTrees
        frame["I3MCTree"] = mctree

        # the table-making module prefers plain I3Particles
        frame["Source"] = daughter

        self.emittedEvents += 1
        
        self.PushFrame(frame)
        
        print self.emittedEvents
        if self.emittedEvents >= self.nEvents:
            self.RequestSuspension()



tray.AddModule("I3InfiniteSource","streams",
               Stream=icetray.I3Frame.DAQ)

tray.AddModule("I3MCEventHeaderGenerator","gen_header",
               Year=2009,
               DAQTime=158100000000000000,
               RunNumber=1,
               EventID=1,
               IncrementEventID=True)

# empty GCD
tray.AddService("I3EmptyStreamsFactory", "empty_streams",
    InstallCalibration=True,
    InstallGeometry=True,
    InstallStatus=True)
tray.AddModule("I3MetaSynth", "synth")


tray.AddModule(MakeParticle, "MakeParticle",
    Zenith = 90.*I3Units.deg,
    ZCoordinate = 0.*I3Units.m,
    Energy = 0.01*I3Units.GeV, # the longitudinal profile parameterization from PPC breaks down at this energy (the cascade will be point-like). The angular parameteriztion will still work okay. This might be exactly what we want for table-making (it should become an option to the PPC parameterization eventually).
    NEvents = 2)


tray.AddSegment(clsim.I3CLSimMakePhotons, "makeCLSimPhotons",
    PhotonSeriesName = "PropagatedPhotons",
    MCTreeName = "I3MCTree",
    MMCTrackListName = None,    # do NOT use MMC
    ParallelEvents = 1,         # only work at one event at a time (it'll take long enough)
    RandomService = randomService,
    UseGPUs=False,              # table-making is not a workload particularly suited to GPUs
    UseCPUs=True,               # it should work fine on CPUs, though
    UnshadowedFraction=1.0,     # no cable shadow
    DOMOversizeFactor=1.0,      # no oversizing (there are no DOMs, so this is pointless anyway)
    StopDetectedPhotons=False,  # do not stop photons on detection (also somewhat pointless without DOMs)
    PhotonHistoryEntries=10000, # record all photon paths
    DoNotParallelize=True,      # no multithreading
    OverrideApproximateNumberOfWorkItems=1, # if you *would* use multi-threading, this would be the maximum number of jobs to run in parallel (OpenCL is free to split them)
    ExtraArgumentsToI3CLSimModule=dict(SaveAllPhotons=True,                 # save all photons, regardless of them hitting anything
                                       SaveAllPhotonsPrescale=1.,           # do not prescale the generated photons
                                       StatisticsName="I3CLSimStatistics",  # save a statistics object (contains the initial number of photons)
                                       FixedNumberOfAbsorptionLengths=46.,  # this is approx. the number used by photonics (it uses -ln(1e-20))
                                      ),
    IceModelLocation=expandvars("$I3_SRC/clsim/resources/ice/photonics_wham/Ice_table.wham.i3coords.cos094.11jul2011.txt"),
    #IceModelLocation=expandvars("$I3_SRC/clsim/resources/ice/spice_mie"),
    )


if os.path.exists(options.OUTFILE):
    os.unlink(options.OUTFILE)
tray.AddModule(tabulator.I3TabulatorModule, 'tabulator',
    Photons='PropagatedPhotons', Source='Source', Statistics='I3CLSimStatistics',
    Filename=options.OUTFILE)


#tray.AddModule("I3Writer","writer",
#    Filename = options.OUTFILE)

tray.AddModule("TrashCan", "the can")

tray.Execute()
tray.Finish()
