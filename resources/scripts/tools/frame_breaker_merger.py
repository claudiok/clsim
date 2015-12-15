#!/usr/bin/env python
import sys
import os
from icecube import icetray, dataclasses, dataio, phys_services, clsim

from I3Tray import *

from optparse import OptionParser
from os.path import expandvars

from icecube.simprod import segments
from icecube import MuonGun

usage = "usage: %prog [options] "
parser = OptionParser(usage)
parser.add_option("-s", "--seed",type="int",
                  default=12345, dest="SEED", 
                  help="Initial seed for the random number generator")
parser.add_option("-g", "--gcd",
                  default="/home/nwandkowsky/workspace/data/GCD/GeoCalibDetectorStatus_IC86.55697_V2.i3",
                  dest="GCDFILE", help="Read geometry from GCDFILE (.i3{.gz} format)")
parser.add_option("-r", "--runnumber", type="int", 
                  default=1, dest="RUNNUMBER", help="The run number for this simulation")
parser.add_option("-n", "--numevents", type="int", 
                  default=2, dest="NUMEVENTS", help="The number of events per run")
parser.add_option("--gamma", type="float", 
                  default=2., dest="GAMMA", help="Power law index to use for generation")
parser.add_option("-o","--output", 
                  default="generator_output.i3", dest="OUTPUT", help="output file name")
parser.add_option("--energy_min", type="float", default=1., 
                  dest="EMIN", help="minimal energy in units of TeV")
parser.add_option("--energy_max", type="float", default=10., 
                  dest="EMAX", help="maximal energy in units of PeV")
parser.add_option("--detector", default="IC86", 
                  dest="DETECTOR",help="detector configuration")
parser.add_option("--low-energy-mode",  action="store_true", 
                  default=False, dest="LOWENERGYMODE", help="simulate at lower energies (applies only to MuonGun)")
parser.add_option("--icemodel",default="SpiceLea", 
                  dest="ICEMODEL", 
                  help="Either Spice1, SpiceMie or SpiceLea")

(options,args) = parser.parse_args()

outfile = options.OUTPUT
# input_files = sys.argv[1:-1]
# output_file = sys.argv[-1]
icetray.set_log_level(icetray.I3LogLevel.LOG_WARN)
@icetray.traysegment
def CLSim(tray, name, InputMCTree='I3MCTree', DropMCTree=100, OutputPESeriesMapName='I3MCPESeriesMap',
    IceModel='SpiceLea', UseGPUs=False, MaxParallelEvents=100, UnshadowedFractions=(0.9,0.99,1.08),
    HybridMode=False, DOMOversizeFactor = 5, UseGEANT4 = False):
    """
    Propagate photons to DOMs, prescaling to multiple DOM efficiencies along the way.
    
    :param DropMCTree: Delete the propagated I3MCTree, leaving only the original particles
                       without secondary losses. If this parameter is set to an integer N, drop
                       on all but 1 out of every N frames, otherwise drop unconditionally.
    """
    from icecube import icetray, clsim
    from os.path import expandvars
    
    tray.AddModule("I3MuonSlicer", name + "_chopMuons",
                   InputMCTreeName=InputMCTree,
                   MMCTrackListName="MMCTrackList",
                   OutputMCTreeName=InputMCTree + "_sliced")

    tray.AddModule(clsim.I3FrameSplitterMerger.I3FrameEnergySplitter, "Splitter",
                   InputMCTreeName = InputMCTree + "_sliced")
                   
    tray.AddModule("I3Writer", "writer2",
        Filename = outfile[:-3] + "_after_split.i3",
        Streams = [icetray.I3Frame.Physics, icetray.I3Frame.DAQ,
                   icetray.I3Frame.Stream("q"),
                   icetray.I3Frame.Stream("S")])
                   
    def GetEnergyLost(frame, mctreename = "I3MCTree"):
        energy_counter = 0
        mctree = frame[mctreename]
        for p in mctree:
            if p.shape == dataclasses.I3Particle.Dark or \
               p.location_type != dataclasses.I3Particle.InIce:
               continue
            if p.type == dataclasses.I3Particle.MuMinus or \
               p.type == dataclasses.I3Particle.MuPlus:
                   energy_counter += 45.*icetray.I3Units.GeV
            else:
                energy_counter += p.energy
        print (energy_counter / (1.*icetray.I3Units.PeV))
            
    tray.AddModule(GetEnergyLost, "check",
                   mctreename = InputMCTree + "_sliced", 
                   Streams = [icetray.I3Frame.Stream("q")])
    
    # tray.Add("Dump")
    
    RandomService = tray.context['I3RandomService']

    table_base = expandvars('$I3_DATA/photon-tables/splines/')
    if IceModel == "Spice1":
        clsimIceModel = expandvars("$I3_SRC/clsim/resources/ice/spice_1")
        table_base += "ems_spice1_z20_a10.%s.fits"
    elif IceModel == "SpiceMie":
        clsimIceModel = expandvars("$I3_SRC/clsim/resources/ice/spice_mie")
        table_base += "ems_mie_z20_a10.%s.fits"
    elif IceModel == "SpiceLea":
        clsimIceModel = expandvars("$I3_SRC/clsim/resources/ice/spice_lea")
        table_base += "ems_lea_z20_a10.%s.fits"
    else:
        raise RuntimeError("Unknown ice model: %s", IceModel)


    # Intermediate objects to be deleted at the end of the segment
    temporaries = []

    if HybridMode:
        tray.AddModule("I3MCTreeHybridSimulationSplitter", name+"_splitMCTree",
            InputMCTreeName=InputMCTree,
            OutputMCTreeNameTracks=InputMCTree+"Tracks",
            OutputMCTreeNameCascades=InputMCTree+"Cascades")
        temporaries += [InputMCTree+"Tracks", InputMCTree+"Cascades"]
        CLSimMCTree = InputMCTree+"Tracks"
    else:
        CLSimMCTree = InputMCTree

    # AMD's OpenCL implemenation starts one thread for each core. If taskset is
    # being used to pin the parent process to a specific CPU, then the Linux
    # scheduler may in some circumstances schedule all threads on the same core,
    # resulting in 100%/N CPU usage rather than 100%. Start a background
    # thread that will reset the CPU affinity once the OpenCL threads are
    # spawned (1 minute should be enough).
    if not UseGPUs:
        from icecube.simprod.segments.HybridPhotonicsCLSim import tasksetInUse, resetTasksetThreads
        from threading import Thread
        if tasksetInUse():
            Thread(target=resetTasksetThreads,args=(os.getpid(),)).start()

    # tray.AddModule("Dump")

    # simulate tracks (with clsim)
    tray.AddSegment(clsim.I3CLSimMakePhotons, name+"_makeCLSimHits",
        PhotonSeriesName=name+"_intermediatePhotons",
        MCTreeName = CLSimMCTree + "_sliced",
        # OutputMCTreeName = CLSimMCTree + "_sliced",
        # MCPESeriesName = OutputPESeriesMapName,
        # MMCTrackListName = "MMCTrackList",
        MMCTrackListName = "",
        ParallelEvents = MaxParallelEvents,
        RandomService = RandomService,
        UnshadowedFraction=max(UnshadowedFractions),
        DoNotParallelize=not UseGPUs, # you may need to turn this on for clusters that assume "1 job == 1 core"
        UseGeant4=UseGEANT4, # never use this with Geant4!
        UseGPUs=UseGPUs,
        UseCPUs=not UseGPUs,
        IceModelLocation=clsimIceModel,
        DisableTilt=HybridMode,
        DOMOversizeFactor=DOMOversizeFactor
        )
    temporaries.append(name+"_intermediatePhotons")
    # temporaries.append(CLSimMCTree+"_sliced")


    # tray.AddModule("Dump")

    # now, prescale photons to make MCPEs for each DOM efficiency
    outputs = []
    for eff in UnshadowedFractions:
        label = "%s_%.3f" % (OutputPESeriesMapName, eff)
        outputs.append(label)
        tray.AddSegment(clsim.I3CLSimMakeHitsFromPhotons, name+"_makePhotons_%.3f" % (eff),
            MCTreeName=CLSimMCTree+"_sliced", PhotonSeriesName=name+"_intermediatePhotons",
            MCPESeriesName=label, RandomService=RandomService, UnshadowedFraction=eff)

    # draw cascade photons from spline tables
    if HybridMode:
        from icecube import photonics_service
        cascade_service = photonics_service.I3PhotoSplineService(
            table_base % "abs", table_base % "prob", 0.)
        for eff, hitlabel in zip(UnshadowedFractions, outputs):
            tray.AddModule("I3PhotonicsHitMaker", name+"_hitsFromTheTable_%f" % eff,
                CascadeService = cascade_service,
                TrackService = None, # tracks are handled by clsim
                UnshadowedFraction = eff,
                Input = InputMCTree+"Cascades",
                Output = hitlabel+"Cascades",
                RandomService = RandomService
                )
            temporaries.append(hitlabel+"Cascades")
            tray.Add("Rename", keys=[hitlabel, hitlabel+"Tracks"])
            tray.AddModule("I3CombineMCPE", name+"_combine_pes_%f" % eff,
                InputResponses = [hitlabel+"Tracks", hitlabel+"Cascades"],
                OutputResponse = hitlabel)

    tray.AddModule("Delete", name+"_cleanup",
        Keys = temporaries)

    if DropMCTree:
        if isinstance(DropMCTree, int):
            prescale = lambda frame: frame.Stop == frame.DAQ and RandomService.uniform(0, DropMCTree) > 1
        else:
            prescale = None
        tray.AddModule("Delete", name+"_mctree_cleanup", Keys=[InputMCTree, "MMCTrackList"], If=prescale)

    return outputs

tray = I3Tray()

# tray.AddModule("I3Reader","reader",
#     FileNameList = input_file)

randomService = phys_services.I3GSLRandomService(seed = 1)
prandomService = phys_services.I3GSLRandomService(seed = 2)


tray.context['I3RandomService'] = randomService

eventGenerationTime = dataclasses.I3Time(2011, 158100000000000000)

# 1. generate muons
print("Generating muons...")
tray.AddSegment(segments.GenerateEmptyEvents, "GenerateEmptyEvents",
    RandomService = randomService,
    GCDFile = options.GCDFILE,
    RunNumber = options.RUNNUMBER,
    NumEvents = options.NUMEVENTS,
    FromTime = eventGenerationTime,
    ToTime = eventGenerationTime,
    )

tray.AddSegment(segments.GenerateSingleMuons, "GenerateCosmicRayMuons",
    #RandomService = randomService, #tray.context['I3RandomService']
    NumEvents = options.NUMEVENTS,
    FromEnergy     = options.EMIN*I3Units.PeV,
    ToEnergy       = options.EMAX*I3Units.PeV,
    GammaIndex = options.GAMMA,
    )

# import PropagateMuons to allow regeneration
from icecube.simprod.segments.PropagateMuons import PropagateMuons

tray.AddSegment(PropagateMuons, RandomService=prandomService)

# tray.Add("Dump")

# tray.AddModule(clsim.I3FrameSplitterMerger.I3FrameStreamChanger, "ChangeStream",
#                NewStream = icetray.I3Frame.Stream('q'),
#                OldStream = icetray.I3Frame.Stream('Q'))
              
# tray.Add("Dump")

# tray.AddModule(clsim.I3FrameSplitterMerger.I3FrameEnergySplitter, "Splitter")

tray.AddSegment(CLSim, InputMCTree='I3MCTree',
    IceModel ='SpiceLea', UseGPUs=True, MaxParallelEvents=1,
       # IceModel =options.ICEMODEL, UseGPUs=False, MaxParallelEvents=1,
       OutputPESeriesMapName='I3MCPESeriesMap',
       # UnshadowedFractions=(0.9,0.99,1.08),
       UnshadowedFractions=(0.792,0.8415,0.9,0.9405,0.99,1.0395,1.089,1.1385, 1.188,1.2375,1.287),
    HybridMode=False, DropMCTree=False)

tray.Add("Dump")


tray.AddModule("I3Writer", "writer1",
    Filename = outfile[:-3] + "_before_merge.i3",
    Streams = [icetray.I3Frame.Physics, icetray.I3Frame.DAQ,
               icetray.I3Frame.Stream("q"),
               icetray.I3Frame.Stream("S")])

tray.AddModule(clsim.I3FrameSplitterMerger.I3FrameMCPEMerger, "Merger",
               InputMCPESeriesMapName = "I3MCPESeriesMap_0.990")


# tray.AddModule(clsim.I3FrameSplitterMerger.I3FrameStreamChanger, "ChangeStreamend1",
#                NewStream = icetray.I3Frame.Stream('Q'),
#                OldStream = icetray.I3Frame.Stream('q'))

# tray.AddModule(clsim.I3FrameSplitterMerger.I3FrameStreamChanger, "ChangeStreamend2",
#                NewStream = icetray.I3Frame.Stream('q'),
#                OldStream = icetray.I3Frame.Stream('Q'))

# tray.Add("Dump")


tray.AddModule("I3Writer", "writer",
    Filename = outfile,
    Streams = [icetray.I3Frame.Physics, icetray.I3Frame.DAQ,
               icetray.I3Frame.Stream("q"),
               icetray.I3Frame.Stream("S")])
 
tray.AddModule("TrashCan", "the can")

tray.Execute()
tray.Finish()


