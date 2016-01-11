# Copyright (c) 2012
# Jakob van Santen <jvansanten@icecube.wisc.edu>
# Claudio Kopper <claudio.kopper@icecube.wisc.edu>
# and the IceCube Collaboration <http://www.icecube.wisc.edu>
#
# Permission to use, copy, modify, and/or distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
# SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
# OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
# CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#
#
# $Id: tabulator.py 140305 2015-12-10 10:54:09Z jvansanten $
#
# @file tabulator.py
# @version $Revision: 140305 $
# @date $Date: 2015-12-10 05:54:09 -0500 (Thu, 10 Dec 2015) $
# @author Jakob van Santen

from __future__ import print_function

from icecube.icetray import I3Units, I3Module, traysegment
from icecube.dataclasses import I3Position, I3Particle, I3MCTree, I3Direction, I3Constants
from icecube.phys_services import I3Calculator, I3GSLRandomService
from icecube.clsim import I3CLSimFunctionConstant
from icecube.clsim import GetIceCubeDOMAcceptance, GetIceCubeDOMAngularSensitivity
from icecube.clsim import FlasherInfoVectToFlasherPulseSeriesConverter, I3CLSimFlasherPulse, I3CLSimFlasherPulseSeries
import numpy, math
from icecube.photospline import numpy_extensions # meshgrid_nd
from icecube.photospline.photonics import FITSTable, Efficiency, Geometry, Parity

def generate_seed():
    import struct
    with open('/dev/random') as rand:
        return struct.unpack('I', rand.read(4))[0]

def makeFlasherPulse(x, y, z, zenith, azimuth, width, brightness, scale):

    pulse = I3CLSimFlasherPulse()
    pulse.type = I3CLSimFlasherPulse.FlasherPulseType.LED405nm
    pulse.pos = I3Position(x, y, z)
    pulse.dir = I3Direction(zenith, azimuth)
    pulse.time = 0.
    pulse.pulseWidth = (float(width)/2.)*I3Units.ns
    lightscale = 1./(32582*5.21) # scale down to match 1 GeV equivalent electromagnetic cascade energy
    pulse.numberOfPhotonsNoBias = 1.17e10*lightscale*scale*(0.0006753+0.00005593*float(brightness))*(float(width)+13.9-(57.5/(1.+float(brightness)/34.4)))

    if numpy.abs(zenith - 90.*I3Units.degree) > 22.5*I3Units.degree:
        tiltedFlasher = True # this is only a rough approximation to describe a tilted flasher
    else:
        tiltedFlasher = False

    pulse.angularEmissionSigmaPolar = FlasherInfoVectToFlasherPulseSeriesConverter.LEDangularEmissionProfile[(pulse.type, tiltedFlasher)][0]
    pulse.angularEmissionSigmaAzimuthal = FlasherInfoVectToFlasherPulseSeriesConverter.LEDangularEmissionProfile[(pulse.type, tiltedFlasher)][1]

    return pulse

def unpin_threads(delay=60):
    """
    When AMD OpenCL fissions the CPU device, it pins each sub-device to a
    a physical core. Since we always use sub-device 0, this means that multiple
    instances of the tabulator on a single machine will compete for core 0.
    Reset thread affinity after *delay* seconds to prevent this from happening.
    """
    import os
    import subprocess
    import threading
    import time
    def taskset(pid,tt=None):
        # get/set the taskset affinity for pid
        # uses a binary number string for the core affinity
        l = ['/bin/taskset','-p']
        if tt:
            l.append(hex(int(tt,2))[2:])
        l.append(str(pid))
        p = subprocess.Popen(l,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        output = p.communicate()[0].split(':')[-1].strip()
        if not tt:
            return bin(int(output,16))[2:]
    
    def resetTasksetThreads(main_pid):
        # reset thread taskset affinity
        time.sleep(delay)
        num_cpus = reduce(lambda b,a: b+int('processor' in a),open('/proc/cpuinfo').readlines(),0)
        tt = '1'*num_cpus
        #tt = taskset(main_pid)
        p = subprocess.Popen(['/bin/ps','-Lo','tid','--no-headers','%d'%main_pid],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        for tid in p.communicate()[0].split():
            tid = tid.strip()
            if tid:
                taskset(tid,tt)
    # only do this on linux
    try:
        open('/proc/cpuinfo')
    except IOError:
        return
    threading.Thread(target=resetTasksetThreads,args=(os.getpid(),)).start()
    

from icecube.clsim import GetDefaultParameterizationList
from icecube.clsim import GetFlasherParameterizationList
from icecube.clsim import GetHybridParameterizationList

from icecube.clsim import AutoSetGeant4Environment

from icecube.clsim.traysegments.common import configureOpenCLDevices, parseIceModel
from icecube import icetray
from os.path import expandvars

@icetray.traysegment
def I3CLSimTabulatePhotons(tray, name,
                       UseCPUs=True,
                       UseGPUs=False,
                       UseOnlyDeviceNumber=None,
                       MCTreeName="I3MCTree",
                       OutputMCTreeName=None,
                       FlasherInfoVectName=None,
                       FlasherPulseSeriesName=None,
                       MMCTrackListName="MMCTrackList",
                       ParallelEvents=1000,
                       RandomService=None,
                       MediumProperties=expandvars("$I3_SRC/clsim/resources/ice/spice_mie"),
                       UseGeant4=False,
                       CrossoverEnergyEM=None,
                       CrossoverEnergyHadron=None,
                       UseCascadeExtension=False,
                       DoNotParallelize=False,
                       DOMOversizeFactor=1,
                       UseHoleIceParameterization=True,
                       OverrideApproximateNumberOfWorkItems=None,
                       ExtraArgumentsToI3CLSimModule=dict(),
                       If=lambda f: True
                       ):
    """Do standard clsim processing up to the I3Photon level.
    These photons still need to be converted to I3MCPEs to be usable
    for further steps in the standard IceCube MC processing chain.
    Reads its particles from an I3MCTree and writes an I3PhotonSeriesMap.

    All available OpenCL GPUs (and optionally CPUs) will
    be used. This will take over your entire machine,
    so make sure to configure your batch jobs correctly
    when using this on a cluster.
    When using nVidia cards, you can set the
    CUDA_VISIBLE_DEVICES environment variable to
    limit GPU visibility. A setting of
    CUDA_VISIBLE_DEVICES="0,3" would only use cards
    #0 and #3 and ignore cards #1 and #2. In case you are
    using a batch system, chances are this variable is already
    set. Unfortunately, there is no corresponding setting
    for the AMD driver.

    This segment assumes that MMC has been applied to the
    I3MCTree and that MMC was *NOT* run using the "-recc" option.

    :param UseCPUs:
        Turn this on to also use CPU-based devices.
        (This may potentially slow down photon generation, which
        is also done on the CPU in parallel.)
    :param UseGPUs:
        Turn this off to not use GPU-based devices.
        This may be useful if your GPU is used for display
        purposes and you don't want it to slow down.
    :param UseOnlyDeviceNumber:
        Use only a single device number, even if there is more than
        one device found matching the required description. The numbering
        starts at 0.
    :param MCTreeName:
        The name of the I3MCTree containing the particles to propagate.
    :param OutputMCTreeName:
        A copy of the (possibly sliced) MCTree will be stored as this name.
    :param FlasherInfoVectName:
        Set this to the name of I3FlasherInfoVect objects in the frame to
        enable flasher simulation. The module will read I3FlasherInfoVect objects
        and generate photons according to assumed parameterizations.
    :param FlasherPulseSeriesName:
        Set this to the name of an I3CLSimFlasherPulseSeries object in the frame to
        enable flasher/Standard Candle simulation.
        This cannot be used at the same time as FlasherInfoVectName.
        (I3CLSimFlasherPulseSeries objects are clsim's internal flasher
        representation, if "FlasherInfoVectName" is used, the I3FlasherInfoVect
        objects are converted to I3CLSimFlasherPulseSeries objects.)
    :param MMCTrackListName:
        Only used if *ChopMuons* is active. Set it to the name
        of the I3MMCTrackList object that contains additional
        muon energy loss information.
    :param PhotonSeriesName:
        Configure this to enable writing an I3PhotonSeriesMap containing
        all photons that reached the DOM surface.
    :param ParallelEvents:
        clsim will work on a couple of events in parallel in order
        not to starve the GPU. Setting this too high will result
        in excessive memory usage (all your frames have to be cached
        in RAM). Setting it too low may impact simulation performance.
        The optimal value depends on your energy distribution/particle type.
    :param RandomService:
        Set this to an instance of a I3RandomService. Alternatively,
        you can specify the name of a configured I3RandomServiceFactory
        added to I3Tray using tray.AddService(). If you don't configure
        this, the default I3RandomServiceFactory will be used.
    :param MediumProperties:
        Set this either to a directory containing a PPC-compatible
        ice description (icemodel.dat, icemodel.par and cfg.txt) or
        to a photonics ice table file. PPC-compatible ice files should
        generally lead to faster execution times on GPUs since it involves
        less interpolation between table entries (the PPC ice-specification
        is parametric w.r.t. wavelength, whereas the photonics specification
        is not).
    :param DisableTilt:
        Do not simulate ice tilt, even if the ice model directory
        provides tilt information. (Photonics-based models will never
        have tilt.)
    :param UnWeightedPhotons:
        Enabling this setting will disable all optimizations. These
        are currently a DOM oversize factor of 5 (with the appropriate
        timing correction) and a biased initial photon spectrum that
        includes the DOM spectral acceptance. Enabling this setting
        essentially means that all photons that would be generated
        in the real detector *will* actually be generated. This will siginificatly
        slow down the simulation, but the optional ``PhotonSeries``
        will contain an unweighted sample of photons that arrive
        at your DOMs. This can be useful for DOM acceptance studies.
    :param StopDetectedPhotons:
        Configures behaviour for photons that hit a DOM. If this is true (the default)
        photons will be stopped once they hit a DOM. If this is false, they continue to
        propagate.
    :param PhotonHistoryEntries:
        The maximum number of scatterings points to be saved for every photon hitting a DOM.
        Only the most recent positions are saved, older positions are overwritten if
        the maximum size is reached.
    :param UseGeant4:
        Enabling this setting will disable all cascade and muon light yield
        parameterizations. All particles will sent to Geant4 for a full
        simulation. This does **not** apply to muons that do have a length
        assigned. These are assumed to have been generated by MMC and
        their light is generated according to the usual parameterization.
    :param CrossoverEnergyEM:
        If set it defines the crossover energy between full Geant4 simulations and 
        light yield parameterizations for electro magnetic cascades. This only works
        when UseGeant4 is set to true. It works in conjunction with CrossoverEnergyHadron.
        If one of both is set to a positiv value greater 0 (GeV), the hybrid simulation
        is used.
        If CrossoverEnergyEM is set to None while CrossoverEnergyHadron is set so
        hybrid mode is working, GEANT4 is used for EM cascades.
        If CrossoverEnergyEM is set to 0 (GeV) while CrossoverEnergyHadron is set
        so hybrid mode is working, leptons and EM cascades will use parameterizations
        for the whole energy range.
    :param CrossoverEnergyHadron:
        If set it defines the crossover energy between full Geant4 simulations and
        light yield parameterizations for hadronic cascades. This only works when
        UseGeant4 is set to true. It works in conjunction with CrossoverEnergyEM.
        If one of both is set to a positiv value greater 0 (GeV), the hybrid simulation
        is used.
        If CrossoverEnergyHadron is set to None while CrossoverEnergyEM is set so
        hybrid mode is working, GEANT4 is used for hadronic cascades.
        If CrossoverEnergyHadron is set to 0 (GeV) while CrossoverEnergyHadron is
        set so hybrid mode is working, hadronic cascades will use parameterizations
        for the whole energy range.
    :param UseCascadeExtension:
        If True, simulate the longitudinal development of cascades. Otherwise,
        simulate cascades as pointlike objects.
    :param DoNotParallelize:
        Try only using a single work item in parallel when running the
        OpenCL simulation. This might be useful if you want to run jobs
        in parallel on a batch system. This will only affect CPUs and
        will be a no-op for GPUs.
    :param DOMOversizeFactor:
        Set the DOM oversize factor. To disable oversizing, set this to 1.
    :param UnshadowedFraction:
        Fraction of photocathode available to receive light (e.g. unshadowed by the cable)
    :param UseHoleIceParameterization:
        Use an angular acceptance correction for hole ice scattering.
    :param OverrideApproximateNumberOfWorkItems:
        Allows to override the auto-detection for the maximum number of parallel work items.
        You should only change this if you know what you are doing.
    :param If:
        Python function to use as conditional execution test for segment modules.        
    """

    from icecube import icetray, dataclasses, phys_services, clsim

    # make sure the geometry is updated to the new granular format (in case it is supported)
    if hasattr(dataclasses, "I3ModuleGeo"):
        tray.AddModule("I3GeometryDecomposer", name + "_decomposeGeometry",
                       If=lambda frame: If(frame) and ("I3OMGeoMap" not in frame))

    if UseGeant4:
        if not clsim.I3CLSimLightSourceToStepConverterGeant4.can_use_geant4:
            raise RuntimeError("You have requested to use Geant 4, but clsim was compiled without Geant 4 support")
    
    # at the moment the Geant4 paths need to be set, even if it isn't used
    # TODO: fix this
    if clsim.I3CLSimLightSourceToStepConverterGeant4.can_use_geant4:
        AutoSetGeant4Environment()

    # some constants
    DOMRadius = 0.16510*icetray.I3Units.m # 13" diameter
    Jitter = 2.*icetray.I3Units.ns

    if MMCTrackListName is None or MMCTrackListName=="":
        # the input does not seem to have been processed by MMC
        ChopMuons = False
    else:
        ChopMuons = True

    if MCTreeName is None or MCTreeName=="":
        clSimMCTreeName=""
        if ChopMuons:
            raise RuntimeError("You cannot have \"MMCTrackListName\" enabled with no MCTree!")
    else:
        clSimMCTreeName=MCTreeName

    if FlasherInfoVectName is None or FlasherInfoVectName=="":
        if (FlasherPulseSeriesName is not None) and (FlasherPulseSeriesName!=""):
            SimulateFlashers=True
            clSimFlasherPulseSeriesName = FlasherPulseSeriesName
            clSimOMKeyMaskName = ""
        else:
            SimulateFlashers=False
            clSimFlasherPulseSeriesName = ""
            clSimOMKeyMaskName = ""
    else:
        if (FlasherPulseSeriesName is not None) and (FlasherPulseSeriesName!=""):
            raise RuntimeError("You cannot use the FlasherPulseSeriesName and FlasherInfoVectName parameters at the same time!")
        
        SimulateFlashers=True
        clSimFlasherPulseSeriesName = FlasherInfoVectName + "_pulses"
        clSimOMKeyMaskName = FlasherInfoVectName + "_OMKeys"
        
        tray.AddModule(clsim.FlasherInfoVectToFlasherPulseSeriesConverter,
                       name + "_FlasherInfoVectToFlasherPulseSeriesConverter",
                       FlasherInfoVectName = FlasherInfoVectName,
                       FlasherPulseSeriesName = clSimFlasherPulseSeriesName,
                       FlasherOMKeyVectName = clSimOMKeyMaskName,
                       If=If)

    # (optional) pre-processing
    if ChopMuons:
        if OutputMCTreeName is not None:
            clSimMCTreeName_new = OutputMCTreeName
        else:
            clSimMCTreeName_new = clSimMCTreeName + "_sliced"
        
        tray.AddModule("I3MuonSlicer", name + "_chopMuons",
                       InputMCTreeName=clSimMCTreeName,
                       MMCTrackListName=MMCTrackListName,
                       OutputMCTreeName=clSimMCTreeName_new,
                       If=If)
        clSimMCTreeName = clSimMCTreeName_new
    else:
        if (OutputMCTreeName is not None) and (OutputMCTreeName != ""):
            # copy the MCTree to the requested output name
            def copyMCTree(frame, inputName, outputName, If_=None):
                if If_ is not None:
                    if not If_(frame): return
                frame[outputName] = frame[inputName]
            tray.AddModule(copyMCTree, name + "_copyMCTree",
                           inputName=clSimMCTreeName,
                           outputName=OutputMCTreeName,
                           Streams=[icetray.I3Frame.DAQ],
                           If_=If)
            clSimMCTreeName = OutputMCTreeName
        else:
            clSimMCTreeName = clSimMCTreeName

    domAcceptance = clsim.GetIceCubeDOMAcceptance(domRadius=DOMOversizeFactor*DOMRadius)
    if UseHoleIceParameterization:
        angularAcceptance = clsim.GetIceCubeDOMAngularSensitivity(UseHoleIceParameterization)
    else:
        icetray.logging.log_warn("Applying *no* angular sensitivity at all. None.")
        angularAcceptance = clsim.I3CLSimFunctionConstant(1.)

    # muon&cascade parameterizations
    ppcConverter = clsim.I3CLSimLightSourceToStepConverterPPC(photonsPerStep=200)
    ppcConverter.SetUseCascadeExtension(UseCascadeExtension)
    if not UseGeant4:
        particleParameterizations = GetDefaultParameterizationList(ppcConverter, muonOnly=False)
    else:
        if CrossoverEnergyEM>0 or CrossoverEnergyHadron>0:
            particleParameterizations = GetHybridParameterizationList(ppcConverter, CrossoverEnergyEM=CrossoverEnergyEM, CrossoverEnergyHadron=CrossoverEnergyHadron)
        elif MMCTrackListName is None or MMCTrackListName=="":
            particleParameterizations = [] # make sure absolutely **no** parameterizations are used
        else:
            # use no parameterizations except for muons with lengths assigned to them
            # (those are assumed to have been generated by MMC)
            particleParameterizations = GetDefaultParameterizationList(ppcConverter, muonOnly=True)
    
    # flasher parameterizations
    if SimulateFlashers:
        # this needs a spectrum table in order to pass spectra to OpenCL
        spectrumTable = clsim.I3CLSimSpectrumTable()
        particleParameterizations += GetFlasherParameterizationList(spectrumTable)
        
        icetray.logging.log_debug("number of spectra (1x Cherenkov + Nx flasher): %d" % len(spectrumTable), unit="clsim")
    else:
        # no spectrum table is necessary when only using the Cherenkov spectrum
        spectrumTable = None

    openCLDevices = configureOpenCLDevices(
        UseGPUs=UseGPUs,
        UseCPUs=UseCPUs,
        OverrideApproximateNumberOfWorkItems=OverrideApproximateNumberOfWorkItems,
        DoNotParallelize=DoNotParallelize,
        UseOnlyDeviceNumber=UseOnlyDeviceNumber
	)
    
    tray.AddModule("I3CLSimTabulatorModule", name + "_clsim",
                   MCTreeName=clSimMCTreeName,
                   RandomService=RandomService,
                   MediumProperties=MediumProperties,
                   SpectrumTable=spectrumTable,
                   FlasherPulseSeriesName=clSimFlasherPulseSeriesName,
                   WavelengthAcceptance=domAcceptance,
                   AngularAcceptance=angularAcceptance,
                   ParameterizationList=particleParameterizations,
                   # MaxNumParallelEvents=ParallelEvents,
                   OpenCLDeviceList=openCLDevices,
                   **ExtraArgumentsToI3CLSimModule
                   )
    
    unpin_threads()

@traysegment
def TabulatePhotonsFromSource(tray, name, PhotonSource="cascade", Zenith=0.*I3Units.degree, Azimuth=0.*I3Units.degree, ZCoordinate=0.*I3Units.m,
    Energy=1.*I3Units.GeV, FlasherWidth=127, FlasherBrightness=127, Seed=12345, NEvents=100,
    IceModel='spice_mie', DisableTilt=False, Filename="", TabulateImpactAngle=False,
    PhotonPrescale=1, Axes=None, Directions=None):
    
    """
    Tabulate the distribution of photoelectron yields on IceCube DOMs from various
    light sources. The light profiles of the sources are computed from the same
    parameterizations used in PPC, but like in the direct propagation mode can
    be computed using GEANT4 instead.

    The mode of tabulation is controlled primarily by the **PhotonSource** parameter.
    
    - *'cascade'* will simulate an electromagnetic cascade of **Energy** GeV at
      (0, 0, **ZCoordinate**), oriented according to **Zenith** and **Azimuth**.
      The default coordinate system is spherical and centered the given vertex,
      with 200 quadratically spaced bins in radius, 36 linear bins in azimuthal
      angle (only from 0 to 180 degrees by default), 100 linear bins in the
      cosine of the polar angle, and 105 quadratic bins in time residual w.r.t
      the direct path from (0, 0, **ZCoordinate**).
    - *'flasher'* will simulate a 405 nm LED flasher pulse with the given
      **FlasherWidth** and **FlasherBrightness** settings. The source position
      and coordinate system are the same as for the 'cascade' case.
    - *'infinite-muon'* will simulate a "bare" muon of infinite length. The
      coordinate system is cylindrical and centered on the axis of the muon.
      Since the muon's position is degenerate with time, the usual parallel
      distance is replaced by the z coordinate of the closest approach to the
      detection position, and the starting positions of the simulated muons are
      sampled randomly (**ZCoordinate** is ignored). There are 100 quadratic
      bins in perpendicular distance to the source axis, 36 linear bins in
      azimuthal angle (0 to :math:`\pi` radians), 100 linear bins in z
      coordinate of closest approach, and 105 quadratic bins in time residual
      w.r.t. the earliest possible Cherenkov photon.

    :param PhotonSource: the type of photon source ('cascade', 'flasher', or 'infinite-muon').
    :param Zenith: the orientation of the source
    :param ZCoordinate: the depth of the source
    :param Energy: the energy of the source (only for cascade tables)
    :param FlasherWidth: the width of the flasher pulse (only for flasher tables)
    :param FlasherBrightness: the brightness of the flasher pulse (only for flasher tables)
    :param Seed: the seed for the random number service
    :param NEvents: the number of events to simulate
    :param IceModel: the path to an ice model in $I3_SRC/clsim/resources/ice. Likely values include:
        'spice_mie' ppc-style SPICE-Mie parametrization
        'photonics_spice_1/Ice_table.spice.i3coords.cos080.10feb2010.txt' Photonics-style SPICE1 table
        'photonics_wham/Ice_table.wham.i3coords.cos094.11jul2011.txt' Photonics-style WHAM! table
    :param DisableTilt: if true, disable tilt in ice model
    :param Filename: the name of the FITS file to write
    :param TabulateImpactAngle: if True, tabulate the impact position of the
           photon on the DOM instead of weighting by the DOM's angular acceptance
    :param Axes: a subclass of :cpp:class:`clsim::tabulator::Axes` that defines the coordinate system.
                 If None, an appropriate default will be chosen based on **PhotonSource**.
    :param Directions: a set of directions to allow table generation for multiple sources.
                 If None, only one direction given by **Zenith** and **Azimuth** is used.
       """

    # check sanity of args
    PhotonSource = PhotonSource.lower()
    if PhotonSource not in ['cascade', 'flasher', 'infinite-muon']:
        raise ValueError("photon source %s is unknown. Please specify either 'cascade', 'flasher', or 'infinite-muon'" % PhotonSource)
    
    from icecube import icetray, dataclasses, dataio, phys_services, sim_services, clsim
    from os.path import expandvars
    
    # a random number generator
    randomService = phys_services.I3GSLRandomService(Seed)
        
    tray.AddModule("I3InfiniteSource",name+"streams",
                   Stream=icetray.I3Frame.DAQ)

    tray.AddModule("I3MCEventHeaderGenerator",name+"gen_header",
                   Year=2009,
                   DAQTime=158100000000000000,
                   RunNumber=1,
                   EventID=1,
                   IncrementEventID=True)

    if Directions is None:
        Directions = numpy.asarray([(Zenith, Azimuth)])

    if PhotonSource == 'cascade' or PhotonSource == 'flasher':

        ptype = I3Particle.ParticleType.EMinus

        def reference_source(zenith, azimuth, scale):
            source = I3Particle()
            source.type = ptype
            source.energy = Energy*scale
            source.pos = I3Position(0., 0., ZCoordinate)
            source.dir = I3Direction(zenith, azimuth)
            source.time = 0.
            source.length = 0.
            source.location_type = I3Particle.LocationType.InIce
        
            return source
    
    elif PhotonSource == 'infinite-muon':
        
        from icecube import MuonGun
        surface = MuonGun.Cylinder(1600, 800)
        
        ptype = I3Particle.ParticleType.MuMinus
        
        def reference_source(zenith, azimuth, scale):
            source = I3Particle()
            source.type = ptype
            source.energy = Energy*scale
            source.dir = I3Direction(zenith, azimuth)
            source.pos = surface.sample_impact_position(source.dir, randomService)
            crossings = surface.intersection(source.pos, source.dir)
            source.length = crossings.second-crossings.first
            source.time = 0.
            source.location_type = I3Particle.LocationType.InIce
            
            return source
    
    import copy
    
    class MakeParticle(icetray.I3Module):
        def __init__(self, ctx):
            super(MakeParticle,self).__init__(ctx)
            self.AddOutBox("OutBox")
            self.AddParameter("SourceFunction", "", lambda : None)
            self.AddParameter("NEvents", "", 100)
        def Configure(self):
            self.reference_source = self.GetParameter("SourceFunction")
            self.nevents = self.GetParameter("NEvents")
            self.emittedEvents = 0
        def DAQ(self, frame):
            if PhotonSource != "flasher":
                primary = I3Particle()
                mctree = I3MCTree()
                mctree.add_primary(primary)
                for zenith, azimuth in Directions:
                    source = self.reference_source(zenith, azimuth, 1./len(Directions))
                    mctree.append_child(primary, source)
                frame["I3MCTree"] = mctree
            else:
                pulseseries = I3CLSimFlasherPulseSeries()
                for zenith, azimuth in Directions:
                    pulse = makeFlasherPulse(0, 0, ZCoordinate, zenith, azimuth, FlasherWidth, FlasherBrightness, 1./len(Directions))
                    pulseseries.append(pulse)
                frame["I3FlasherPulseSeriesMap"] = pulseseries

            # use the primary particle as a geometrical reference
            frame["ReferenceParticle"] = self.reference_source(Zenith, Azimuth, 1.)
            
            self.PushFrame(frame)
            
            self.emittedEvents += 1
            if self.emittedEvents >= self.nevents:
                self.RequestSuspension()

    tray.AddModule(MakeParticle, SourceFunction=reference_source, NEvents=NEvents)

    if PhotonSource == "flasher":
        flasherpulse = "I3FlasherPulseSeriesMap"
        mctree = None
    else:
        flasherpulse = None
        mctree = "I3MCTree"
    
    header = dict(FITSTable.empty_header)
    header['zenith'] = Zenith/I3Units.degree
    header['azimuth'] = Azimuth/I3Units.degree
    header['z'] = ZCoordinate
    header['energy'] = Energy
    header['type'] = int(ptype)
    header['efficiency'] = Efficiency.RECEIVER | Efficiency.WAVELENGTH
    
    if Axes is None:
        if PhotonSource != "infinite-muon":
            dims = [
                clsim.tabulator.PowerAxis(0, 580, 200, 2),
                clsim.tabulator.LinearAxis(0, 180, 36),
                clsim.tabulator.LinearAxis(-1, 1, 100),
                clsim.tabulator.PowerAxis(0, 7e3, 105, 2),
            ]
            geo = clsim.tabulator.SphericalAxes
        else:
            dims = [
                clsim.tabulator.PowerAxis(0, 580, 100, 2),
                clsim.tabulator.LinearAxis(0, numpy.pi, 36),
                clsim.tabulator.LinearAxis(-8e2, 8e2, 80),
                clsim.tabulator.PowerAxis(0, 7e3, 105, 2),
            ]
            geo = clsim.tabulator.CylindricalAxes
        # Add a dimension for the impact angle
        if TabulateImpactAngle:
            dims.append(clsim.tabulator.LinearAxis(-1, 1, 20))
        Axes = geo(dims)

    if PhotonSource == "flasher":
        header['flasherwidth'] = FlasherWidth
        header['flasherbrightness'] = FlasherBrightness

    tray.AddSegment(I3CLSimTabulatePhotons, name+"makeCLSimPhotons",
        MCTreeName = mctree,                        # if source is a cascade this will point to the I3MCTree
        FlasherPulseSeriesName = flasherpulse,      # if source is a flasher this will point to the I3CLSimFlasherPulseSeries
        MMCTrackListName = None,                    # do NOT use MMC
        ParallelEvents = 1,                         # only work at one event at a time (it'll take long enough)
        RandomService = randomService,
        # UnWeightedPhotons=True,
        UseGPUs=False,                              # table-making is not a workload particularly suited to GPUs
        UseCPUs=True,                               # it should work fine on CPUs, though
        DOMOversizeFactor=math.sqrt(PhotonPrescale),
        UseHoleIceParameterization=True,
        DoNotParallelize=True,                      # no multithreading
        UseGeant4=False,
        OverrideApproximateNumberOfWorkItems=1,     # if you *would* use multi-threading, this would be the maximum number of jobs to run in parallel (OpenCL is free to split them)
        ExtraArgumentsToI3CLSimModule=dict(Filename=Filename, TableHeader=header,
            Axes=Axes, PhotonsPerBunch=200, EntriesPerPhoton=5000),
        MediumProperties=parseIceModel(expandvars("$I3_SRC/clsim/resources/ice/" + IceModel), disableTilt=DisableTilt),
    )
