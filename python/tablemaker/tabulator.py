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
# $Id$
#
# @file tabulator.py
# @version $Revision$
# @date $Date$
# @author Jakob van Santen

from __future__ import print_function

from icecube.icetray import I3Units, I3Module, traysegment
from icecube.dataclasses import I3Position, I3Particle, I3MCTree, I3Direction, I3Constants
from icecube.phys_services import I3Calculator, I3GSLRandomService
from icecube.clsim import I3Photon, I3CLSimTabulator, GetIceCubeDOMAcceptance, GetIceCubeDOMAngularSensitivity
from icecube.clsim import FlasherInfoVectToFlasherPulseSeriesConverter, I3CLSimFlasherPulse, I3CLSimFlasherPulseSeries
from icecube.clsim.traysegments.common import parseIceModel
from icecube.clsim.util import GetRefractiveIndexRange
import numpy, math
from icecube.photospline import numpy_extensions # meshgrid_nd
from icecube.photospline.photonics import FITSTable, Efficiency, Geometry, Parity

def generate_seed():
    import struct
    with open('/dev/random') as rand:
        return struct.unpack('I', rand.read(4))[0]

class MakeParticle(I3Module):
    def __init__(self, context):
        I3Module.__init__(self, context)
        self.AddParameter("I3RandomService", "the service", None)
        self.AddParameter("PhotonSource", "", "CASCADE")
        self.AddParameter("ParticleType", "", I3Particle.ParticleType.EMinus)
        self.AddParameter("Energy", "", 1.*I3Units.GeV)
        self.AddParameter("Zenith", "", 90.*I3Units.degree)
        self.AddParameter("Azimuth", "", 0.*I3Units.degree)
        self.AddParameter("ZCoordinate", "", 0.*I3Units.m)
        self.AddParameter("FlasherWidth", "", 127)
        self.AddParameter("FlasherBrightness", "", 127)
        self.AddParameter("NEvents", "", 100)

        self.AddOutBox("OutBox")        

    def Configure(self):
        self.rs = self.GetParameter("I3RandomService")
        self.photonSource = self.GetParameter("PhotonSource")
        self.particleType = self.GetParameter("ParticleType")
        self.energy = self.GetParameter("Energy")
        self.zenith = self.GetParameter("Zenith")
        self.azimuth = self.GetParameter("Azimuth")
        self.zCoordinate = self.GetParameter("ZCoordinate")
        self.flasherWidth = self.GetParameter("FlasherWidth")
        self.flasherBrightness = self.GetParameter("FlasherBrightness")
        self.nEvents = self.GetParameter("NEvents")
        
        self.emittedEvents=0

    def DAQ(self, frame):

        # create source particle
        source = I3Particle()
        source.type = self.particleType
        source.energy = self.energy
        source.pos = I3Position(0., 0., self.zCoordinate)
        source.dir = I3Direction(self.zenith, self.azimuth)
        source.time = 0.
        source.length = 0.
        source.location_type = I3Particle.LocationType.InIce

        # the table-making module prefers plain I3Particles
        frame["Source"] = source

        if self.photonSource.upper() == "CASCADE":

            primary = I3Particle()
            primary.type = I3Particle.ParticleType.NuE
            primary.energy = self.energy
            primary.pos = I3Position(0., 0., self.zCoordinate)
            primary.dir = I3Direction(self.zenith, self.azimuth)
            primary.time = 0.
            primary.location_type = I3Particle.LocationType.Anywhere

            mctree = I3MCTree()
            mctree.add_primary(primary)
            mctree.append_child(primary, source)
    
            # clsim likes I3MCTrees
            frame["I3MCTree"] = mctree
            
        elif self.photonSource.upper() == "FLASHER":

            outputSeries = I3CLSimFlasherPulseSeries()
            flasherPulse = self.getFlasherPulse(0., 0., self.zCoordinate, self.zenith, self.azimuth)
            outputSeries.append(flasherPulse)

            # clsim likes I3CLSimFlasherPulseSeries
            frame["I3FlasherPulseSeriesMap"] =  outputSeries

        self.emittedEvents += 1
        
        self.PushFrame(frame)
        
        print(self.emittedEvents)
        if self.emittedEvents >= self.nEvents:
            self.RequestSuspension()
    
def makeFlasherPulse(x, y, z, zenith, azimuth, width, brightness):

    pulse = I3CLSimFlasherPulse()
    pulse.type = I3CLSimFlasherPulse.FlasherPulseType.LED405nm
    pulse.pos = I3Position(x, y, z)
    pulse.dir = I3Direction(zenith, azimuth)
    pulse.time = 0.
    pulse.pulseWidth = (float(width)/2.)*I3Units.ns
    scale = 32582*5.21 # scale down to match 1 GeV equivalent electromagnetic cascade energy
    pulse.numberOfPhotonsNoBias = 1.17e10/scale*(0.0006753+0.00005593*float(brightness))*(float(width)+13.9-(57.5/(1.+float(brightness)/34.4)))

    if numpy.abs(zenith - 90.*I3Units.degree) > 22.5*I3Units.degree:
        tiltedFlasher = True # this is only a rough approximation to describe a tilted flasher
    else:
        tiltedFlasher = False

    pulse.angularEmissionSigmaPolar = FlasherInfoVectToFlasherPulseSeriesConverter.LEDangularEmissionProfile[(pulse.type, tiltedFlasher)][0]
    pulse.angularEmissionSigmaAzimuthal = FlasherInfoVectToFlasherPulseSeriesConverter.LEDangularEmissionProfile[(pulse.type, tiltedFlasher)][1]

    return pulse

from icecube.clsim import GetDefaultParameterizationList
from icecube.clsim import GetFlasherParameterizationList
from icecube.clsim import GetHybridParameterizationList

from icecube.clsim import AutoSetGeant4Environment

from icecube.clsim.traysegments.common import configureOpenCLDevices, parseIceModel
from icecube import icetray
from os.path import expandvars

@icetray.traysegment
def I3CLSimTabulatePhotons(tray, name,
                       UseCPUs=False,
                       UseGPUs=True,
                       UseOnlyDeviceNumber=None,
                       MCTreeName="I3MCTree",
                       OutputMCTreeName=None,
                       FlasherInfoVectName=None,
                       FlasherPulseSeriesName=None,
                       MMCTrackListName="MMCTrackList",
                       PhotonSeriesName="PhotonSeriesMap",
                       ParallelEvents=1000,
                       RandomService=None,
                       IceModelLocation=expandvars("$I3_SRC/clsim/resources/ice/spice_mie"),
                       DisableTilt=False,
                       UseGeant4=False,
                       CrossoverEnergyEM=None,
                       CrossoverEnergyHadron=None,
                       StopDetectedPhotons=True,
                       PhotonHistoryEntries=0,
                       DoNotParallelize=False,
                       UnshadowedFraction=0.9,
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
    :param IceModelLocation:
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

    # ice properties
    if isinstance(IceModelLocation, str):
        mediumProperties = parseIceModel(IceModelLocation, disableTilt=DisableTilt)
    else:
        # get ice model directly if not a string
        mediumProperties = IceModelLocation
    
    # The photonics convention for photon table normalization is "per Cherenkov
    # photon between 400 and 600 nm," of which there are 32582.0 per meter of
    # track. Since we generate over a wider range, re-normalize to compensate.
    lightScale = 32582.0/clsim.NumberOfPhotonsPerMeter(mediumProperties.GetPhaseRefractiveIndex(0), clsim.I3CLSimFunctionConstant(1.),
        mediumProperties.MinWavelength, mediumProperties.MaxWavelength)

    domAcceptance = clsim.GetIceCubeDOMAcceptance(efficiency=lightScale)
    angularAcceptance = clsim.GetIceCubeDOMAngularSensitivity(UseHoleIceParameterization)

    # muon&cascade parameterizations
    ppcConverter = clsim.I3CLSimLightSourceToStepConverterPPC(photonsPerStep=200)
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
        
        print("number of spectra (1x Cherenkov + Nx flasher):", len(spectrumTable))
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
                   # MCTreeName=clSimMCTreeName,
                   # PhotonSeriesMapName=PhotonSeriesName,
                   # DOMRadius = DOMRadius,
                   # DOMOversizeFactor = DOMOversizeFactor,
                   # DOMPancakeFactor = DOMOversizeFactor, # you will probably want this to be the same as DOMOversizeFactor
                   RandomService=RandomService,
                   MediumProperties=mediumProperties,
                   # SpectrumTable=spectrumTable,
                   # FlasherPulseSeriesName=clSimFlasherPulseSeriesName,
                   # OMKeyMaskName=clSimOMKeyMaskName,
                   # ignore IceTop
                   # IgnoreSubdetectors = ["IceTop"],
                   #IgnoreNonIceCubeOMNumbers=False,
                   # GenerateCherenkovPhotonsWithoutDispersion=False,
                   WavelengthAcceptance=domAcceptance,
                   AngularAcceptance=angularAcceptance,
                   ParameterizationList=particleParameterizations,
                   # MaxNumParallelEvents=ParallelEvents,
                   OpenCLDeviceList=openCLDevices,
                   #UseHardcodedDeepCoreSubdetector=False, # setting this to true saves GPU constant memory but will reduce performance
                   # StopDetectedPhotons=StopDetectedPhotons,
                   # PhotonHistoryEntries=PhotonHistoryEntries,
                   # If=If,
                   **ExtraArgumentsToI3CLSimModule
                   )

@traysegment
def CombinedPhotonGenerator(tray, name, PhotonSource="CASCADE", Zenith=90.*I3Units.degree, Azimuth=0*I3Units.degree, ZCoordinate=0.*I3Units.m,
    Energy=1.*I3Units.GeV, FlasherWidth=127, FlasherBrightness=127, Seed=12345, NEvents=100,
    IceModel='spice_mie', DisableTilt=False, Filename=""):
    
    """

    :param PhotonSource: the type of photon source (CASCADE or FLASHER)
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
    
    """

    # check sanity of args
    if PhotonSource.upper() != "CASCADE" and PhotonSource.upper() != "FLASHER":
        print("photon source %s is unknown. Please specify either CASCADE or FLASHER!" % PhotonSource)
        sys.exit(1)
    
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

    ptype = I3Particle.ParticleType.EMinus

    def reference_source():
        source = I3Particle()
        source.type = ptype
        source.energy = Energy
        source.pos = I3Position(0., 0., ZCoordinate)
        source.dir = I3Direction(Zenith, Azimuth)
        source.time = 0.
        source.length = 0.
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
            source = self.reference_source()
            primary = I3Particle()
            # primary.id = I3Particle().id
            # primary.type = I3Particle.NuE
            # primary.location_type = primary.Anywhere
            
            mctree = I3MCTree()
            mctree.add_primary(primary)
            mctree.append_child(primary, source)
    
            # clsim likes I3MCTrees
            frame["I3MCTree"] = mctree
            
            self.PushFrame(frame)
            
            self.emittedEvents += 1
            if self.emittedEvents >= self.nevents:
                self.RequestSuspension()

    tray.AddModule(MakeParticle, SourceFunction=reference_source, NEvents=NEvents)

    if PhotonSource.upper() == "FLASHER":
        flasherpulse = "I3FlasherPulseSeriesMap"
        mctree = None
    elif PhotonSource.upper() == "CASCADE":
        flasherpulse = None
        mctree = "I3MCTree"
    
    header = dict(FITSTable.empty_header)
    header['zenith'] = Zenith/I3Units.degree
    header['z'] = ZCoordinate
    header['energy'] = Energy
    header['type'] = int(ptype)
    header['efficiency'] = Efficiency.RECEIVER | Efficiency.WAVELENGTH
    if PhotonSource.upper() == "FLASHER":
        header['flasherwidth'] = FlasherWidth
        header['flasherbrightness'] = FlasherBrightness
    
    tray.AddSegment(I3CLSimTabulatePhotons, name+"makeCLSimPhotons",
        PhotonSeriesName = "PropagatedPhotons",
        MCTreeName = mctree,                        # if source is a cascade this will point to the I3MCTree
        FlasherPulseSeriesName = flasherpulse,      # if source is a flasher this will point to the I3CLSimFlasherPulseSeries
        MMCTrackListName = None,                    # do NOT use MMC
        ParallelEvents = 1,                         # only work at one event at a time (it'll take long enough)
        RandomService = randomService,
        # UnWeightedPhotons=True,
        UseGPUs=False,                              # table-making is not a workload particularly suited to GPUs
        UseCPUs=True,                               # it should work fine on CPUs, though
        UnshadowedFraction=1.0,                     # no cable shadow
        StopDetectedPhotons=False,                  # do not stop photons on detection (also somewhat pointless without DOMs)
        PhotonHistoryEntries=10000,                 # record all photon paths
        DoNotParallelize=True,                      # no multithreading
        UseGeant4=False,
        OverrideApproximateNumberOfWorkItems=1,     # if you *would* use multi-threading, this would be the maximum number of jobs to run in parallel (OpenCL is free to split them)
        ExtraArgumentsToI3CLSimModule=dict(Filename=Filename, TableHeader=header,
            ReferenceSource=reference_source()),
        IceModelLocation=expandvars("$I3_SRC/clsim/resources/ice/" + IceModel),
        DisableTilt=DisableTilt,
    )
