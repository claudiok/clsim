import string
from os.path import expandvars, exists, isdir, isfile

from .. import icetray, dataclasses
from . import I3CLSimParticleParameterization


def genDefaultParameterizationList(theConverter, muonOnly=False):
    fromEnergy=0.
    toEnergy=float('Inf')

    muons    = [dataclasses.I3Particle.MuMinus,
                dataclasses.I3Particle.MuPlus]

    cascades = [dataclasses.I3Particle.Neutron,
                dataclasses.I3Particle.Hadrons,
                dataclasses.I3Particle.Pi0,
                dataclasses.I3Particle.PiPlus,
                dataclasses.I3Particle.PiMinus,
                dataclasses.I3Particle.K0_Long,
                dataclasses.I3Particle.KPlus,
                dataclasses.I3Particle.KMinus,
                dataclasses.I3Particle.PPlus,
                dataclasses.I3Particle.PMinus,
                dataclasses.I3Particle.K0_Short,
                dataclasses.I3Particle.EMinus,
                dataclasses.I3Particle.EPlus,
                dataclasses.I3Particle.Gamma,
                dataclasses.I3Particle.Brems,
                dataclasses.I3Particle.DeltaE,
                dataclasses.I3Particle.PairProd,
                dataclasses.I3Particle.NuclInt]

    parameterizationsMuon = []
    for type in muons:
        converter = \
          I3CLSimParticleParameterization(
            converter=theConverter,
            forParticleType=type,
            fromEnergy=fromEnergy,
            toEnergy=toEnergy, 
            needsLength=True)
        parameterizationsMuon.append(converter)

    if muonOnly:
        return parameterizationsMuon

    parameterizationsOther = []
    for type in cascades:
        converter = \
          I3CLSimParticleParameterization(
            converter=theConverter,
            forParticleType=type,
            fromEnergy=fromEnergy,
            toEnergy=toEnergy, 
            needsLength=False)
        parameterizationsOther.append(converter)

    return parameterizationsMuon+parameterizationsOther

@icetray.traysegment
def I3CLSimMakeHits(tray, name,
                    UseCPUs=False,
                    UseGPUs=True,
                    MCTreeName="I3MCTree",
                    MMCTrackListName="MMCTrackList",
                    MCHitSeriesName="MCHitSeriesMap",
                    PhotonSeriesName=None,
                    SimulateAfterPulses=False,
                    ParallelEvents=1000,
                    RandomService=None,
                    IceModelLocation=expandvars("$I3_SRC/clsim/resources/ice/spice_mie"),
                    UnWeightedPhotons=False,
                    UseGeant4=False,
                    DoNotParallelize=False
                    ):
    """Do standard clsim processing, compatible to hit-maker/PPC.
    Reads its particles from an I3MCTree and writes an I3MCHitSeriesMap.

    Only the PPC-compatible parameterizations are used,
    Geant4 is not.

    All available OpenCL GPUs (and optionally CPUs) will
    be used. This will take over your entire machine,
    so make sure to configure your batch jobs correctly
    when using this on a cluster. 

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
    :param MCTreeName:
        The name of the I3MCTree containing the particles to propagate.
    :param MMCTrackListName:
        Only used if *ChopMuons* is active. Set it to the name
        of the I3MMCTrackList object that contains additional
        muon energy loss information.
    :param MCHitSeriesName:
        Name of the output I3MCHitSeriesMap written by the module.
    :param PhotonSeriesName:
        Configure this to enable writing an I3PhotonSeriesMap containing
        all photons that reached the DOM surface.
    :param SimulateAfterPulses:
        Use an algorithm from hit-maker to simulate after-pulses.
        Turn this on to be compatble to hit-maker with afer-pulse
        simulation. Do not use this when using external after-pulse
        simulation modules.
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
    :param UseGeant4:
        Enabling this setting will disable all cascade and muon light yield
        parameterizations. All particles will sent to Geant4 for a full
        simulation. This does **not** apply to muons that do have a length
        assigned. These are assumed to have been generated by MMC and
        their light is generated according to the usual parameterization.
    :param DoNotParallelize:
        Try only using a single work item in parallel when running the
        OpenCL simulation. This might be useful if you want to run jobs
        in parallel on a batch system.
    """

    from icecube import icetray, dataclasses, clsim

    # some constants
    if not UnWeightedPhotons:
        RadiusOverSizeFactor = 5.
    else:
        RadiusOverSizeFactor = 1.
    DOMRadius = 0.16510*icetray.I3Units.m # 13" diameter
    Jitter = 2.*icetray.I3Units.ns

    if MMCTrackListName is None or MMCTrackListName=="":
        # the input does not seem to have been processed by MMC
        ChopMuons = False
    else:
        ChopMuons = True

    if PhotonSeriesName is not None:
        photonsName=PhotonSeriesName
    else:
        photonsName=name + "____intermediatePhotons"

    # (optional) pre-processing
    if ChopMuons:
        tray.AddModule("I3MuonSlicer", name + "_chopMuons",
                       InputMCTreeName=MCTreeName,
                       MMCTrackListName=MMCTrackListName,
                       OutputMCTreeName=MCTreeName + "_sliced")
        clSimMCTreeName = MCTreeName + "_sliced"
    else:
        clSimMCTreeName = MCTreeName

    # ice properties
    # Did the user configure a diretory or a file?
    if not exists(IceModelLocation):
        raise RuntimeError("The specified ice model path \"%s\" does not exist" % IceModelLocation)
    
    if isdir(IceModelLocation):
        # it's a PPC ice description directory
        mediumProperties = clsim.MakeIceCubeMediumProperties(iceDataDirectory=IceModelLocation)
    elif isfile(IceModelLocation):
        # it's a photonics ice description file
        mediumProperties = clsim.MakeIceCubeMediumPropertiesPhotonics(tableFile=IceModelLocation)
    else:
        raise RuntimeError("The specified ice model path \"%s\" is neither a directory nor a file." % IceModelLocation)


    # detector properties
    domAcceptance = clsim.GetIceCubeDOMAcceptance(domRadius = DOMRadius*RadiusOverSizeFactor)
    domAngularSensitivity = clsim.GetIceCubeDOMAngularSensitivity(holeIce=True)

    # photon generation wavelength bias
    if not UnWeightedPhotons:
        wavelengthGenerationBias = domAcceptance
    else:
        wavelengthGenerationBias = None

    # muon&cascade parameterizations
    cascadeConverter = clsim.I3CLSimParticleToStepConverterPPC(photonsPerStep=200)
    if not UseGeant4:
        particleParameterizations = clsim.genDefaultParameterizationList(cascadeConverter, muonOnly=False)
    else:
        if MMCTrackListName is None or MMCTrackListName=="":
            particleParameterizations = [] # make sure absolutely **no** parameterizations are used
        else:
            # use no parameterizations except for muons with lengths assigned to them
            # (those are assumed to have been generated by MMC)
            particleParameterizations = clsim.genDefaultParameterizationList(cascadeConverter, muonOnly=True)

    # after-pulse simulation
    if SimulateAfterPulses:
        pmtPhotonSimulator = clsim.I3CLSimPMTPhotonSimulatorIceCube(jitter=Jitter)
    else:
        pmtPhotonSimulator = clsim.I3CLSimPMTPhotonSimulatorIceCube(jitter=Jitter,
                                                                    pre_pulse_probability=0.,
                                                                    late_pulse_probability=0.,
                                                                    after_pulse_probability=0.)

    # get OpenCL devices
    openCLDevices = [device for device in clsim.I3CLSimOpenCLDevice.GetAllDevices() if (device.gpu and UseGPUs) or (device.cpu and UseCPUs)]

    # (auto-)configure OpenCL devices
    for device in openCLDevices:
        if string.count(device.device, 'Tesla') > 0 or string.count(device.device, 'GTX') > 0:
            # assume these are "fast", all others are "slow"
            device.useNativeMath=True
            device.approximateNumberOfWorkItems=1024000
        else:
            device.useNativeMath=False
            device.approximateNumberOfWorkItems=10240

        if DoNotParallelize:
            device.approximateNumberOfWorkItems=1

    tray.AddModule("I3CLSimModule", name + "_clsim",
                   MCTreeName=clSimMCTreeName,
                   PhotonSeriesMapName=photonsName,
                   DOMRadius = DOMRadius*RadiusOverSizeFactor, 
                   RandomService=RandomService,
                   MediumProperties=mediumProperties,
                   # ignore IceTop
                   IgnoreSubdetectors = ["IceTop"],
                   IgnoreNonIceCubeOMNumbers=False, 
                   GenerateCherenkovPhotonsWithoutDispersion=False,
                   WavelengthGenerationBias=wavelengthGenerationBias,
                   ParameterizationList=particleParameterizations,
                   MaxNumParallelEvents=ParallelEvents,
                   OpenCLDeviceList=openCLDevices,
                   UseHardcodedDeepCoreSubdetector=True
                   )

    tray.AddModule("I3PhotonToMCHitConverter", name + "_clsim_make_hits",
                   RandomService = RandomService,
                   MCTreeName = clSimMCTreeName,
                   InputPhotonSeriesMapName = photonsName,
                   OutputMCHitSeriesMapName = MCHitSeriesName,
                   DOMRadiusWithoutOversize=DOMRadius,
                   DOMOversizeFactor = RadiusOverSizeFactor,
                   WavelengthAcceptance = domAcceptance,
                   AngularAcceptance = domAngularSensitivity,
                   PMTPhotonSimulator = pmtPhotonSimulator,
                   IgnoreDOMsWithoutDetectorStatusEntry = True)

    if PhotonSeriesName is None:
        tray.AddModule("Delete", "delete_photons",
            Keys = [photonsName])


