import string
from os.path import expandvars, exists, isdir, isfile

from .. import icetray, dataclasses
from . import I3CLSimLightSourceParameterization
from . import I3CLSimFlasherPulse, I3CLSimLightSourceToStepConverterFlasher
from . import GetIceCubeFlasherSpectrum

def genFlasherParameterizationList(spectrumTable):
    spectrumTypes = [I3CLSimFlasherPulse.FlasherPulseType.LED340nm,
                     I3CLSimFlasherPulse.FlasherPulseType.LED370nm,
                     I3CLSimFlasherPulse.FlasherPulseType.LED405nm,
                     I3CLSimFlasherPulse.FlasherPulseType.LED450nm,
                     I3CLSimFlasherPulse.FlasherPulseType.LED505nm]
    
    parameterizations = []
    for flasherSpectrumType in spectrumTypes:
        theSpectrum = GetIceCubeFlasherSpectrum(spectrumType=flasherSpectrumType)
        theConverter = I3CLSimLightSourceToStepConverterFlasher(flasherSpectrumNoBias=theSpectrum, spectrumTable=spectrumTable)
        parameterization = I3CLSimLightSourceParameterization(converter=theConverter, forFlasherPulseType=flasherSpectrumType)
        parameterizations.append(parameterization)
        
    return parameterizations

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
          I3CLSimLightSourceParameterization(
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
          I3CLSimLightSourceParameterization(
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
                    FlasherInfoVectName=None,
                    MMCTrackListName="MMCTrackList",
                    MCHitSeriesName="MCHitSeriesMap",
                    PhotonSeriesName=None,
                    SimulateAfterPulses=False,
                    ParallelEvents=1000,
                    RandomService=None,
                    IceModelLocation=expandvars("$I3_SRC/clsim/resources/ice/spice_mie"),
                    UnWeightedPhotons=False,
                    UseGeant4=False,
                    DoNotParallelize=False,
                    DOMOversizeFactor=5.,
                    If=lambda f: True
                    ):
    """Do standard clsim processing, compatible to hit-maker/PPC.
    Reads its particles from an I3MCTree and writes an I3MCHitSeriesMap.

    Only the PPC-compatible parameterizations are used,
    Geant4 is not.

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
    :param MCTreeName:
        The name of the I3MCTree containing the particles to propagate.
    :param FlasherPulseSeriesName:
        Set this to the name of I3FlasherInfoVect objects in the frame to
        enable flasher simulation. The module will read I3FlasherInfoVect objects
        and generate photons according to assumed parameterizations.
    :param MMCTrackListName:
        Only used if *ChopMuons* is active. Set it to the name
        of the I3MMCTrackList object that contains additional
        muon energy loss information.
    :param MCHitSeriesName:
        Name of the output I3MCHitSeriesMap written by the module.
        Set this to None to prevent generating MCHits from
        Photons.
    :param PhotonSeriesName:
        Configure this to enable writing an I3PhotonSeriesMap containing
        all photons that reached the DOM surface.
    :param SimulateAfterPulses:
        Use an algorithm from hit-maker to simulate after-pulses.
        Turn this on to be compatble to hit-maker with afer-pulse
        simulation. This also includes PMT jitter simulation (2ns).
        Do not use this when using external after-pulse simulation modules.
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
    :param DOMOversizeFactor:
        Set the DOM oversize factor. To disable oversizing, set this to 1.
    :param If:
        Python function to use as conditional execution test for segment modules.        
    """

    from icecube import icetray, dataclasses, clsim

    # simple sanity check
    if (PhotonSeriesName is None) and (MCHitSeriesName is None):
        raise RuntimeError("You need to set at least one of the \"PhotonSeriesName\" or \"MCHitSeriesName\" arguments to something other than None!")

    # some constants
    if UnWeightedPhotons:
        # no oversizing if the user requested unweighted photons
        RadiusOverSizeFactor = 1.
    else:
        RadiusOverSizeFactor = DOMOversizeFactor
    DOMRadius = 0.16510*icetray.I3Units.m # 13" diameter
    Jitter = 2.*icetray.I3Units.ns
    numFlasherPhotonsAtFullBrightness = 8.0e9

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

    if PhotonSeriesName is not None:
        photonsName=PhotonSeriesName
    else:
        photonsName=name + "____intermediatePhotons"

    if FlasherInfoVectName is None or FlasherInfoVectName=="":
        SimulateFlashers=False
        clSimFlasherPulseSeriesName = ""
        clSimOMKeyMaskName = ""
    else:
        SimulateFlashers=True
        clSimFlasherPulseSeriesName = FlasherInfoVectName + "_pulses"
        clSimOMKeyMaskName = FlasherInfoVectName + "_OMKeys"
        
        tray.AddModule(clsim.FlasherInfoVectToFlasherPulseSeriesConverter,
                       "FlasherInfoVectToFlasherPulseSeriesConverter",
                       FlasherInfoVectName = FlasherInfoVectName,
                       FlasherPulseSeriesName = clSimFlasherPulseSeriesName,
                       FlasherOMKeyVectName = clSimOMKeyMaskName,
                       NumberOfPhotonsAtMaxBrightness = numFlasherPhotonsAtFullBrightness,
                       If=If)

    # (optional) pre-processing
    if ChopMuons:
        tray.AddModule("I3MuonSlicer", name + "_chopMuons",
                       InputMCTreeName=clSimMCTreeName,
                       MMCTrackListName=MMCTrackListName,
                       OutputMCTreeName=clSimMCTreeName + "_sliced",
                       If=If)
        clSimMCTreeName = clSimMCTreeName + "_sliced"
    else:
        clSimMCTreeName = clSimMCTreeName

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
    ppcConverter = clsim.I3CLSimLightSourceToStepConverterPPC(photonsPerStep=200)
    if not UseGeant4:
        particleParameterizations = clsim.genDefaultParameterizationList(ppcConverter, muonOnly=False)
    else:
        if MMCTrackListName is None or MMCTrackListName=="":
            particleParameterizations = [] # make sure absolutely **no** parameterizations are used
        else:
            # use no parameterizations except for muons with lengths assigned to them
            # (those are assumed to have been generated by MMC)
            particleParameterizations = clsim.genDefaultParameterizationList(ppcConverter, muonOnly=True)

    # flasher parameterizations
    if SimulateFlashers:
        # this needs a spectrum table in order to pass spectra to OpenCL
        spectrumTable = clsim.I3CLSimSpectrumTable()
        particleParameterizations += clsim.genFlasherParameterizationList(spectrumTable)
        
        print "number of spectra (1x Cherenkov + Nx flasher):", len(spectrumTable)
    else:
        # no spectrum table is necessary when only using the Cherenkov spectrum
        spectrumTable = None

    # after-pulse simulation
    has_I3CLSimPMTPhotonSimulatorIceCube = "I3CLSimPMTPhotonSimulatorIceCube" in clsim.__dict__
    if SimulateAfterPulses:
        if not has_I3CLSimPMTPhotonSimulatorIceCube:
            raise RuntimeError("cannot simulate jitter/after-pulses because hit-maker is not installed")
        pmtPhotonSimulator = clsim.I3CLSimPMTPhotonSimulatorIceCube(jitter=Jitter)
    else:
        if has_I3CLSimPMTPhotonSimulatorIceCube:
            pmtPhotonSimulator = clsim.I3CLSimPMTPhotonSimulatorIceCube(jitter=0.,
                                                                        pre_pulse_probability=0.,
                                                                        late_pulse_probability=0.,
                                                                        after_pulse_probability=0.)
        else:
            pmtPhotonSimulator = None

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
                   SpectrumTable=spectrumTable,
                   FlasherPulseSeriesName=clSimFlasherPulseSeriesName,
                   OMKeyMaskName=clSimOMKeyMaskName,
                   # ignore IceTop
                   IgnoreSubdetectors = ["IceTop"],
                   IgnoreNonIceCubeOMNumbers=False, 
                   GenerateCherenkovPhotonsWithoutDispersion=False,
                   WavelengthGenerationBias=wavelengthGenerationBias,
                   ParameterizationList=particleParameterizations,
                   MaxNumParallelEvents=ParallelEvents,
                   OpenCLDeviceList=openCLDevices,
                   UseHardcodedDeepCoreSubdetector=False, # setting this to true saves GPU constant memory but will reduce performance
                   If=If
                   )

    if MCHitSeriesName is not None:
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
                       IgnoreDOMsWithoutDetectorStatusEntry = True,
                       If=If)

    if PhotonSeriesName is None:
        tray.AddModule("Delete", "delete_photons",
            Keys = [photonsName])


