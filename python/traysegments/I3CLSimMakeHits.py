import string
from os.path import expandvars, exists, isdir, isfile

from icecube import icetray, dataclasses

from icecube.clsim import GetDefaultParameterizationList
from icecube.clsim import GetFlasherParameterizationList

from I3CLSimMakePhotons import I3CLSimMakePhotons
from I3CLSimMakeHitsFromPhotons import I3CLSimMakeHitsFromPhotons

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
                    StopDetectedPhotons=True,
                    PhotonHistoryEntries=0,
                    DoNotParallelize=False,
                    DOMOversizeFactor=5.,
                    If=lambda f: True
                    ):
    """Do standard clsim processing, compatible to hit-maker/PPC.
    Reads its particles from an I3MCTree and writes an I3MCHitSeriesMap.

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

    # warn the user in case they might have done something they probably don't want
    if UnWeightedPhotons and (DOMOversizeFactor != 1.):
        print "********************"
        print "Enabling the clsim.I3CLSimMakeHits() \"UnWeightedPhotons=True\" option without setting"
        print "\"DOMOversizeFactor=1.\" will still apply a constant weighting factor of DOMOversizeFactor**2."
        print "If this is what you want, you can safely ignore this warning."
        print "********************"

    if PhotonSeriesName is not None:
        photonsName=PhotonSeriesName
    else:
        photonsName=name + "____intermediatePhotons"

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
        if ChopMuons:
            clSimMCTreeName=MCTreeName+"_sliced"
        else:
            clSimMCTreeName=MCTreeName+"_clsim"


    tray.AddSegment(I3CLSimMakePhotons, name + "_makePhotons",
                    UseCPUs=UseCPUs,
                    UseGPUs=UseGPUs,
                    MCTreeName=MCTreeName,
                    OutputMCTreeName=clSimMCTreeName,
                    FlasherInfoVectName=FlasherInfoVectName,
                    MMCTrackListName=MMCTrackListName,
                    PhotonSeriesName=photonsName,
                    ParallelEvents=ParallelEvents,
                    RandomService=RandomService,
                    IceModelLocation=IceModelLocation,
                    UnWeightedPhotons=UnWeightedPhotons,
                    UseGeant4=UseGeant4,
                    StopDetectedPhotons=StopDetectedPhotons,
                    PhotonHistoryEntries=PhotonHistoryEntries,
                    DoNotParallelize=DoNotParallelize,
                    DOMOversizeFactor=DOMOversizeFactor,
                    If=If
                    )

    if MCHitSeriesName is not None:
        tray.AddSegment(I3CLSimMakeHitsFromPhotons, name + "_makeHitsFromPhotons",
                        MCTreeName=clSimMCTreeName,
                        PhotonSeriesName=photonsName,
                        MCHitSeriesName=MCHitSeriesName,
                        SimulateAfterPulses=SimulateAfterPulses,
                        RandomService=RandomService,
                        DOMOversizeFactor=DOMOversizeFactor,
                        If=If
                        )

    if PhotonSeriesName is None:
        tray.AddModule("Delete", name + "_deletePhotons",
            Keys = [photonsName],
            If=If)


