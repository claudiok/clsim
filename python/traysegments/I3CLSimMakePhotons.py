#
# Copyright (c) 2011, 2012
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
# @file I3CLSimMakePhotons.py
# @version $Revision$
# @date $Date$
# @author Claudio Kopper
#

from __future__ import print_function

import string
from os.path import expandvars, exists, isdir, isfile

from icecube import icetray, dataclasses

from icecube.clsim import GetDefaultParameterizationList
from icecube.clsim import GetFlasherParameterizationList

from icecube.clsim import AutoSetGeant4Environment

from icecube.clsim.traysegments.common import configureOpenCLDevices, parseIceModel

# use this instead of a simple "@icetray.traysegment" to support
# ancient versions of IceTray that do not have tray segments.
def unchanged(func): return func
my_traysegment = icetray.traysegment if hasattr(icetray, "traysegment") else unchanged
@my_traysegment
def I3CLSimMakePhotons(tray, name,
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
                       UnWeightedPhotons=False,
                       UseGeant4=False,
                       StopDetectedPhotons=True,
                       PhotonHistoryEntries=0,
                       DoNotParallelize=False,
                       DOMOversizeFactor=5.,
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

    # warn the user in case they might have done something they probably don't want
    if UnWeightedPhotons and (DOMOversizeFactor != 1.):
        print("********************")
        print("Enabling the clsim.I3CLSimMakeHits() \"UnWeightedPhotons=True\" option without setting")
        print("\"DOMOversizeFactor=1.\" will still apply a constant weighting factor of DOMOversizeFactor**2.")
        print("If this is what you want, you can safely ignore this warning.")
        print("********************")

    # some constants
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
                       NumberOfPhotonsAtMaxBrightness = numFlasherPhotonsAtFullBrightness,
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

    # detector properties
    if UseHoleIceParameterization:
        # the hole ice acceptance curve peaks at 0.75 instead of 1
        domEfficiencyCorrection = UnshadowedFraction*0.75*1.35 * 1.01 # DeepCore DOMs have a relative efficiency of 1.35 plus security margin of +1%
    else:
        domEfficiencyCorrection = UnshadowedFraction*1.35      * 1.01 # security margin of +1%
    domAcceptance = clsim.GetIceCubeDOMAcceptance(domRadius = DOMRadius*DOMOversizeFactor, efficiency=domEfficiencyCorrection)

    # photon generation wavelength bias
    #if isinstance(UnWeightedPhotons, float) or isinstance(UnWeightedPhotons, int):
    #    print("***** running unweighted simulation with a photon pre-scaling of", UnWeightedPhotons)
    #    wavelengthGenerationBias = clsim.I3CLSimFunctionConstant(UnWeightedPhotons)
    #else:
    if not UnWeightedPhotons:
        wavelengthGenerationBias = domAcceptance
    else:
        wavelengthGenerationBias = None

    # muon&cascade parameterizations
    ppcConverter = clsim.I3CLSimLightSourceToStepConverterPPC(photonsPerStep=200)
    if not UseGeant4:
        particleParameterizations = GetDefaultParameterizationList(ppcConverter, muonOnly=False)
    else:
        if MMCTrackListName is None or MMCTrackListName=="":
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

    tray.AddModule("I3CLSimModule", name + "_clsim",
                   MCTreeName=clSimMCTreeName,
                   PhotonSeriesMapName=PhotonSeriesName,
                   DOMRadius = DOMRadius,
                   DOMOversizeFactor = DOMOversizeFactor,
                   DOMPancakeFactor = DOMOversizeFactor, # you will probably want this to be the same as DOMOversizeFactor
                   RandomService=RandomService,
                   MediumProperties=mediumProperties,
                   SpectrumTable=spectrumTable,
                   FlasherPulseSeriesName=clSimFlasherPulseSeriesName,
                   OMKeyMaskName=clSimOMKeyMaskName,
                   # ignore IceTop
                   IgnoreSubdetectors = ["IceTop"],
                   #IgnoreNonIceCubeOMNumbers=False,
                   GenerateCherenkovPhotonsWithoutDispersion=False,
                   WavelengthGenerationBias=wavelengthGenerationBias,
                   ParameterizationList=particleParameterizations,
                   MaxNumParallelEvents=ParallelEvents,
                   OpenCLDeviceList=openCLDevices,
                   #UseHardcodedDeepCoreSubdetector=False, # setting this to true saves GPU constant memory but will reduce performance
                   StopDetectedPhotons=StopDetectedPhotons,
                   PhotonHistoryEntries=PhotonHistoryEntries,
                   If=If,
                   **ExtraArgumentsToI3CLSimModule
                   )

