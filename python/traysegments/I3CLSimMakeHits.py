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
# $Id: I3CLSimMakeHits.py 130869 2015-04-01 21:37:59Z benedikt.riedel $
# 
# @file I3CLSimMakeHits.py
# @version $Revision: 130869 $
# @date $Date: 2015-04-01 17:37:59 -0400 (Wed, 01 Apr 2015) $
# @author Claudio Kopper
#

from __future__ import print_function

import string
from os.path import expandvars, exists, isdir, isfile

from icecube import icetray, dataclasses

from icecube.clsim import GetDefaultParameterizationList
from icecube.clsim import GetFlasherParameterizationList
from icecube.clsim import GetHybridParameterizationList

from .I3CLSimMakePhotons import I3CLSimMakePhotons
from .I3CLSimMakeHitsFromPhotons import I3CLSimMakeHitsFromPhotons


# use this instead of a simple "@icetray.traysegment" to support
# ancient versions of IceTray that do not have tray segments.
def unchanged(func): return func
my_traysegment = icetray.traysegment if hasattr(icetray, "traysegment") else unchanged
@my_traysegment
def I3CLSimMakeHits(tray, name,
                    UseCPUs=False,
                    UseGPUs=True,
                    UseOnlyDeviceNumber=None,
                    MCTreeName="I3MCTree",
                    OutputMCTreeName=None,
                    FlasherInfoVectName=None,
                    FlasherPulseSeriesName=None,
                    MMCTrackListName="MMCTrackList",
                    MCPESeriesName="MCPESeriesMap",
                    PhotonSeriesName=None,
                    ParallelEvents=1000,
                    TotalEnergyToProcess=0.,
                    RandomService=None,
                    IceModelLocation=expandvars("$I3_SRC/clsim/resources/ice/spice_mie"),
                    DisableTilt=False,
                    UnWeightedPhotons=False,
                    UseGeant4=False,
                    CrossoverEnergyEM=None,
                    CrossoverEnergyHadron=None,
                    StopDetectedPhotons=True,
                    PhotonHistoryEntries=0,
                    DoNotParallelize=False,
                    DOMOversizeFactor=5.,
                    UnshadowedFraction=0.9,
                    UseHoleIceParameterization=True,
                    ExtraArgumentsToI3CLSimModule=dict(),
                    If=lambda f: True
                    ):
    """Do standard clsim processing, compatible to hit-maker/PPC.
    Reads its particles from an I3MCTree and writes an I3MCPESeriesMap.

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
    :param FlasherInfoVectName:
        Set this to the name of an I3FlasherInfoVect object in the frame to
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
    :param MCPESeriesName:
        Name of the output I3MCPESeriesMap written by the module.
        Set this to None to prevent generating MCPEs from
        Photons.
    :param PhotonSeriesName:
        Configure this to enable writing an I3PhotonSeriesMap containing
        all photons that reached the DOM surface.
    :param ParallelEvents:
        clsim will work on a couple of events in parallel in order
        not to starve the GPU. Setting this too high will result
        in excessive memory usage (all your frames have to be cached
        in RAM). Setting it too low may impact simulation performance.
        The optimal value depends on your energy distribution/particle type.
    :param TotalEnergyToProcess:
       clsim will work on a couple of events in parallel in order
       not to starve the GPU. With this setting clsim will figure out
       how many frames to accumulate as to not starve the GPU based on 
       the energy of the light sources that are producing the photons 
       in the detector. Setting this too high will result
       in excessive memory usage (all your frames have to be cached
       in RAM). Setting it too low may impact simulation performance. 
       This cannot be used in flasher mode, since we cannot measure
       the energy of the light sources.
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
        in parallel on a batch system.
    :param DOMOversizeFactor:
        Set the DOM oversize factor. To disable oversizing, set this to 1.
    :param UnshadowedFraction:
        Fraction of photocathode available to receive light (e.g. unshadowed by the cable)
    :param UseHoleIceParameterization:
        Use an angular acceptance correction for hole ice scattering.
    :param If:
        Python function to use as conditional execution test for segment modules.        
    """

    from icecube import icetray, dataclasses, clsim

    # simple sanity check
    if (PhotonSeriesName is None) and (MCPESeriesName is None):
        raise RuntimeError("You need to set at least one of the \"PhotonSeriesName\" or \"MCPESeriesName\" arguments to something other than None!")

    # warn the user in case they might have done something they probably don't want
    if UnWeightedPhotons and (DOMOversizeFactor != 1.):
        print("********************")
        print("Enabling the clsim.I3CLSimMakeHits() \"UnWeightedPhotons=True\" option without setting")
        print("\"DOMOversizeFactor=1.\" will still apply a constant weighting factor of DOMOversizeFactor**2.")
        print("If this is what you want, you can safely ignore this warning.")
        print("********************")

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
        if OutputMCTreeName is not None:
            raise RuntimeError("cannot have an output MCTree without an input MCTree")

        clSimMCTreeName=""
        if ChopMuons:
            raise RuntimeError("You cannot have \"MMCTrackListName\" enabled with no MCTree!")
    else:
        if OutputMCTreeName is None or OutputMCTreeName=="":
            if ChopMuons:
                clSimMCTreeName=MCTreeName+"_sliced"
            else:
                clSimMCTreeName=MCTreeName+"_clsim"
        else:
            clSimMCTreeName = OutputMCTreeName

    I3CLSimMakePhotons_kwargs = dict(UseCPUs=UseCPUs,
                                     UseGPUs=UseGPUs,
                                     UseOnlyDeviceNumber=UseOnlyDeviceNumber,
                                     MCTreeName=MCTreeName,
                                     OutputMCTreeName=clSimMCTreeName,
                                     FlasherInfoVectName=FlasherInfoVectName,
                                     FlasherPulseSeriesName=FlasherPulseSeriesName,
                                     MMCTrackListName=MMCTrackListName,
                                     PhotonSeriesName=photonsName,
                                     ParallelEvents=ParallelEvents,
                                     TotalEnergyToProcess=TotalEnergyToProcess,
                                     RandomService=RandomService,
                                     IceModelLocation=IceModelLocation,
                                     DisableTilt=DisableTilt,
                                     UnWeightedPhotons=UnWeightedPhotons,
                                     UseGeant4=UseGeant4,
                                     CrossoverEnergyEM=CrossoverEnergyEM,
                                     CrossoverEnergyHadron=CrossoverEnergyHadron,
                                     StopDetectedPhotons=StopDetectedPhotons,
                                     PhotonHistoryEntries=PhotonHistoryEntries,
                                     DoNotParallelize=DoNotParallelize,
                                     DOMOversizeFactor=DOMOversizeFactor,
                                     UnshadowedFraction=UnshadowedFraction,
                                     UseHoleIceParameterization=UseHoleIceParameterization,
                                     ExtraArgumentsToI3CLSimModule=ExtraArgumentsToI3CLSimModule,
                                     If=If)

    if hasattr(icetray, "traysegment"):
        tray.AddSegment(I3CLSimMakePhotons, name + "_makePhotons",
                        **I3CLSimMakePhotons_kwargs)
    else:
        # if there is no tray segment support in IceTray, use this instead:
        I3CLSimMakePhotons(tray, name + "_makePhotons",
                           **I3CLSimMakePhotons_kwargs)

    if MCPESeriesName is not None:
        I3CLSimMakeHitsFromPhotons_kwargs = dict(MCTreeName=clSimMCTreeName,
                                                 PhotonSeriesName=photonsName,
                                                 MCPESeriesName=MCPESeriesName,
                                                 RandomService=RandomService,
                                                 DOMOversizeFactor=DOMOversizeFactor,
                                                 UnshadowedFraction=UnshadowedFraction,
                                                 UseHoleIceParameterization=UseHoleIceParameterization,
                                                 If=If)
        
        if hasattr(icetray, "traysegment"):
            tray.AddSegment(I3CLSimMakeHitsFromPhotons, name + "_makeHitsFromPhotons",
                            **I3CLSimMakeHitsFromPhotons_kwargs
                            )
        else:
            # if there is no tray segment support in IceTray, use this instead:
            I3CLSimMakeHitsFromPhotons(tray, name + "_makeHitsFromPhotons",
                            **I3CLSimMakeHitsFromPhotons_kwargs
                            )
            
    if PhotonSeriesName is None:
        tray.AddModule("Delete", name + "_deletePhotons",
            Keys = [photonsName],
            If=If)


