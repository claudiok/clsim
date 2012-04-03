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
# @file I3CLSimMakeHitsFromPhotons.py
# @version $Revision$
# @date $Date$
# @author Claudio Kopper
#

import string
from os.path import expandvars, exists, isdir, isfile

from icecube import icetray, dataclasses

# use this instead of a simple "@icetray.traysegment" to support
# ancient versions of IceTray that do not have tray segments.
def unchanged(func): return func
my_traysegment = icetray.traysegment if hasattr(icetray, "traysegment") else unchanged
@my_traysegment
def I3CLSimMakeHitsFromPhotons(tray, name,
                               MCTreeName="I3MCTree_sliced",
                               PhotonSeriesName="PhotonSeriesMap",
                               MCHitSeriesName="MCHitSeriesMap",
                               SimulateAfterPulses=False,
                               RandomService=None,
                               DOMOversizeFactor=5.,
                               UnshadowedFraction=0.9,
                               UseHoleIceParameterization=True,
                               If=lambda f: True
                               ):
    """
    Convert I3Photons into I3MCHits. This applies the DOM
    angular acceptance (and wavenelgth acceptance in case
    you are using the unbiased photon propagation mode.)

    :param MCTreeName:
        The name of the I3MCTree containing the particles to propagate.
    :param PhotonSeriesName:
        Name of the input I3PhotonSeriesMap to be converted.
    :param MCHitSeriesName:
        Name of the output I3MCHitSeriesMap written by the module.
        Set this to None to prevent generating MCHits from
        Photons.
    :param SimulateAfterPulses:
        Use an algorithm from hit-maker to simulate after-pulses.
        Turn this on to be compatble to hit-maker with afer-pulse
        simulation. This also includes PMT jitter simulation (2ns).
        Do not use this when using external after-pulse simulation modules.
    :param RandomService:
        Set this to an instance of a I3RandomService. Alternatively,
        you can specify the name of a configured I3RandomServiceFactory
        added to I3Tray using tray.AddService(). If you don't configure
        this, the default I3RandomServiceFactory will be used.
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

    # some constants
    DOMRadius = 0.16510*icetray.I3Units.m # 13" diameter
    Jitter = 2.*icetray.I3Units.ns

    # detector properties
    domAcceptance = clsim.GetIceCubeDOMAcceptance(domRadius = DOMRadius*DOMOversizeFactor, efficiency=UnshadowedFraction)
    domAngularSensitivity = clsim.GetIceCubeDOMAngularSensitivity(holeIce=UseHoleIceParameterization)

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

    tray.AddModule("I3PhotonToMCHitConverter", name + "_clsim_make_hits",
                   RandomService = RandomService,
                   MCTreeName = MCTreeName,
                   InputPhotonSeriesMapName = PhotonSeriesName,
                   OutputMCHitSeriesMapName = MCHitSeriesName,
                   DOMRadiusWithoutOversize=DOMRadius,
                   DOMOversizeFactor = DOMOversizeFactor,
                   WavelengthAcceptance = domAcceptance,
                   AngularAcceptance = domAngularSensitivity,
                   PMTPhotonSimulator = pmtPhotonSimulator,
                   IgnoreDOMsWithoutDetectorStatusEntry = True,
                   If=If)

