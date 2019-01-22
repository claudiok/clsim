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
# $Id: I3CLSimMakeHitsFromPhotons.py 108199 2013-07-12 21:33:08Z nwhitehorn $
# 
# @file I3CLSimMakeHitsFromPhotons.py
# @version $Revision: 108199 $
# @date $Date: 2013-07-12 17:33:08 -0400 (Fri, 12 Jul 2013) $
# @author Claudio Kopper
#

import string
from os.path import expandvars, exists, isdir, isfile

from icecube import icetray, dataclasses, dataio

# use this instead of a simple "@icetray.traysegment" to support
# ancient versions of IceTray that do not have tray segments.
def unchanged(func): return func

def get_compensation_factor(gcd_file):
    fr = gcd_file.pop_frame()
    while "I3Calibration" not in fr :
        fr = gcd_file.pop_frame()
    domcal = fr.Get("I3Calibration").dom_cal
    max_compensation_factor = 0
    for om, i3domcal in domcal.items():
        if not hasattr(i3domcal,"combined_spe_charge_distribution"):
            logging.log_error("I3Calibration has no element 'combined_spe_charge_distribution'!!!")
            logging.log_warn("Setting max_compensation_factor to 1.0")
            max_compensation_factor = 1.0
            break
        max_compensation_factor = max(max_compensation_factor,
            i3domcal.combined_spe_charge_distribution.compensation_factor)
    icetray.logging.log_info("compensation factor %s" % max_compensation_factor)
    return max_compensation_factor

my_traysegment = icetray.traysegment if hasattr(icetray, "traysegment") else unchanged
@my_traysegment
def I3CLSimMakeHitsFromPhotons(tray, name,
                               PhotonSeriesName="PhotonSeriesMap",
                               MCPESeriesName="MCPESeriesMap",
                               RandomService=None,
                               DOMOversizeFactor=5.,
                               UnshadowedFraction=1.0,
                               IceModelLocation=None, #Needed for icemodel-dependent efficiency
                               HoleIceParameterization=expandvars("$I3_BUILD/ice-models/resources/models/angsens/as.h2-50cm"),
                               MergeHits=False,
                               GCDFile= None,
                               If=lambda f: True
                               ):
    """
    Convert I3Photons into I3MCPEs. This applies the DOM
    angular acceptance (and wavenelgth acceptance in case
    you are using the unbiased photon propagation mode.)

    :param PhotonSeriesName:
        Name of the input I3PhotonSeriesMap to be converted.
    :param MCPESeriesName:
        Name of the output I3MCPESeriesMap written by the module.
        Set this to None to prevent generating MCPEs from
        Photons.
    :param RandomService:
        Set this to an instance of a I3RandomService. Alternatively,
        you can specify the name of a configured I3RandomServiceFactory
        added to I3Tray using tray.AddService(). If you don't configure
        this, the default I3RandomServiceFactory will be used.
    :param DOMOversizeFactor:
        Set the DOM oversize factor. To disable oversizing, set this to 1.
    :param UnshadowedFraction:
        Fraction of photocathode available to receive light (e.g. unshadowed by the cable)
    :param HoleIceParameterization:
        Set this to a hole ice parameterization file. The default file contains the 
        coefficients for nominal angular acceptance correction due to hole ice (ice-models 
        project is required). Use file $I3_BUILD/ice-models/resources/models/angsens/as.nominal 
        for no hole ice parameterization. 
    :param MergeHits:
    	Set to true to perform time merging on the MCPE as they are produced. This is useful for 
    	reducing the memory and disk space used by high energy (bright) events, and can allow 
    	detector simulation to run much more quickly. This causes parent particle information to be 
    	stored in an additional frame object. 
    :param If:
        Python function to use as conditional execution test for segment modules.        
    """

    from icecube import icetray, dataclasses, clsim
    from icecube.clsim.traysegments.common import parseIceModel
    
    icemodel_efficiency_factor = 1.0
    # ice properties
    if isinstance(IceModelLocation, str):
        mediumProperties = parseIceModel(IceModelLocation)
    else:
        # get ice model directly if not a string
        mediumProperties = IceModelLocation


    if mediumProperties is not None:
        # IceModel-dependent efficiency
        icemodel_efficiency_factor = mediumProperties.efficiency 
    else:
        icetray.logging.log_warn("No IceModel configured. Using DOMefficiency of %f" % icemodel_efficiency_factor)
        

    if UnshadowedFraction<=0:
        raise RuntimeError("UnshadowedFraction must be a positive number")

    # some constants
    DOMRadius = 0.16510*icetray.I3Units.m # 13" diameter
    Jitter = 2.*icetray.I3Units.ns

    # detector properties
   
    max_compensation_factor = 1.0
    if GCDFile is None:
        icetray.logging.log_warn("No GCD file given. Setting compensation factor to 1.0!!!!")
    else:
        max_compensation_factor = get_compensation_factor(dataio.I3File(GCDFile))


    icetray.logging.log_info("Net DOM efficiency (with compensation factor): %s" % (UnshadowedFraction*icemodel_efficiency_factor*max_compensation_factor))

    domAcceptance = clsim.GetIceCubeDOMAcceptance(
                                domRadius = DOMRadius*DOMOversizeFactor, 
                                efficiency=icemodel_efficiency_factor*UnshadowedFraction)
    domAngularSensitivity = clsim.GetIceCubeDOMAngularSensitivity(holeIce=HoleIceParameterization)

    tray.AddModule("I3PhotonToMCPEConverter", name + "_clsim_make_hits",
                   RandomService = RandomService,
                   InputPhotonSeriesMapName = PhotonSeriesName,
                   OutputMCPESeriesMapName = MCPESeriesName,
                   DOMRadiusWithoutOversize=DOMRadius,
                   DOMOversizeFactor = DOMOversizeFactor,
                   DOMPancakeFactor = DOMOversizeFactor,
                   WavelengthAcceptance = domAcceptance,
                   AngularAcceptance = domAngularSensitivity,
                   IgnoreDOMsWithoutDetectorStatusEntry = False, # in icesim4 it is the job of the DOM simulation tools to cut out these DOMs
                   MergeHits=MergeHits,
                   If=If)

