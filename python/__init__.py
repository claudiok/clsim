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
# @file __init__.py
# @version $Revision$
# @date $Date$
# @author Claudio Kopper
#

from icecube.load_pybindings import load_pybindings
from icecube import icetray, dataclasses # be nice and pull in our dependencies
load_pybindings(__name__,__path__)


from MakeAntaresMediumProperties import GetPetzoldScatteringCosAngleDistribution, GetAntaresScatteringCosAngleDistribution, MakeAntaresMediumProperties
from MakeIceCubeMediumProperties import MakeIceCubeMediumProperties
from MakeIceCubeMediumPropertiesPhotonics import MakeIceCubeMediumPropertiesPhotonics
from GetIceCubeDOMAcceptance import GetIceCubeDOMAcceptance
from GetIceCubeDOMAngularSensitivity import GetIceCubeDOMAngularSensitivity
from GetIceCubeFlasherSpectrum import GetIceCubeFlasherSpectrum

from GetAntaresOMAcceptance import GetAntaresOMAcceptance
from GetAntaresOMAngularSensitivity import GetAntaresOMAngularSensitivity

from GetKM3NeTDOMAcceptance import GetKM3NeTDOMAcceptance

from FlasherInfoVectToFlasherPulseSeriesConverter import FlasherInfoVectToFlasherPulseSeriesConverter
from FakeFlasherInfoGenerator import FakeFlasherInfoGenerator

from GetDefaultParameterizationList import GetDefaultParameterizationList
from GetFlasherParameterizationList import GetFlasherParameterizationList
from AsyncTap import AsyncTap
from AutoSetGeant4Environment import AutoSetGeant4Environment

from I3CLSimRandomValueIceCubeFlasherTimeProfile import I3CLSimRandomValueIceCubeFlasherTimeProfile

import util

# import tray segments (if available)
from traysegments import I3CLSimMakeHits, I3CLSimMakePhotons, I3CLSimMakeHitsFromPhotons

# clean up the clsim namespace
del icetray
del dataclasses

