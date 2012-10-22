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
# @file GetFlasherParameterizationList.py
# @version $Revision$
# @date $Date$
# @author Claudio Kopper
#

from icecube import icetray, dataclasses
from icecube.clsim import I3CLSimLightSourceParameterization
from icecube.clsim import I3CLSimFlasherPulse, I3CLSimLightSourceToStepConverterPointSource
from icecube.clsim import GetIceCubeFlasherSpectrum

def GetPointSourceParameterizationList():
    spectrumTypes = [I3CLSimFlasherPulse.FlasherPulseType.Unknown]
    
    parameterizations = []
    for flasherSpectrumType in spectrumTypes:
        theConverter = I3CLSimLightSourceToStepConverterPointSource()
        parameterization = I3CLSimLightSourceParameterization(converter=theConverter, forFlasherPulseType=flasherSpectrumType)
        parameterizations.append(parameterization)
        
    return parameterizations