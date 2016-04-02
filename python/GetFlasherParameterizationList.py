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
# $Id: GetFlasherParameterizationList.py 108199 2013-07-12 21:33:08Z nwhitehorn $
# 
# @file GetFlasherParameterizationList.py
# @version $Revision: 108199 $
# @date $Date: 2013-07-12 17:33:08 -0400 (Fri, 12 Jul 2013) $
# @author Claudio Kopper
#

from __future__ import print_function

from icecube import icetray, dataclasses
from icecube.clsim import I3CLSimLightSourceParameterization
from icecube.clsim import I3CLSimFlasherPulse, I3CLSimLightSourceToStepConverterFlasher
from icecube.clsim import GetIceCubeFlasherSpectrum
from icecube.clsim import I3CLSimRandomValueNormalDistribution
from icecube.clsim import I3CLSimRandomValueFixParameter
from icecube.clsim import I3CLSimRandomValueConstant
from icecube.clsim import I3CLSimRandomValueUniform

from . import I3CLSimRandomValueIceCubeFlasherTimeProfile

from I3Tray import I3Units

# for now, all flasher types get the same time delay profile
__theFlasherTimeDelayDistribution = I3CLSimRandomValueIceCubeFlasherTimeProfile.I3CLSimRandomValueIceCubeFlasherTimeProfile()

def GetFlasherParameterizationList(spectrumTable):
    spectrumTypes = [I3CLSimFlasherPulse.FlasherPulseType.LED340nm,
                     I3CLSimFlasherPulse.FlasherPulseType.LED370nm,
                     I3CLSimFlasherPulse.FlasherPulseType.LED405nm,
                     I3CLSimFlasherPulse.FlasherPulseType.LED450nm,
                     I3CLSimFlasherPulse.FlasherPulseType.LED505nm]

    spectrumTypesSC = [I3CLSimFlasherPulse.FlasherPulseType.SC1,
                       I3CLSimFlasherPulse.FlasherPulseType.SC2]
    
    # all flasher types have the same angular smearing profiles (a gaussian
    # with its width set as a runtime parameter [read from I3CLSimFlasherPulse])
    normalDistribution = I3CLSimRandomValueFixParameter(I3CLSimRandomValueNormalDistribution(), 0, 0.) # the mean (parameter #0) is fixed to 0.
    
    # Standard Candle angular distributions
    # (mean[parameter #1] fixed to 2ns)
    standardCandleTimeDelayDistribution =  I3CLSimRandomValueFixParameter(I3CLSimRandomValueNormalDistribution(), 0, 2.*I3Units.ns)
    standardCandlePolarDistribution = I3CLSimRandomValueConstant()
    standardCandleAzimuthalDistribution =  I3CLSimRandomValueUniform(0., float('NaN')) # make the "to" parameter a run-time setting
    
    
    # generate the parameterizations
    parameterizations = []
    
    # (for flashers)
    for flasherSpectrumType in spectrumTypes:
        theSpectrum = GetIceCubeFlasherSpectrum(spectrumType=flasherSpectrumType)
        theConverter = I3CLSimLightSourceToStepConverterFlasher(flasherSpectrumNoBias=theSpectrum,
                                                                spectrumTable=spectrumTable,
                                                                angularProfileDistributionPolar=normalDistribution,
                                                                angularProfileDistributionAzimuthal=normalDistribution,
                                                                timeDelayDistribution=__theFlasherTimeDelayDistribution,
                                                                interpretAngularDistributionsInPolarCoordinates=False)
        parameterization = I3CLSimLightSourceParameterization(converter=theConverter, forFlasherPulseType=flasherSpectrumType)
        parameterizations.append(parameterization)

    # for the standard candles
    for flasherSpectrumType in spectrumTypesSC:
        theSpectrum = GetIceCubeFlasherSpectrum(spectrumType=flasherSpectrumType)
        theConverter = I3CLSimLightSourceToStepConverterFlasher(flasherSpectrumNoBias=theSpectrum,
                                                                spectrumTable=spectrumTable,
                                                                angularProfileDistributionPolar=standardCandlePolarDistribution,
                                                                angularProfileDistributionAzimuthal=standardCandleAzimuthalDistribution,
                                                                timeDelayDistribution=standardCandleTimeDelayDistribution,
                                                                interpretAngularDistributionsInPolarCoordinates=True)
        parameterization = I3CLSimLightSourceParameterization(converter=theConverter, forFlasherPulseType=flasherSpectrumType)
        parameterizations.append(parameterization)
        
    return parameterizations
