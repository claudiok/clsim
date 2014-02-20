#
# Copyright (c) 2012
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
# @file GetRefractiveIndexRange.py
# @version $Revision$
# @date $Date$
# @author Claudio Kopper
#

import numpy

def getGroupRefIndex_derivative(wavelength):
    n_inv = 1./getPhaseRefIndex(wavelength)
    y = getDispersionPhase(wavelength);
    return 1./((1.0 + y*wavelength*n_inv) * n_inv)

def __MakeGroupRefIndexFunctionFromPhaseRefIndex(phaseRefIndex):
    if not phaseRefIndex.HasDerivative():
        raise RuntimeError("Phase refractive index does not have the option to calculate the derivative. Cannot calculate group refractive index.")

    return numpy.vectorize(lambda wlen: phaseRefIndex.GetValue(wlen)/(1. + phaseRefIndex.GetDerivative(wlen)*wlen/phaseRefIndex.GetValue(wlen)))

def GetGroupRefractiveIndexRange(mediumProperties):
    """
    Returns the minimum and maximum group refractive indices as tuple
    over all ice layers for a given mediumProperties object.
    """
    
    minIndex = None
    maxIndex = None
    #atWlen = None
    
    overallMinWlen = mediumProperties.GetMinWavelength()
    overallMaxWlen = mediumProperties.GetMaxWavelength()
    
    for layer in range(mediumProperties.LayersNum):
        minWlen = overallMinWlen
        maxWlen = overallMaxWlen
        
        
        groupRefIndex = mediumProperties.GetGroupRefractiveIndexOverride(layer)
        if groupRefIndex:
            minWlen = max(minWlen, groupRefIndex.GetMinWlen())
            maxWlen = min(maxWlen, groupRefIndex.GetMaxWlen())
            groupRefIndexFunc = numpy.vectorize(lambda wlen: groupRefIndex.GetValue(wlen))
        else:
            phaseRefIndex = mediumProperties.GetPhaseRefractiveIndex(layer)
            minWlen = max(minWlen, phaseRefIndex.GetMinWlen())
            maxWlen = min(maxWlen, phaseRefIndex.GetMaxWlen())
            groupRefIndexFunc = __MakeGroupRefIndexFunctionFromPhaseRefIndex(phaseRefIndex)
        
        testWlens = numpy.linspace(minWlen, maxWlen, 10000, endpoint=True)
        evaluatedFunc = groupRefIndexFunc(testWlens)

        thisLayerMinIndex_arg = numpy.argmin(evaluatedFunc)
        thisLayerMinIndex = evaluatedFunc[thisLayerMinIndex_arg]
        thisLayerMaxIndex_arg = numpy.argmax(evaluatedFunc)
        thisLayerMaxIndex = evaluatedFunc[thisLayerMaxIndex_arg]
        #thisLayerAtWlen = testWlens[thisLayerMaxIndex_arg]

        if (minIndex is None) or (thisLayerMinIndex < minIndex):
            minIndex = thisLayerMinIndex
        if (maxIndex is None) or (thisLayerMaxIndex > maxIndex):
            maxIndex = thisLayerMaxIndex
            #atWlen = thisLayerAtWlen
    
    #return (maxIndex, atWlen, overallMinWlen, overallMaxWlen)
    return (minIndex, maxIndex)

def GetPhaseRefractiveIndexRange(mediumProperties):
    """
    Returns the minimum and maximum phase refractive indices as tuple
    over all ice layers for a given mediumProperties object.
    """
    
    minIndex = None
    maxIndex = None
    
    overallMinWlen = mediumProperties.GetMinWavelength()
    overallMaxWlen = mediumProperties.GetMaxWavelength()
    
    for layer in range(mediumProperties.LayersNum):
        minWlen = overallMinWlen
        maxWlen = overallMaxWlen
        
        phaseRefIndex = mediumProperties.GetPhaseRefractiveIndex(layer)
        minWlen = max(minWlen, phaseRefIndex.GetMinWlen())
        maxWlen = min(maxWlen, phaseRefIndex.GetMaxWlen())
        phaseRefIndexFunc = numpy.vectorize(lambda wlen: phaseRefIndex.GetValue(wlen))
        
        testWlens = numpy.linspace(minWlen, maxWlen, 10000, endpoint=True)
        evaluatedFunc = phaseRefIndexFunc(testWlens)
        
        thisLayerMinIndex_arg = numpy.argmin(evaluatedFunc)
        thisLayerMinIndex = evaluatedFunc[thisLayerMinIndex_arg]
        thisLayerMaxIndex_arg = numpy.argmax(evaluatedFunc)
        thisLayerMaxIndex = evaluatedFunc[thisLayerMaxIndex_arg]

        if (minIndex is None) or (thisLayerMinIndex < minIndex):
            minIndex = thisLayerMinIndex
        if (maxIndex is None) or (thisLayerMaxIndex > maxIndex):
            maxIndex = thisLayerMaxIndex

    return (minIndex, maxIndex)

