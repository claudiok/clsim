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
# $Id: I3CLSimRandomValueIceCubeFlasherTimeProfile.py 108199 2013-07-12 21:33:08Z nwhitehorn $
# 
# @file I3CLSimRandomValueIceCubeFlasherTimeProfile.py
# @version $Revision: 108199 $
# @date $Date: 2013-07-12 17:33:08 -0400 (Fri, 12 Jul 2013) $
# @author Claudio Kopper
#

from icecube import icetray, dataclasses

from icecube.clsim import I3CLSimRandomValue
from icecube.clsim import I3CLSimRandomValueInterpolatedDistribution
from icecube.clsim.util.interpolate import interp1d

from I3Tray import I3Units

import numpy, math

class I3CLSimRandomValueIceCubeFlasherTimeProfile(I3CLSimRandomValue):
    """
    Samples from a IceCube flasher time distribution.
    The parameter is the configured width of a flasher specified in
    units of time. This is the time delay from the start of the flasher
    light pulse.
    
    To convert from an I3FlasherInfo width (1-127) [in units of 0.5ns], use:
    
    p = width*0.5*I3Units.ns
    
    Distributions should be similar to what is found on
    http://wiki.icecube.wisc.edu/index.php/LED_output_time_profile
    """
    
    ## static data
    _pulse_FB_WIDTH15 = numpy.array(
          [[  0.00000000e+00,   1.00000000e+00,   2.00000000e+00,
              3.00000000e+00,   4.00000000e+00,   5.00000000e+00,
              6.00000000e+00,   7.00000000e+00,   8.00000000e+00,
              9.00000000e+00,   1.00000000e+01,   1.10000000e+01,
              1.20000000e+01,   1.30000000e+01,   1.40000000e+01,
              1.50000000e+01,   1.60000000e+01,   1.70000000e+01,
              1.80000000e+01,   1.90000000e+01,   2.00000000e+01,
              2.10000000e+01,   2.20000000e+01,   2.30000000e+01,
              2.40000000e+01,   2.50000000e+01,   2.60000000e+01,
              2.70000000e+01,   2.80000000e+01,   2.90000000e+01,
              3.00000000e+01,   3.10000000e+01,   3.20000000e+01,
              3.30000000e+01,   3.40000000e+01,   3.50000000e+01,
              3.60000000e+01,   3.70000000e+01,   3.80000000e+01,
              3.90000000e+01,   4.00000000e+01,   4.10000000e+01,
              4.20000000e+01,   4.30000000e+01,   4.40000000e+01,
              4.50000000e+01,   4.60000000e+01,   4.70000000e+01,
              4.80000000e+01,   4.90000000e+01,   5.00000000e+01],
           [  1.18000000e-03,   2.76900000e-02,   1.25170000e-01,
              2.14840000e-01,   3.20890000e-01,   4.32390000e-01,
              4.64370000e-01,   5.00230000e-01,   4.31610000e-01,
              3.16210000e-01,   2.29650000e-01,   1.37640000e-01,
              8.77400000e-02,   7.21400000e-02,   5.96600000e-02,
              4.79700000e-02,   4.09500000e-02,   2.92500000e-02,
              3.08100000e-02,   2.84700000e-02,   2.61300000e-02,
              1.83400000e-02,   1.83400000e-02,   1.99000000e-02,
              1.28800000e-02,   1.28800000e-02,   1.28800000e-02,
              1.60000000e-02,   1.44400000e-02,   1.67800000e-02,
              7.42000000e-03,   6.64000000e-03,   9.76000000e-03,
              1.13200000e-02,   7.42000000e-03,   9.76000000e-03,
              4.30000000e-03,   5.86000000e-03,   7.42000000e-03,
              4.30000000e-03,   8.20000000e-03,   5.86000000e-03,
              3.52000000e-03,   1.96000000e-03,   2.74000000e-03,
              4.30000000e-03,   5.08000000e-03,   2.74000000e-03,
              3.52000000e-03,   4.30000000e-03,   2.74000000e-03]])
    #adjust zero offset and re-scale to one
    _pulse_FB_WIDTH15[1]=(_pulse_FB_WIDTH15[1]-0.00118)/0.49905
    _pulse_narrow = interp1d(_pulse_FB_WIDTH15[0], _pulse_FB_WIDTH15[1], kind='linear', bounds_error=False, fill_value=0.)
    
    @staticmethod
    def _pulse_rising_edge(x, width):
        templateWidth=7.
        scaledX = templateWidth*x/width
    
        scaledX = numpy.where(scaledX>templateWidth,templateWidth,scaledX)
        scaledX = numpy.where(scaledX<0.,0.,scaledX)
    
        return I3CLSimRandomValueIceCubeFlasherTimeProfile._pulse_narrow(scaledX)
    
    @staticmethod
    def _pulse_falling_edge(x):
        templateStart=7.
    
        scaledX = x+templateStart
        scaledX = numpy.where(scaledX<templateStart,templateStart,scaledX)
    
        return I3CLSimRandomValueIceCubeFlasherTimeProfile._pulse_narrow(scaledX)
    
    @staticmethod
    def _pulse_plateau_width(FB_WIDTH):
        return (FB_WIDTH-15.)*59.5/(124.-15.)
    
    @staticmethod
    def _pulse_rising_edge_width(FB_WIDTH):
        return numpy.log((FB_WIDTH-12.))*1.91+5.
    
    @staticmethod
    def _the_pulse(x, FB_WIDTH):
        if isinstance(x, float):
            useX = numpy.array([x])
        else:
            useX = x
        
        if FB_WIDTH <= 15:
            return I3CLSimRandomValueIceCubeFlasherTimeProfile._pulse_narrow(useX*(15./FB_WIDTH))
        else:
            plateau_width = I3CLSimRandomValueIceCubeFlasherTimeProfile._pulse_plateau_width(FB_WIDTH)
            rising_edge_width = I3CLSimRandomValueIceCubeFlasherTimeProfile._pulse_rising_edge_width(FB_WIDTH)
            return numpy.where(useX<=rising_edge_width,I3CLSimRandomValueIceCubeFlasherTimeProfile._pulse_rising_edge(useX,rising_edge_width), numpy.where(useX<=rising_edge_width+plateau_width,numpy.ones(len(useX)), I3CLSimRandomValueIceCubeFlasherTimeProfile._pulse_falling_edge(useX-rising_edge_width-plateau_width) ) )
    
    
    def __init__(self):
        I3CLSimRandomValue.__init__(self)
        self.distributionCache=dict()
        
    def SampleFromDistribution(self, random, parameters):
        if len(parameters) != self.NumberOfParameters():
            raise RuntimeError("Expected %u parameters but got %u." % (self.NumberOfParameters(), len(parameters)))
        width = parameters[0]/I3Units.ns
        
        # This is kind of evil, but it'll work as long as we don't want
        # to use this from python: create a new random distribution for each
        # requested width and sample from it. Cache a few distributions.
        if width not in self.distributionCache:
            # flush the cache in case it gets too large
            if len(self.distributionCache) > 200:
                del self.distributionCache
                self.distributionCache=dict()
            
            maxDuration=120.
            xVals = numpy.linspace(0.,maxDuration,int(maxDuration*2),endpoint=False)
            dist = I3CLSimRandomValueInterpolatedDistribution(0., 0.5, I3CLSimRandomValueIceCubeFlasherTimeProfile._the_pulse(xVals, width*2.))
            self.distributionCache[width] = dist
            
        return self.distributionCache[width].SampleFromDistribution(random, [])
    
    def NumberOfParameters(self):
        return 1
    
    def OpenCLFunctionWillOnlyUseASingleRandomNumber(self):
        return False
    
    def GetOpenCLFunction(self, functionName, functionArgs, functionArgsToCall, uniformRandomCall_co, uniformRandomCall_oc):
        raise RuntimeError("Generating OpenCL code is currently not implemented for I3CLSimRandomValueIceCubeFlasherTimeProfile")
    
    def CompareTo(self, other):
        if not isinstance(other, I3CLSimRandomValueIceCubeFlasherTimeProfile): return False
        return True
