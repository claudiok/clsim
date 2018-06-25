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
# $Id: GetIceCubeDOMAngularSensitivity.py 108199 2013-07-12 21:33:08Z nwhitehorn $
# 
# @file GetIceCubeDOMAngularSensitivity.py
# @version $Revision: 108199 $
# @date $Date: 2013-07-12 17:33:08 -0400 (Fri, 12 Jul 2013) $
# @author Claudio Kopper
#

from icecube import icetray, dataclasses
from icecube.clsim import I3CLSimFunctionPolynomial

from I3Tray import I3Units

import numpy, math
from os.path import expandvars


def GetIceCubeDOMAngularSensitivity(holeIce=expandvars("$I3_BUILD/ice-models/resources/models/angsens/as.nominal")):
    """
    The relative collection efficiency of the DOM as a polynomial in the
    cosine of the photon's impact angle with respect to the DOM orientation
    (for IceCube, straight down).
    """
    
    coefficients = numpy.loadtxt(holeIce)[1:].tolist()
    
    return I3CLSimFunctionPolynomial(coefficients)
