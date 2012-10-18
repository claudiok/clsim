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
# @file GetIceCubeDOMAngularSensitivity.py
# @version $Revision$
# @date $Date$
# @author Claudio Kopper
#

from icecube import icetray, dataclasses
from icecube.clsim import I3CLSimFunctionPolynomial

from I3Tray import I3Units

import numpy, math
from os.path import expandvars


def GetIceCubeDOMAngularSensitivity(holeIce=True):
    """
    The relative collection efficiency of the DOM as a polynomial in the
    cosine of the photon's impact angle with respect to the DOM orientation
    (for IceCube, straight down).
    """
    if holeIce:
        coefficients = [ 0.32813,
                         0.63899,
                         0.20049,
                        -1.2250,
                        -0.14470,
                         4.1695,
                         0.76898,
                        -5.8690,
                        -2.0939,
                         2.3834,
                         1.0435,
                       ]
    else:
        coefficients = [ 0.26266,
                         0.47659,
                         0.15480,
                        -0.14588,
                         0.17316,
                         1.3070,
                         0.44441,
                        -2.3538,
                        -1.3564,
                         1.2098,
                         0.81569,
                       ]
    
    return I3CLSimFunctionPolynomial(coefficients)
