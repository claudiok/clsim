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
# $Id: GetAntaresOMAngularSensitivity.py 108199 2013-07-12 21:33:08Z nwhitehorn $
# 
# @file GetAntaresOMAngularSensitivity.py
# @version $Revision: 108199 $
# @date $Date: 2013-07-12 17:33:08 -0400 (Fri, 12 Jul 2013) $
# @author Claudio Kopper
#

from __future__ import print_function

import numpy, math, sys
from os.path import expandvars

from icecube import icetray, dataclasses
from icecube.clsim import I3CLSimFunctionPolynomial

from I3Tray import I3Units


def GetAntaresOMAngularSensitivity(name='NotSet'):
    """
    A python function to return the Antares OM angular acceptance
    stored in an instance of I3CLSimWlenDependetValuePolynomial
    """
    
    if name == 'NotSet':
        print("GetAntaresOMAngularSensitivity ERROR: Please specify an angulare sensitivity type! Possible are 'Spring09', 'Genova', 'NIM' and 'old'.")
        sys.exit()
        
    elif name == 'Spring09':
        # Parameterisation of angular acceptance, "spring 2009"
        # AUTHOR: Heide Costantini, December, 2008,
        # also named dic08
        #
        # Fit from Genova - GEANT4
        # GE6 param, cutoff at costheta=-0.65
        
        coefficients = [ 0.3265,
                         0.6144,
                        -0.0343,
                        -0.0641,
                         0.2988,
                        -0.1422
                       ]
    
        return I3CLSimFunctionPolynomial(coefficients, -0.65, 1., 0., 1.)

    elif name == 'Genova':
        # Genova parameterisation of angular acceptance
        # AUTHOR: M. Anghinolfi, H. Costantini, October, 2007
        # Values from H. Constantini's talk in Sinaia
        #
        # cutoff at costheta=-0.80
        
        coefficients = [ 0.349,
                         0.547,
                         0.063,
                        -0.036,
                         0.077
                       ]
                       
        return I3CLSimFunctionPolynomial(coefficients, -0.80, 1., 0., 1.)
        
    elif name == 'NIM':
        # NIM parameterisation of angular acceptance
        # AUTHOR: M. Spurio, Apr. 26th, 2007
        #
        # cutoff at costheta=-0.65
        
        coefficients = [ 0.2549,
                         0.6093,
                         0.2556,
                        -0.1231
                       ]
                       
        return I3CLSimFunctionPolynomial(coefficients, -0.65, 1., 0., 1.)

    elif name == 'old':
        # calculate weight for photon which hits PM under
        # angle th (degree) w.r.t. axis with cos(th) = ct 
        #
        # The coefficients in this case are a Taylor expansion derived from
        # the following original code. From the taylor expansion 30 terms are 
        # taken into account to provide an accurancy that is sufficient enough.
        #
        # const double a0=59.115;
        # const double a1=0.52258;
        # const double a2=0.60944E-02;
        # const double a3=-0.16955E-03;
        # const double a4=0.60929E-06;
        #
        # th = std::acos(ct)*57.29578 + 57.75;
        # double wang = a0+ a1*th + a2*th*th + a3*th*th*th + a4*th*th*th*th;
        # wang/=84.;
        #
        # The cutoff is at -0.36
        
        coefficients = [ 0.153099,
                         0.627246,
                         0.41998,
                        -0.322113,
                         0.218163,
                        -0.166283,
                         0.126776,
                        -0.10355,
                         0.0844767,
                        -0.0720585,
                         0.0612634,
                        -0.0537683,
                         0.0469892,
                        -0.042072,
                         0.0374956,
                        -0.0340695,
                         0.0308118,
                        -0.0283139,
                         0.0258992,
                        -0.0240126,
                         0.0221646,
                        -0.0206989,
                         0.0192477,
                        -0.0180824,
                         0.0169184,
                        -0.0159738,
                         0.0150234,
                        -0.0142452,
                         0.0134573,
                        -0.0128072,
                         0.0121454
                        ]
                        

        return I3CLSimFunctionPolynomial(coefficients, -0.36, 1., 0., 1.)
    
    else:
        print("GetAntaresOMAngularSensitivity ERROR: '"+name+"' does not name an angulare sensitivity type! Possible are 'Spring09', 'Genova', 'NIM' and 'old'.")
        sys.exit()



