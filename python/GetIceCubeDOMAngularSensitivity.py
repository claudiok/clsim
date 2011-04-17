from icecube import icetray, dataclasses
from icecube.clsim import I3CLSimWlenDependentValuePolynomial

from I3Tray import I3Units

import numpy, math
from os.path import expandvars


def GetIceCubeDOMAngularSensitivity(holeIce=True):
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
    
    return I3CLSimWlenDependentValuePolynomial(coefficients)
