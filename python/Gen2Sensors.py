
from icecube.icetray import I3Units
from icecube.clsim import I3CLSimFunctionFromTable, I3CLSimFunctionPolynomial

def GetWOMAcceptance(active_fraction=1.):
    """
    :param active_fraction: the fraction of the head-on geometric area that
       is photosensitive (i.e. the ratio of the area of the wavelength-shifting)
       tube to the area of the housing.
    """
    # probability that Cherenkov photon directly incident on the wavelength-
    # shifting paint is absorbed, and that the shifted photon is captured in
    # the plexiglass tube.
    # Pers. comm., D. Hebecker, April 2016
    capture_efficiency = [0.0,
                          0.34587,
                          0.45655,
                          0.48452,
                          0.46706,
                          0.47998,
                          0.48761,
                          0.48948,
                          0.49017,
                          0.4905,
                          0.49127,
                          0.49325,
                          0.4966,
                          0.49651,
                          0.4857,
                          0.40011,
                          0.15273,
                          0.00779,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0]
    # KM3NeT PMT QE, weighted with the emission spectrum of the wavelength
    # shifter
    recapture_efficiency = 0.2403
    return I3CLSimFunctionFromTable(245*I3Units.nanometer, 10*I3Units.nanometer,
        [a*recapture_efficiency*active_fraction for a in capture_efficiency])

def GetWOMAngularSensitivity():
    # transmission efficiency from ice into a quartz glass tube, averaged over
    # the cross-section, and multiplied by sin(eta) to account for projected
    # area
    coefficients = [0.70161228651625462,
                    0.0,
                    -0.78196095712541591,
                    0.0,
                    1.9327345553744812,
                    0.0,
                    -14.801481314906798,
                    0.0,
                    37.180692649664785,
                    0.0,
                    -34.627444106282297]
    return I3CLSimFunctionPolynomial(coefficients, -1./1.33, 1/1.33, 0, 0)
