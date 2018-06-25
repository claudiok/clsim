
from icecube.icetray import I3Units, logging
from icecube.dataclasses import I3Constants
from icecube.clsim import I3CLSimFunctionFromTable, I3CLSimFunctionPolynomial

from os.path import expandvars

def GetDEggAcceptance(active_fraction=1.):
    """
    :param active_fraction: the fraction of the head-on geometric area that
       is photosensitive (i.e. the ratio of the photocathode area to the
       geometric area)
    """
    # Combined efficiency for D-Egg glass (10 mm), high-UV transparency gel
    # (5 mm), and Hamamatsu R5912-100 at the center of the photocathode
    # Pers. comm., Lu Lu, April 2016

    logging.log_warn("DeprecationWarning: The numbers here are outdated! Also this functionality is now part of the 'mDOM-WOM-simulation' project. Please check it out!")

    center_efficiency = [0.0,
                         0.0,
                         0.0,
                         0.0005,
                         0.0093,
                         0.058,
                         0.1473,
                         0.2358,
                         0.2904,
                         0.3139,
                         0.3237,
                         0.3336,
                         0.339,
                         0.3373,
                         0.3292,
                         0.3195,
                         0.3087,
                         0.3017,
                         0.2873,
                         0.2717,
                         0.2532,
                         0.2305,
                         0.2119,
                         0.1962,
                         0.1832,
                         0.1708,
                         0.1523,
                         0.1227,
                         0.0928,
                         0.0728,
                         0.0597,
                         0.0494,
                         0.0404,
                         0.0318,
                         0.0241,
                         0.0174,
                         0.0118,
                         0.0076,
                         0.0047,
                         0.0027,
                         0.0,
                         0.0,
                         0.0,
                         0.0]
    # Hamamatsu quotes a minimum photocathode diameter of 190 mm. Pending a real
    # 2D average of the capture efficiency, approximate with 90% of the center
    # efficiency times the minimum photocathode area.
    active_fraction *= 0.9*(190./300.)**2
    return I3CLSimFunctionFromTable(250*I3Units.nanometer, 10*I3Units.nanometer,
        [a*active_fraction for a in center_efficiency])

def GetDEggAngularSensitivity(pmt='both'):
    
    logging.log_warn("DeprecationWarning: The numbers here are outdated! Also this functionality is now part of the 'mDOM-WOM-simulation' project. Please check it out!")

    import numpy
    from icecube.clsim import GetIceCubeDOMAngularSensitivity
    
    angularAcceptance = GetIceCubeDOMAngularSensitivity(holeIce=expandvars("$I3_BUILD/ice-models/resources/models/angsens/as.nominal"))
    
    # mirror the function in cos(eta) by inverting the odd components
    coeffs = numpy.array(angularAcceptance.GetCoefficients())
    coeffs[numpy.arange(coeffs.size) % 2 == 1] *= -1
    
    if pmt.lower() == 'down':
        return angularAcceptance
    elif pmt.lower() == 'up':
        return I3CLSimFunctionPolynomial(coeffs)
    elif pmt.lower() == 'both':
        return I3CLSimFunctionPolynomial(numpy.array(angularAcceptance.GetCoefficients()) + coeffs)
    else:
        raise ValueError("Unknown PMT orientation '%s'" % pmt)

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

