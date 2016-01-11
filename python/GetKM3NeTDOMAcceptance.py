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
# $Id: GetKM3NeTDOMAcceptance.py 108199 2013-07-12 21:33:08Z nwhitehorn $
# 
# @file GetKM3NeTDOMAcceptance.py
# @version $Revision: 108199 $
# @date $Date: 2013-07-12 17:33:08 -0400 (Fri, 12 Jul 2013) $
# @author Claudio Kopper
#

"""
A python script to return the wavelength-dependant acceptance
of the KM3NeT DOM. 

To get the whole acceptance, call the function
GetKM3NeTDOMAcceptance(omRadius) that calls all other 
implemented function in this file.

The angular acceptance is not taken into account here.

The table of photo-electron acceptance of the OM is 
calculated at an injection angle of 0 deg.
"""


from icecube import icetray, dataclasses
from icecube.clsim import I3CLSimFunctionFromTable

from I3Tray import I3Units

import numpy, math
from os.path import expandvars

#################################################################
# Some getters to store KM3NeT specific constants
#################################################################
def GetKM3NeTPMTCollectionEfficiency():
    return 0.9

def GetKM3NeTOMGlassThickness():
    return 1.5*I3Units.cm
    
def GetKM3NeTOMGelThickness():
    return 1.*I3Units.cm
    

    
    
    
def GetKM3NeTOMQuantumEfficiency(peakQE=None, wpdQE=False):
    """
    A function to return the quantum efficiency as instance of
    I3CLSimFunctionFromTable
    """
    if peakQE is None:
        if wpdQE: peakQE = 0.304
        else: peakQE = 0.32

    if wpdQE:
        # this is taken straight from the WPD document
        QEscale = peakQE/0.304
        
        QEvals = [ 0.0,  0.0,  0.5,  3.1,  9.8, 17.5, 23.2, 26.5, 28.1, 28.1,
                  29.1, 30.1, 30.4, 30.1, 29.9, 29.3, 28.6, 27.5, 26.5, 25.0,
                  23.2, 21.1, 19.6, 18.5, 17.2, 15.4, 12.1,  9.3,  7.2,  6.2,
                   4.6,  3.6,  2.8,  2.1,  1.3,  0.8,  0.5,  0.3,  0.0,  0.0]
        QEvals = numpy.array(QEvals)*0.01 * QEscale
        
        return I3CLSimFunctionFromTable(260.*I3Units.nanometer, 10.*I3Units.nanometer, QEvals)
        
    else:
        # this is the version we used before the WPD document
        QEscale = peakQE
        
        QEvals = [0.00, 0.87, 1.00, 0.94, 0.78, 0.49, 0.24, 0.09, 0.02, 0.00]
        QEvals = numpy.array(QEvals)*QEscale
        
        return I3CLSimFunctionFromTable(250.*I3Units.nanometer, 50.*I3Units.nanometer, QEvals)




def GetKM3NeTOMGlassAbsorptionLength():
    """
    A function to return the absoprtion length of
    the glass sphere of an KM3NeT OM
    """

    # Data copied from the km3 file hit-ini_optic.f
    # as measured by Pavel
    # the last 3 bins contain measurements from Saclay (6/3/98)
    al_glass_reverse = [148.37, # at 610 nm
                        142.87, # at 600 nm
                        135.64,
                        134.58,
                        138.27,
                        142.40,
                        147.16,
                        151.80,
                        150.88,
                        145.68,
                        139.70,
                        126.55,
                        118.86,
                        113.90,
                        116.08,
                        109.23,
                        81.63,
                        65.66,
                        77.30,
                        73.02,
                        81.25,
                        128.04,
                        61.84,
                        19.23,
                        27.21, 
                        18.09,  
                        8.41,  
                        3.92,  
                        1.82,  
                        0.84,  
                        0.39, 
                        0.17 # at 300 nm
                       ]
                       
    # Apply units
    al_glass_reverse = [ (i * I3Units.cm) for i in al_glass_reverse]
    al_glass_reverse.reverse() # reverse the list (in-place)

    return I3CLSimFunctionFromTable(300.*I3Units.nanometer, 10.*I3Units.nanometer, al_glass_reverse)            





def GetKM3NeTOMGelAbsorptionLength():
    """
    A function to return the absorption length
    the gel of an KM3NeT OM
    Note: The file hit-ini_optic.f has three different 
    datasets for this absorption length!
    However in the file hit.f it always is initialized with the
    same (gel_id=1). Thus this one is implemented here.
    """
    
    # Data copied from the km3 file hit-ini_optic.f
    # GEL WACKER (default)
    al_gel_default_reverse = [100.81, # at 610 nm
                              99.94, # at 600 nm
                              99.89, 
                              96.90, 
                              96.42, 
                              94.36, 
                              89.09, 
                              90.10,
                              86.95, 
                              85.88, 
                              84.49, 
                              81.08, 
                              78.18, 
                              76.48, 
                              74.55, 
                              72.31,
                              68.05, 
                              66.91, 
                              64.48, 
                              62.53, 
                              59.38, 
                              56.64, 
                              53.29, 
                              48.96,
                              45.71, 
                              41.88, 
                              37.14, 
                              30.49, 
                              23.08, 
                              15.60,  
                              8.00,  
                              0.00 # at 300 nm
                             ]

    # Apply units
    al_gel_default_reverse = [ (i * I3Units.cm) for i in al_gel_default_reverse]
    al_gel_default_reverse.reverse() # reverse the list (in-place)
    
    return I3CLSimFunctionFromTable(300.*I3Units.nanometer, 10.*I3Units.nanometer, al_gel_default_reverse)







#################################################################
# The main function to return the effective area
# of the KM3NeT DOM
#################################################################
def GetKM3NeTDOMAcceptance(domRadius = (17./2.) * 0.0254*I3Units.m, peakQE=None, wpdQE=False, withWinstonCone=False): # 17 inch diameter
    """
    The main function to return the effective area
    of the KM3NeT DOM
    """
    
    # the multiPMT simulation expects photons on a sphere of 17" diameter.
    # The result of this function is supposed to be used for emission
    # spectrum biassing only, so only allow 17" spheres just to be safe.
    if domRadius != (17./2.) * 0.0254*I3Units.m:
        raise RuntimeError("GetKM3NeTDOMAcceptance is currently only valid for a diameter of 17inch.")
    
    # Load the constants
    glass_width = GetKM3NeTOMGlassThickness()
    gel_width = GetKM3NeTOMGelThickness()
    pmt_collection_efficiency = GetKM3NeTPMTCollectionEfficiency()
    
    # Get the tables from above
    q_eff = GetKM3NeTOMQuantumEfficiency(peakQE=peakQE, wpdQE=wpdQE)
    abs_glass = GetKM3NeTOMGlassAbsorptionLength()
    abs_gel = GetKM3NeTOMGelAbsorptionLength()
    
    if withWinstonCone:
        coneScaler = 2.0 # correction factor at cos(theta)=1
    else:
        coneScaler = 1.
        
    # Each of the tables above has 32 bins in the wavelength range of 300nm - 610nm,
    # where the exact value is at the beginning of each bin, means
    # value of bin  0 belongs to 300nm
    # value of bin  1 belongs to 310nm
    # value of bin 31 belongs to 610nm
    # Now combine them
    om_eff_area = [0.] # use a single entry at 290nm to have the same range as other functions(wlen)
    for wavelength in range(300, 611, 10):
        this_abs_glass = abs_glass.GetValue(wavelength*I3Units.nanometer)
        this_abs_gel = abs_gel.GetValue(wavelength*I3Units.nanometer)
        
        if (this_abs_glass <= 0.) or (this_abs_gel <= 0.):
            current_om_eff_area = 0.
        else:
            current_om_eff_area = pmt_collection_efficiency * \
                                  q_eff.GetValue(wavelength*I3Units.nanometer) * \
                                  coneScaler
                                  # do not include these two (they assume a fixed photon path length through glass which
                                  # currently is not the smallest possible length..
                                  #math.exp( -( glass_width / this_abs_glass ) ) * \
                                  #math.exp( -( gel_width / this_abs_gel ) )

        om_eff_area.append(current_om_eff_area)
    
    
    
    return I3CLSimFunctionFromTable(290.*I3Units.nanometer, 10.*I3Units.nanometer, numpy.array(om_eff_area))





