#################################################################
# A python script to return the wavelength-dependant acceptance
# of the ANTARES OM. 
#
# To get the whole acceptance, call the function
# GetAntaresOMAcceptance(omRadius) that calls all other 
# implemented function in this file.
#
# The angular acceptance is not taken into account here.
# That purpose is evaluated by the script
# GetAntaresOMAngularSensitivity.py
#
# The following code has been ported from the Fortran code
# of km3 version v4r0
#
# The table of photo-electron acceptance of the OM is 
# calculated at an injection angle of 0 deg.
#
# The acceptance is calculated the following way
# eff_area = PM collection efficiency
#            * PM quantum efficiency (lambda)
#            * glass + gel transmission probability (lambda)
#
# May 2011
# Florian Folger
# florian.folger@physik.uni-erlangen.de
#################################################################


from icecube import icetray, dataclasses
from icecube.clsim import I3CLSimFunctionFromTable

from I3Tray import I3Units

import numpy, math
from os.path import expandvars

    #################################################################
    # Some getters to store ANTARES specific constants
    # copied from km3 (hit-eff_area_pmt.f and hit-transmit.f)
    #################################################################
def GetAntaresPMTCollectionEfficiency():
    return 0.9

def GetAntaresOMGlassThickness():
    return 1.5*I3Units.cm
    
def GetAntaresOMGelThickness():
    return 1.*I3Units.cm
    
def GetAntaresPMTDiameter():
    return 9.3 * 0.0254*I3Units.m # 9.3 inch PMT
    

    
    
    
#################################################################
# A function to return the quantum efficiency as instance of
# I3CLSimFunctionFromTable
#################################################################
def GetAntaresOMQuantumEfficiency():
    # Data copied from the km3 file hit-ini_optic.f
    # BB5912 from Hamamatsu
    q_eff_0 = 0.01
    q_eff_reverse = [1.988, # at 610 nm
                     2.714, # at 600 nm
                     3.496,
                     4.347,
                     5.166,
                     6.004,
                     6.885,
                     8.105,
                     10.13,
                     13.03,
                     15.29,
                     16.37,
                     17.11,
                     17.86,
                     18.95,
                     20.22,
                     21.26,
                     22.10,
                     22.65,
                     23.07,
                     23.14,
                     23.34,
                     22.95,
                     22.95,
                     22.74,
                     23.48,
                     22.59,
                     20.61,
                     17.68,
                     13.18,
                     7.443,
                     2.526  # at 300 nm
                    ]
                    
    q_eff_reverse = [ (q_eff_0 * i) for i in q_eff_reverse ]
    q_eff_reverse.reverse() # reverse the list (in-place)
    
    return I3CLSimFunctionFromTable(300.*I3Units.nanometer, 10.*I3Units.nanometer, q_eff_reverse)







#################################################################
# A function to return the absoprtion length of
# the glass sphere of an ANTARES OM
#################################################################
def GetAntaresOMGlassAbsorptionLength():

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







#################################################################
# A function to return the absorption length
# the gel of an ANTARES OM
# Note: The file hit-ini_optic.f has three different 
# datasets for this absorption length!
# However in the file hit.f it always is initialized with the
# same (gel_id=1). Thus this one is implemented here.
#################################################################
def GetAntaresOMGelAbsorptionLength():
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
# of the Antares OM
#################################################################
def GetAntaresOMAcceptance(domRadius = 0.2159*I3Units.m): # 17 inch diameter
    
    # Load the constants
    glass_width = GetAntaresOMGlassThickness()
    gel_width = GetAntaresOMGelThickness()
    pmt_collection_efficiency = GetAntaresPMTCollectionEfficiency()
    pmt_diameter = GetAntaresPMTDiameter()
    
    # Geometrical area of the om profile
    pmt_area = math.pi * (pmt_diameter/2.)**2 #is im square meters
    om_area = math.pi*domRadius**2.
    
    # Get the tables from above
    q_eff = GetAntaresOMQuantumEfficiency()
    abs_glass = GetAntaresOMGlassAbsorptionLength()
    abs_gel = GetAntaresOMGelAbsorptionLength()
    
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
            current_om_eff_area = pmt_area * \
                                  pmt_collection_efficiency * \
                                  q_eff.GetValue(wavelength*I3Units.nanometer) * \
                                  math.exp( -( glass_width / this_abs_glass ) ) * \
                                  math.exp( -( gel_width / this_abs_gel ) )

        om_eff_area.append(current_om_eff_area)
    
    
    
    return I3CLSimFunctionFromTable(290.*I3Units.nanometer, 10.*I3Units.nanometer, numpy.array(om_eff_area)/om_area)





