#!/usr/bin/env python

#--------------------------------------------
# plot_antares_om_acceptance.py
#
# A script to plot the wavelength dependent
# acceptance of the ANTARES OM as implemented
# in python/GetAntaresOMAcceptance.py
#
# In this script the acceptance is a
# combination of light absoprtion in the 
# glass and gel of the OM and the
# quantum efficiency of the PM.
#--------------------------------------------

import matplotlib
matplotlib.use("PDF")

import matplotlib.pylab as plt

import numpy

from icecube import icetray, dataclasses
from icecube.clsim import I3CLSimFunctionFromTable
from icecube.clsim.GetAntaresOMAcceptance import *

params = {'backend': 'pdf',
        'axes.labelsize': 10,
        'text.fontsize': 10,
        'legend.fontsize': 10,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'text.usetex': True}
matplotlib.rcParams.update(params)
matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern']})

#Get things from GetAntaresOMAcceptance
quantum_efficiency = GetAntaresOMQuantumEfficiency()
glass_absorption_length = GetAntaresOMGlassAbsorptionLength()
gel_absorption_length = GetAntaresOMGelAbsorptionLength()
om_acceptance = GetAntaresOMAcceptance()

#Load the constants
pm_collection_efficiency = GetAntaresPMTCollectionEfficiency()
glass_width = GetAntaresOMGlassThickness()
gel_width = GetAntaresOMGelThickness()

pmt_diameter = 9.3 * 0.0254*I3Units.m
pmt_area = math.pi * (pmt_diameter/2.)**2
domRadius = 0.2159*I3Units.m
om_area = math.pi*domRadius**2.
pm_area_factor = pmt_area/om_area

#Evaluate the functions
x =          [float(i)/10. for i in range(2800,6201)]
qe =         [quantum_efficiency.GetValue(i*I3Units.nanometer) for i in x]
glass_abs =  [glass_absorption_length.GetValue(i*I3Units.nanometer) for i in x]
gel_abs =    [gel_absorption_length.GetValue(i*I3Units.nanometer) for i in x]
om_acc =     [om_acceptance.GetValue(i*I3Units.nanometer) for i in x]
pm_coll =    [pm_collection_efficiency for i in x]
area_factor =    [pm_area_factor for i in x]

#Calculate some more
glass_trans = [plt.exp(-(glass_width / i)) if i != 0 else 0 for i in glass_abs]
gel_trans =   [plt.exp(-(gel_width / i)) if i != 0 else 0 for i in gel_abs] 

#Make the main plot
fig = plt.figure(1, figsize=[8,10])
fig.canvas.set_window_title("Acceptance properties of an ANTARES OM")

plt.subplot(311)
plt.xlabel("Wavelength $\\lambda$ [nm]")
plt.ylabel("Absorption length [m]")
#plt.title("Absorption lengths in an ANTARES OM as function of wavelength")
plt.plot(x, glass_abs, label=str(glass_width*100.)+' cm Glass')
plt.plot(x, gel_abs, label=str(gel_width*100.)+' cm Gel')
plt.legend(loc='upper left')
plt.grid()

plt.subplot(312)
plt.xlabel("wavelength $\\lambda$ [nm]")
plt.ylabel("Efficiency")
#plt.title("Properties of the ANTARES OM acceptance")
plt.plot(x, glass_trans, label='OM glass transmissibility ('+str(glass_width*100.)+' cm)')
plt.plot(x, gel_trans, label='OM gel transmissibility ('+str(gel_width*100.)+' cm)')
plt.plot(x, qe, label='PM Quantum efficiency')
plt.plot(x, pm_coll, label='PMT collection efficiency')
plt.plot(x, area_factor, label=r'PMT to OM area factor (17inch $\rightarrow$ 9.3inch)')
plt.plot(x, om_acc, label='overall OM acceptance', linewidth=3, color='y')
plt.legend(loc='upper right')
plt.grid()


plt.subplot(313)
plt.xlabel("wavelength $\\lambda$ [nm]")
plt.ylabel("Efficiency")
#plt.title("Properties of the ANTARES OM acceptance")
plt.plot(x, om_acc, label='overall OM acceptance', linewidth=3, color='y')
plt.plot(x, numpy.array(om_acc)*(17.**2)/(13.**2), label='overall OM acceptance (scaled to IceCube DOM radius)', linewidth=1)
plt.ylim(0.,0.14)
plt.legend(loc='upper right')
plt.grid()

fig.savefig("antares_om_acceptance.pdf")


