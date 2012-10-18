#!/usr/bin/env python

#--------------------------------------------
# plot_icecube_dom_acceptance.py
#
# A script to plot the wavelength dependent
# acceptance of the IceCube DOM as implemented
# in python/GetIceCubeDOMAcceptance.py
#--------------------------------------------

import matplotlib
matplotlib.use("PDF")

import matplotlib.pylab as plt

from icecube import icetray, dataclasses
from icecube.clsim import I3CLSimFunctionFromTable
from icecube.clsim.GetIceCubeDOMAcceptance import *

params = {'backend': 'pdf',
        'axes.labelsize': 10,
        'text.fontsize': 10,
        'legend.fontsize': 10,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'text.usetex': True}
matplotlib.rcParams.update(params)
matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern']})

#Get the acceptance
dom_acceptance = GetIceCubeDOMAcceptance()

#Evaluate the function
x =          [float(i)/10. for i in range(2500,7001)]
dom_acc =    [dom_acceptance.GetValue(i*I3Units.nanometer) for i in x]

#Make the plot
fig = plt.figure(1, figsize=[10,7])
fig.canvas.set_window_title("Acceptance of an IceCube DOM")

plt.subplot(111)
plt.xlabel("Wavelength $\\lambda$ [nm]")
plt.ylabel("Acceptance")
plt.title("Acceptance of an IceCube DOM as function of wavelength")
plt.plot(x, dom_acc)
plt.grid()



fig.savefig("icecube_dom_acceptance.pdf")




