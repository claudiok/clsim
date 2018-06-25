#!/usr/bin/env python

#--------------------------------------------
# icecube_dom_angular_sensitivity
#
# A script to plot the icecube
# dom angular sensitivity (holeIce/not holeIce)
#--------------------------------------------

import matplotlib
matplotlib.use("PDF")

import matplotlib.pylab as plt

from os.path import expandvars
from icecube import icetray, dataclasses
from icecube.clsim import I3CLSimFunctionPolynomial
from icecube.clsim.GetIceCubeDOMAngularSensitivity import *

params = {'backend': 'pdf',
        'axes.labelsize': 10,
        'text.fontsize': 10,
        'legend.fontsize': 10,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'text.usetex': True}
matplotlib.rcParams.update(params)
matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern']})


#Get the implementations of the acceptance
acceptance = GetIceCubeDOMAngularSensitivity(holeIce=expandvars("$I3_BUILD/ice-models/resources/models/angsens/as.nominal"))
acceptance_holeIce = GetIceCubeDOMAngularSensitivity(holeIce=expandvars("$I3_BUILD/ice-models/resources/models/angsens/as.h2-50cm"))



#Evaluate the functions
x =          [float(i)/1000. for i in range(-1000,1001)]
y =          [acceptance.GetValue(i) for i in x]
y_holeIce =  [acceptance_holeIce.GetValue(i) for i in x]



#Make the main plot
fig = plt.figure(1, figsize=[8,10])
fig.canvas.set_window_title("Angular sensitivity of an IceCube DOM")

plt.subplot(211)
plt.xlabel("$cos(\\theta)$ [rad]")
plt.ylabel("sensitivity")
plt.title("Angular sensitivity of an IceCube DOM")
plt.plot(x, y, label='without holeIce')
plt.plot(x, y_holeIce, label='with holeIce')
plt.legend(loc='upper left')
plt.grid()

plt.subplot(212)
plt.xlabel("$cos(\\theta)$ [rad]")
plt.ylabel("log sensitivity")
plt.title("Angular sensitivity of an IceCube DOM")
plt.plot(x, y, label='without holeIce')
plt.plot(x, y_holeIce, label='with holeIce')
plt.legend(loc='upper left')
plt.grid()
plt.yscale('log')

fig.savefig("icecube_dom_angular_sensitivity.pdf")



#Show
#plt.show()

