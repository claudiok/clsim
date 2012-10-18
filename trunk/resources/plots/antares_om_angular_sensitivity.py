#!/usr/bin/env python

#--------------------------------------------
# plot_antares_om_angular_sensitivity
#
# A script to plot the 4 four possible angular
# acceptances for an ANTARES OM.
#
# Further a comparison plot is done to show 
# the precision of the Taylor expansion for the
# 'old' acceptance. In this case a Taylor expansion
# was neccessairy as the acceptance has been
# parametrized as polynomial of arccos, that can not
# be stored in the I3CLSimFunctionPolynomial class 
#--------------------------------------------

import matplotlib
matplotlib.use("PDF")

import matplotlib.pylab as plt

from icecube import icetray, dataclasses
from icecube.clsim import I3CLSimFunctionPolynomial
from icecube.clsim.GetAntaresOMAngularSensitivity import *

params = {'backend': 'pdf',
        'axes.labelsize': 10,
        'text.fontsize': 10,
        'legend.fontsize': 10,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'text.usetex': True}
matplotlib.rcParams.update(params)
matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern']})


# A function definition for the original implementation of the 'old' acceptance
def OriginalWang(x):

        a0=59.115
        a1=0.52258
        a2=0.60944E-02
        a3=-0.16955E-03
        a4=0.60929E-06
        
        if x < -0.36:       
            return 0.
        elif x >= 1.:
            return 1.
        else:
            th = plt.arccos(x) * 57.29578 + 57.75
            wang = a0 + a1*th + a2*th*th + a3*th*th*th + a4*th*th*th*th
            wang /= 84.
            return wang


#Get the implementations of the acceptance
acceptance_old = GetAntaresOMAngularSensitivity('old')
acceptance_NIM = GetAntaresOMAngularSensitivity('NIM')
acceptance_Genova = GetAntaresOMAngularSensitivity('Genova')
acceptance_Spring09 = GetAntaresOMAngularSensitivity('Spring09')


#Evaluate the functions
x =          [float(i)/1000. for i in range(-1000,1001)]
y_old =      [acceptance_old.GetValue(i) for i in x]
y_NIM =      [acceptance_NIM.GetValue(i) for i in x]
y_Genova =   [acceptance_Genova.GetValue(i) for i in x]
y_Spring09 = [acceptance_Spring09.GetValue(i) for i in x]
y_origwang = [OriginalWang(i) for i in x]


#Make the main plot
fig = plt.figure(1, figsize=[8,10])
fig.canvas.set_window_title("Angular sensitivity of an ANTARES OM from different measurements")

plt.subplot(211)
plt.xlabel("$cos(\\theta)$ [rad]")
plt.ylabel("sensitivity")
plt.title("Angular sensitivity of an ANTARES OM\nfrom different measurements")
plt.plot(x, y_old, label='old')
plt.plot(x, y_NIM, label='NIM')
plt.plot(x, y_Genova, label='Genova')
plt.plot(x, y_Spring09, label='Spring09')
plt.legend(loc='upper left')
plt.grid()

plt.subplot(212)
plt.xlabel("$cos(\\theta)$ [rad]")
plt.ylabel("log sensitivity")
plt.plot(x, y_old, label='old')
plt.plot(x, y_NIM, label='NIM')
plt.plot(x, y_Genova, label='Genova')
plt.plot(x, y_Spring09, label='Spring09')
plt.legend(loc='upper left')
plt.grid()
plt.yscale('log')

fig.savefig("antares_om_angular_sensitivity.pdf")

#Make the difference plot for the 'old' acceptance
diff = [y_old[i] - y_origwang[i] for i in range(y_old.__len__())]
fig = plt.figure(2, figsize=[8,10])
fig.canvas.set_window_title("The 'old' ANTARES OM angular sensitivity vs. its Taylor expansion")

plt.subplot(211)
plt.xlabel("$cos(\\theta)$ [rad]")
plt.ylabel("sensitivity")
plt.title("'old' angular sensitivity of an ANTARES OM")
plt.plot(x, y_origwang, label='Original implementation')
plt.plot(x, y_old, label='30 orders Taylor expansion')
plt.legend(loc='upper left')
plt.grid()

plt.subplot(212)
plt.xlabel("$cos(\\theta)$ [rad]")
plt.ylabel("Difference (Taylor - original) in sensitivity")
plt.title("Difference between the original implementation and its Taylor expansion")
plt.plot(x, diff)
plt.grid()

fig.savefig("TaylorExpansionOfOldAntaresOMAngularSensitivity.pdf")

#Show
#plt.show()

