#!/usr/bin/env python

from __future__ import print_function

import math
import numpy

import matplotlib
matplotlib.use("PDF")

fig_size = [8.3,11.7] # din A4
params = {'backend': 'pdf',
        'axes.labelsize': 10,
        'text.fontsize': 10,
        'legend.fontsize': 10,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'text.usetex': True,
        'figure.figsize': fig_size}
matplotlib.rcParams.update(params)
matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern']})

import pylab
import scipy
import scipy.interpolate
import scipy.integrate

nm=1.
m=1.
wlens=numpy.array(
      [290.*nm, 310.*nm, 330.*nm, 350.*nm, 370.*nm,
       390.*nm, 410.*nm, 430.*nm, 450.*nm, 470.*nm,
       490.*nm, 510.*nm, 530.*nm, 550.*nm, 570.*nm,
       590.*nm, 610.*nm]
    )
scatLength=[16.67194612*m, 20.24988356*m, 23.828246*m, 27.60753133*m, 31.54474622*m,
            35.6150723*m, 39.79782704*m, 44.07227854*m, 48.42615012*m, 52.84574328*m, 
            57.31644409*m, 61.83527084*m, 66.38783775*m, 70.97232079*m, 75.58007709*m, 
            80.20532563*m, 84.84642797*m]
refIndex=[1.374123775, 1.368496907, 1.364102384, 1.360596772, 1.357746292,
         1.355388160, 1.353406686, 1.351718123, 1.350260806, 1.348988618,
         1.347866574, 1.346867782, 1.345971300, 1.345160644, 1.344422686,
         1.343746868, 1.343124618]

# very old model, do not use!
absLength=[4.750413286*m, 7.004812306*m, 9.259259259*m, 14.92537313*m, 20.00000000*m, 
          26.31578947*m, 34.48275862*m, 43.47826087*m, 50.00000000*m, 62.50000000*m, 
          58.82352941*m, 50.00000000*m, 29.41176471*m, 17.85714286*m, 16.12903226*m, 
          8.849557522*m, 4.504504505*m]

# taken from Geasim (used by km3 by default)
# original comment says:
# "mix from Antares and Smith-Baker (smallest value for each bin)"
absCoeffGeasim=[0.2890,0.2440,0.1570,0.1080,0.0799,0.0708
               ,0.0638,0.0558,0.0507,0.0477,0.0357,0.0257,0.0196
               ,0.0182,0.0182,0.0191,0.0200,0.0218,0.0237,0.0255
               ,0.0291,0.0325,0.0363,0.0415,0.0473,0.0528,0.0629
               ,0.0710,0.0792,0.0946,0.1090,0.1390,0.215]
wlenGeasim=    [610.,  600.,  590.,  580.,  570.,  560.,
                550.,  540.,  530.,  520.,  510.,  500.,  490.,
                480.,  470.,  460.,  450.,  440.,  430.,  420.,
                410.,  400.,  390.,  380.,  370.,  360.,  350.,
                340.,  330.,  320.,  310.,  300.,  290.]
interpolatedAbsorptionCoeffsGeasim = scipy.interpolate.interp1d(wlenGeasim[::-1], absCoeffGeasim[::-1], bounds_error=False)


AbsorptionCoeffsSmithBaker = numpy.loadtxt("ExistingData/AbsorptionCoefficients_SmithBaker.txt", unpack=True)
interpolatedAbsorptionCoeffsSmithBaker = scipy.interpolate.interp1d(AbsorptionCoeffsSmithBaker[0], AbsorptionCoeffsSmithBaker[1])
AbsorptionCoeffsPopeFry = numpy.loadtxt("ExistingData/AbsorptionCoefficients_PopeFry.txt", unpack=True)
interpolatedAbsorptionCoeffsPopeFry = scipy.interpolate.interp1d(AbsorptionCoeffsPopeFry[0], AbsorptionCoeffsPopeFry[1],bounds_error=False)

AbsorptionAntaresMeasurement = numpy.loadtxt("ExistingData/ANTARES_measured_absorption.txt", unpack=True)
AbsorptionAntaresMeasurement_Test3_Saclay = numpy.loadtxt("ExistingData/ANTARES_measured_absorption_Test3_Saclay.txt", unpack=True)


def absCoeffGaussians(wavelen, temperature=22.):
    #overtone level(vs,vb) (*,*)     (*,*)    (4,0)   (*,*)   (*,*)    (4,1)   (*,*)    (5,0)   (5,1)   (6,0)   (6,1)   (7,0)   (7,1)   (8,0)    (8,1)
    M =                    [47.48,   23.33,   35.07,  1.794,  9.216,   4.955,  2.341,   3.574,  1.310,  0.3359, 0.2010, 0.1161, 0.0138, 0.03839, 0.2219]
    lambda0 =              [795.,    775.,    744.,   740.,   697.,    669.,   638.,    610.,   558.,   517.,   485.,   449.,   415.,   396.,    370.]
    sigma =                [29.87,   24.79,   20.28,  5.48,   28.22,   24.78,  20.08,   18.40,  22.84,  13.52,  19.27,  18.86,  15.79,  20.88,   21.09]
    MT =                   [-0.0010, -0.0010, 0.0062, 0.0045, -0.0010, 0.0020, -0.0040, 0.0045, 0.0020, 0.0045, 0.0020, 0.0045, 0.0020, 0.0045,  0.0020]
    refTemp = 22.

    resultAbsCoeff = 0.
    resultTempCorr = 0.
    for i in range(len(M)):
        resultAbsCoeff = resultAbsCoeff +       (M[i]/sigma[i])*numpy.exp(-((wavelen-lambda0[i])**2.)/(2.*sigma[i]**2.))
        resultTempCorr = resultTempCorr + (MT[i]*M[i]/sigma[i])*numpy.exp(-((wavelen-lambda0[i])**2.)/(2.*sigma[i]**2.))

    result = resultAbsCoeff + resultTempCorr*(temperature-refTemp)

    return numpy.where(wavelen<380., wavelen*float('NaN'), result)

def absLenIceCubeExample(wavelen):
    alpha=     0.898608505726  #+-   0.027638472617
    kappa=     1.084106802940  #+-   0.014470303431
    A=      6954.090332031250  #+- 973.426452636719
    B=      6617.754394531250  #+-  71.282264709473
    D=        71.402900695801  #+-  12.159952163696
    E=         2.566572427750  #+-   0.584202528000
    
    # for example: 2558.47m
    be400=    0.0413266
    adust400= 0.0676581
    deltat=   20.501200
    
    astar400 = D*adust400+E

    absCoeff = astar400 * (wavelen)**(-kappa) + A*numpy.exp(-B/wavelen)*(1.+0.01*deltat)
    return 1./absCoeff


# interpolated versions (data copied from km3)
interpolatedPhaseRefIndex = scipy.interpolate.InterpolatedUnivariateSpline(wlens, refIndex)
#interpolatedPhaseRefIndex = scipy.interpolate.interp1d(wlens, refIndex)
def interpolatedGroupRefIndex(wlen): # calculation of the group refractive index from the phase refractive index
    interpolatedPhaseRefIndexDerivative = lambda x: scipy.misc.derivative(interpolatedPhaseRefIndex, x)
    np = interpolatedPhaseRefIndex(wlen)
    return np/(1.+(wlen/np)*interpolatedPhaseRefIndexDerivative(wlen))
interpolatedScatLen = scipy.interpolate.InterpolatedUnivariateSpline(wlens, scatLength)
interpolatedAbsLen = scipy.interpolate.interp1d(wlens, absLength, bounds_error=False)

#### Oxford analysis scanned values
OxfordSC_lowerbound = numpy.loadtxt("ExistingData/OxfordSC_lowerbound.txt", unpack=True)
interpolatedOxfordSC_lowerbound = scipy.interpolate.interp1d(OxfordSC_lowerbound[0], OxfordSC_lowerbound[1])
OxfordSC_upperbound = numpy.loadtxt("ExistingData/OxfordSC_upperbound.txt", unpack=True)
interpolatedOxfordSC_upperbound = scipy.interpolate.interp1d(OxfordSC_upperbound[0], OxfordSC_upperbound[1])

OxfordLC_lowerbound = numpy.loadtxt("ExistingData/OxfordLC_lowerbound.txt", unpack=True)
interpolatedOxfordLC_lowerbound = scipy.interpolate.interp1d(OxfordLC_lowerbound[0], OxfordLC_lowerbound[1])
OxfordLC_upperbound = numpy.loadtxt("ExistingData/OxfordLC_upperbound.txt", unpack=True)
interpolatedOxfordLC_upperbound = scipy.interpolate.interp1d(OxfordLC_upperbound[0], OxfordLC_upperbound[1])



# parametric versions
####
# Quan&Fry (taken from W. Schuster's thesis):
refind_S  = 38.44    # salinity in ppt
refind_T  = 13.1     # temperature in degC
refind_P  = 213.0    # ambient pressure [atm]   # 213 bar in comb. with the previous salinity and temp. seems to approximate the km3 tables very closely

refind_n0 = 1.31405  # offset
refind_n1 = 1.45e-5
refind_n2 = 1.779e-4
refind_n3 = 1.05e-6
refind_n4 = 1.6e-8
refind_n5 = 2.02e-6
refind_n6 = 15.868
refind_n7 = 0.01155
refind_n8 = 0.00423
refind_n9 = 4382.
refind_n10 = 1.1455e6
# these get used in the calculation:
refind_a0 = refind_n0+(refind_n2-refind_n3*refind_T+refind_n4*refind_T*refind_T)*refind_S-refind_n5*refind_T*refind_T
refind_a1 = refind_n1
refind_a2 = refind_n6+refind_n7*refind_S-refind_n8*refind_T
refind_a3 = -refind_n9
refind_a4 = refind_n10

def getPhaseRefIndex(wavelength):
    x = 1./wavelength
    return (refind_a0  +  refind_a1*refind_P  +  x*(refind_a2 + x*(refind_a3 + x*refind_a4)))

def getDispersionPhase(wavelength):
    x = 1./wavelength
    return -x*x*(refind_a2 + x*(2.0*refind_a3 + x*3.0*refind_a4))

def getGroupRefIndex(wavelength):
    #c_light = 0.299792458
    n_inv = 1./getPhaseRefIndex(wavelength)
    y = getDispersionPhase(wavelength);
    return 1./((1.0 + y*wavelength*n_inv) * n_inv)

def Cherenkov_dN_dXdwlen(wlen, beta=1.):
    return (2.*math.pi/(137.*(wlen**2.)))*(1. - 1./((beta*getPhaseRefIndex(wlen))**2.))


def getPhaseRefIndex_IceCube(wavelength):
    x = wavelength/1000.# wavelength in micrometer
    return 1.55749 - 1.57988*x + 3.99993*x**2. - 4.68271*x**3. + 2.09354*x**4.

def getGroupRefIndex_IceCube(wavelength):
    x = wavelength/1000.# wavelength in micrometer
    np = getPhaseRefIndex_IceCube(wavelength)
    return np * (1. + 0.227106 - 0.954648*x + 1.42568*x**2. - 0.711832*x**3.)
    
def Cherenkov_dN_dXdwlen_IceCube(wlen, beta=1.):
    return (2.*math.pi/(137.*(wlen**2.)))*(1. - 1./((beta*getPhaseRefIndex_IceCube(wlen))**2.))

print(Cherenkov_dN_dXdwlen_IceCube(470.))

numberOfPhotonsPerNanometer, err = scipy.integrate.quadrature(Cherenkov_dN_dXdwlen, 290., 610.)
#numberOfPhotonsPerNanometer, err = scipy.integrate.quadrature(Cherenkov_dN_dXdwlen_IceCube, 265., 675.)
print(err)
numberOfPhotonsPerMeter = numberOfPhotonsPerNanometer*1e9

print("photons per meter between [290..610]nm =", numberOfPhotonsPerMeter)

def getScatteringLengthSCOxford(wavelength):
    fixedWlen = 374.5
    fixedScatlen = 137.
    exponent = 4.32
    return fixedScatlen * ((wavelength/fixedWlen)**(exponent))

def getScatteringLengthSCOxfordUpper(wavelength):
    fixedWlen = 374.5
    fixedScatlen = 137.+math.sqrt(6.**2. + 16.**2.)
    exponent = 4.32+0.31*numpy.sign(wavelength-fixedWlen)
    return fixedScatlen * ((wavelength/fixedWlen)**(exponent))

def getScatteringLengthSCOxfordLower(wavelength):
    fixedWlen = 374.5
    fixedScatlen = 137.-math.sqrt(6.**2. + 14.**2.)
    exponent = 4.32-0.31*numpy.sign(wavelength-fixedWlen)
    return fixedScatlen * ((wavelength/fixedWlen)**(exponent))

def getScatteringLengthLCOxford(wavelength):
    fixedWlen = 374.5
    fixedScatlen = 173.
    exponent = 1.0
    return fixedScatlen * ((wavelength/fixedWlen)**(exponent))

def getScatteringLengthLCOxfordUpper(wavelength):
    fixedWlen = 374.5
    fixedScatlen = 173.+math.sqrt(10.**2. + 15.**2.)
    exponent = 1.0+0.7*numpy.sign(wavelength-fixedWlen)
    return fixedScatlen * ((wavelength/fixedWlen)**(exponent))

def getScatteringLengthLCOxfordLower(wavelength):
    fixedWlen = 374.5
    fixedScatlen = 173.-math.sqrt(10.**2. + 16.**2.)
    exponent = 1.0-0.7*numpy.sign(wavelength-fixedWlen)
    return fixedScatlen * ((wavelength/fixedWlen)**(exponent))

def getScatteringLengthOxford(wavelength):
    return 1./((1./getScatteringLengthSCOxford(wavelength))+(1./getScatteringLengthLCOxford(wavelength)))

def getEtaOxford():
    return 173./(173.+137.)

#### Oxford Model "B" (fit assuming fixed g_LC==0.92)
def getScatteringLengthSCOxfordModelB(wavelength):
    fixedWlen = 374.5
    fixedScatlen = 115. # +-5
    exponent = 4.32
    return fixedScatlen * ((wavelength/fixedWlen)**(exponent))

def getScatteringLengthLCOxfordModelB(wavelength):
    fixedWlen = 374.5
    fixedScatlen = 47. # +-7
    exponent = 1.0
    return fixedScatlen * ((wavelength/fixedWlen)**(exponent))

def getScatteringLengthOxfordModelB(wavelength):
    return 1./((1./getScatteringLengthSCOxfordModelB(wavelength))+(1./getScatteringLengthLCOxfordModelB(wavelength)))

def getEtaOxfordModelB():
    return 47./(47.+115.)

#### Saclay model

def getScatteringLengthSCSaclay(wavelength):
    fixedWlen = 374.5
    fixedScatlen = 32.46/0.17
    exponent = 4.32
    return fixedScatlen * ((wavelength/fixedWlen)**(exponent))

def getScatteringLengthLCSaclay(wavelength):
    fixedWlen = 374.5
    fixedScatlen = 32.46/(1.-0.17)
    exponent = 1.0
    return fixedScatlen * ((wavelength/fixedWlen)**(exponent))

def getScatteringLengthSaclay(wavelength):
    return 1./((1./getScatteringLengthSCSaclay(wavelength))+(1./getScatteringLengthLCSaclay(wavelength)))


def getScatteringLengthKopelevich(wavelength,
                                  volumeConcentrationSmallParticles=0.0075,   # in ppm
                                  volumeConcentrationLargeParticles=0.0075):  # in ppm
    refWlen = 550.
    
    x = refWlen/wavelength
    scatCoeff = 0.0017 * x**4.3 + 1.34 * volumeConcentrationSmallParticles * x**1.7 + 0.312 * volumeConcentrationLargeParticles * x**0.3
    
    return 1./scatCoeff


#print getScatteringLengthSaclay(374.5)
#print interpolatedScatLen(374.5)

####


fig = pylab.figure(3)
fig.subplots_adjust(left=0.09, bottom=0.05, top=0.95, right=0.98)

ax = fig.add_subplot(3, 1, 1)
bx = fig.add_subplot(3, 1, 2)
cx = fig.add_subplot(3, 1, 3)

wlens=numpy.linspace(290.,610.,num=10000)

ax.plot(wlens, interpolatedPhaseRefIndex(wlens), linewidth=3., color='k', linestyle='solid', label="$n_\\mathrm{phase}$ (km3 table, spline interp.)")
l, = ax.plot(wlens, interpolatedGroupRefIndex(wlens), linewidth=3., color='k', label="$n_\\mathrm{group}$ (km3 table, spline interp.)")
l.set_dashes([5,5])

ax.plot(wlens, getPhaseRefIndex(wlens), linewidth=1., color='r', linestyle='solid', label=r"$n_\mathrm{phase}$ (Quan\&Fry @ $%.0f\mathrm{bar}$)" % (refind_P))
l, = ax.plot(wlens, getGroupRefIndex(wlens), linewidth=1., color='r', label=r"$n_\mathrm{group}$ (Quan\&Fry @ $%.0f\mathrm{bar}$)" % (refind_P))
l.set_dashes([5,5])


ax.plot(wlens, getPhaseRefIndex_IceCube(wlens), linewidth=1., color='b', linestyle='solid', label=r"$n_\mathrm{phase}$ (IceCube)")
l, = ax.plot(wlens, getGroupRefIndex_IceCube(wlens), linewidth=1., color='b', label=r"$n_\mathrm{group}$ (IceCube)")
l.set_dashes([5,5])


ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(20))
ax.set_xlim(290.,610.)
ax.legend()
ax.grid(True)
ax.set_xlabel("wavelength $\\lambda [\\mathrm{nm}]$")
ax.set_ylabel("refractive index $n$")



bx.plot(wlens, 1./interpolatedAbsorptionCoeffsGeasim(wlens), linewidth=2, color='k', label="table from geasim (used in km3)")
bx.plot(wlens, interpolatedAbsLen(wlens), linewidth=2, color='k', label="very old table, used to be in km3")
bx.plot(wlens, 1./interpolatedAbsorptionCoeffsSmithBaker(wlens), linewidth=2, linestyle='-', color='0.5', label=r"Smith\&Baker (pure water)")
bx.plot(wlens, absLenIceCubeExample(wlens), linewidth=1, linestyle='-', color='r', label=r"IceCube")
bx.plot(wlens, 1./interpolatedAbsorptionCoeffsPopeFry(wlens), linewidth=1, linestyle='-', color='k', label=r"Pope\&Fry (pure water)")
bx.plot(wlens, 1./absCoeffGaussians(wlens), linewidth=1, linestyle='-', color='b', label=r"gaussians")


bx.errorbar(AbsorptionAntaresMeasurement[0],
            1./AbsorptionAntaresMeasurement[1], 
            yerr=[(1./AbsorptionAntaresMeasurement[1]-1./(AbsorptionAntaresMeasurement[1]-AbsorptionAntaresMeasurement[2])),
                  1./(AbsorptionAntaresMeasurement[1]+AbsorptionAntaresMeasurement[2])-1./AbsorptionAntaresMeasurement[1]],
            fmt='x', color='b',
            label="Antares site measurements (July 2002)")
bx.errorbar(AbsorptionAntaresMeasurement_Test3_Saclay[0],
            AbsorptionAntaresMeasurement_Test3_Saclay[1],
            yerr=AbsorptionAntaresMeasurement_Test3_Saclay[2],
            fmt='x', color='r',
            label="Test 3' analysis (Saclay)")
bx.errorbar([374.5],
            [25.9],
            yerr=[[math.sqrt(0.5**2. + 1.1**2.)], [math.sqrt(0.5**2. + 1.3**2.)]],
            fmt='x', color='g',
            label="Test 3' analysis (Oxford)")


bx.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(20))
bx.set_xlim(290.,610.)
#bx.set_ylim(0.,100.)
bx.legend(loc="upper right")
bx.grid(True)
bx.set_xlabel("wavelength $\\lambda [\\mathrm{nm}]$")
bx.set_ylabel("absorption length $[\\mathrm{m}]$")

# the table is the same as the kopelevich model with nu_s==nu_l==0.0075
cx.plot(wlens, interpolatedScatLen(wlens), linewidth=3, color='k', label="used in km3")
cx.plot(wlens, 
        getScatteringLengthKopelevich(wlens,
                                      volumeConcentrationSmallParticles=0.0075,
                                      volumeConcentrationLargeParticles=0.0075), 
        linewidth=1, color='r', label=r"Kopelevich $\nu_s=0.0075, \nu_l=0.0075$")

wlen_blue = 473.0
wlen_UV = 374.5
#print "km3 @ blue = ", interpolatedScatLen(wlen_blue)
#print "km3 @ UV   = ", interpolatedScatLen(wlen_UV)

#cx.plot(wlens, getScatteringLengthSCOxford(wlens), linewidth=2, color='r', label="Oxford SC")
#cx.fill_between(wlens, getScatteringLengthSCOxfordLower(wlens), getScatteringLengthSCOxfordUpper(wlens), linewidth=1, linestyle='-', color=(1.0,0.5,0.5), label="Oxford SC")
#
#cx.plot(wlens, getScatteringLengthLCOxford(wlens), linewidth=2, color='g', label="Oxford LC")
#cx.fill_between(wlens, getScatteringLengthLCOxfordLower(wlens), getScatteringLengthLCOxfordUpper(wlens), linewidth=1, linestyle='-', color=(0.5,1.0,0.5), label="Oxford LC")

cx.plot(wlens, getScatteringLengthOxford(wlens), linewidth=2, linestyle='-',  color='b', label=r"Oxford model (``Model A'') ($\eta_{374.5\mathrm{nm}}=%4.2f$)" % (getEtaOxford()))
cx.plot(wlens, getScatteringLengthOxfordModelB(wlens), linewidth=2, linestyle='--', color='b', label=r"Oxford model (``Model B'') ($\eta_{374.5\mathrm{nm}}=%4.2f$)" % (getEtaOxfordModelB()))


cx.plot(wlens, getScatteringLengthSaclay(wlens), linewidth=2, color='m', label=r"Saclay model ($\eta_{374.5\mathrm{nm}}=0.17$)")





#wlens2=numpy.linspace(354.,547.,num=10000)
#cx.fill_between(wlens2, interpolatedOxfordSC_lowerbound(wlens2), interpolatedOxfordSC_upperbound(wlens2), linewidth=1, linestyle='-', color=(1.0,0.5,0.5), label="Oxford SC (scanned)")
#cx.fill_between(wlens2, interpolatedOxfordLC_lowerbound(wlens2), interpolatedOxfordLC_upperbound(wlens2), linewidth=1, linestyle='-', color=(1.0,0.5,0.5), label="Oxford LC (scanned)")



cx.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(20))
cx.set_xlim(290.,610.)
cx.grid(True)
cx.legend(loc='upper left')
cx.set_xlabel(r"wavelength $\lambda [\mathrm{nm}]$")
cx.set_ylabel(r"scattering length $\lambda_\mathrm{scat;geom}$ $[\mathrm{m}]$")

pylab.savefig("water_properties.pdf", transparent=True)

