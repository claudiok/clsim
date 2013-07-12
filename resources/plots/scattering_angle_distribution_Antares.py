#!/usr/bin/env python

from __future__ import print_function

import math
import numpy
import random

import matplotlib
matplotlib.use("PDF")

prng = numpy.random.RandomState()


fig_size = [11.7,8.3] # din A4
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

deg=1.
data_ang = [0.100*deg,   0.126*deg,   0.158*deg,   0.200*deg,   0.251*deg,
        0.316*deg,   0.398*deg,   0.501*deg,   0.631*deg,   0.794*deg,
        1.000*deg,   1.259*deg,   1.585*deg,   1.995*deg,   2.512*deg,
        3.162*deg,   3.981*deg,   5.012*deg,   6.310*deg,   7.943*deg,
        10.000*deg,  15.000*deg,  20.000*deg,  25.000*deg,  30.000*deg,
        35.000*deg,  40.000*deg,  45.000*deg,  50.000*deg,  55.000*deg,
        60.000*deg,  65.000*deg,  70.000*deg,  75.000*deg,  80.000*deg,
        85.000*deg,  90.000*deg,  95.000*deg, 100.000*deg, 105.000*deg,
        110.000*deg, 115.000*deg, 120.000*deg, 125.000*deg, 130.000*deg,
        135.000*deg, 140.000*deg, 145.000*deg, 150.000*deg, 155.000*deg,
        160.000*deg, 165.000*deg, 170.000*deg, 175.000*deg, 180.000*deg]

data_val=[1.767E+03, 1.296E+03, 9.502E+02, 6.991E+02, 5.140E+02,
        3.764E+02, 2.763E+02, 2.188E+02, 1.444E+02, 1.022E+02,
        7.161E+01, 4.958E+01, 3.395E+01, 2.281E+01, 1.516E+01,
        1.002E+01, 6.580E+00, 4.295E+00, 2.807E+00, 1.819E+00,
        1.153E+00, 4.893E-01, 2.444E-01, 1.472E-01, 8.609E-02,
        5.931E-02, 4.210E-02, 3.067E-02, 2.275E-02, 1.699E-02,
        1.313E-02, 1.046E-02, 8.488E-03, 6.976E-03, 5.842E-03,
        4.953E-03, 4.292E-03, 3.782E-03, 3.404E-03, 3.116E-03,
        2.912E-03, 2.797E-03, 2.686E-03, 2.571E-03, 2.476E-03,
        2.377E-03, 2.329E-03, 2.313E-03, 2.365E-03, 2.506E-03,
        2.662E-03, 2.835E-03, 3.031E-03, 3.092E-03, 3.154E-03]

PhaseFunctionKopelevic = numpy.loadtxt("ExistingData/PhaseFunctionKopelevicTwoComponent.txt", unpack=True)
#interpolatedPhaseFunctionKopelevicSmall = scipy.interpolate.interp1d(PhaseFunctionKopelevic[0], PhaseFunctionKopelevic[1])
#interpolatedPhaseFunctionKopelevicLarge = scipy.interpolate.interp1d(PhaseFunctionKopelevic[0], PhaseFunctionKopelevic[2])

def interpolatedPhaseFunctionKopelevic(angle, wavelength):
    pSmall = 0.0075
    pLarge = 0.0075
    factorSmall = pSmall*(wavelength/550.)**(-1.7)
    factorLarge = pLarge*(wavelength/550.)**(-0.3)
    phaseFunctionValues = factorSmall*PhaseFunctionKopelevic[1] + factorLarge*PhaseFunctionKopelevic[2]
    
    return numpy.exp(scipy.interpolate.interp1d(numpy.log(PhaseFunctionKopelevic[0]), numpy.log(phaseFunctionValues),bounds_error=False)(numpy.log(angle)))


factor_p00075 = 0.17
factor_Oxford = 0.55

def prepareTables(data_ang, data_val):
    if len(data_ang) != len(data_val):
        raise Exception
        
    # Convert to radians part_ang and
    # normalize part_beta according to 3.8 (Light and Water, Mobley)
    partic_ang = []
    partic_beta = []
    for i in range(len(data_ang)):
        angleInRad = data_ang[i]*math.pi/180.
        partic_ang.append(angleInRad)
        partic_beta.append(2.*math.pi*data_val[i]*math.sin(angleInRad))
    
    # Compute first value (1e-9, practically zero) assuming
    # beta = theta**m (m = -1.346)
    partic_ang_m1 = 1e-9
    partic_beta_m1 = partic_beta[0] * ((partic_ang_m1/partic_ang[0]) ** (-1.346))
    partic_beta_m1 = 2. * math.pi * partic_beta_m1 * math.sin(partic_ang_m1);
    
    
    # Integrate angular distribution (trapezoidal rule)
    # int = h*(f0+f1)/2
    partic_acu = numpy.zeros(len(partic_ang))
    partic_acu[0] = (partic_ang[0]-partic_ang_m1)*(partic_beta_m1+partic_beta[0])/2.
    for j in range(1,len(partic_ang)):
        partic_acu[j] = partic_acu[j-1] + (partic_ang[j]-partic_ang[j-1]) * (partic_beta[j]+partic_beta[j-1]) / 2.;
    
    # Normalize
    for j in range(len(partic_ang)):
        partic_acu[j] = partic_acu[j]/partic_acu[len(partic_ang)-1];
    
    return (partic_ang, partic_acu)

partic_ang, partic_acu = prepareTables(data_ang, data_val)

#print "partic_acu", partic_acu
#print "partic_ang", partic_ang


def generateParticCosAngle(randomVariable, partic_ang, partic_acu):
    r = randomVariable
    
    k=0
    while (r > partic_acu[k]):
        k+=1;
    
    if k==0:
        angulo = r * partic_ang[0] / partic_acu[0]
    else:
        angulo = partic_ang[k-1] + (r-partic_acu[k-1])* \
        (partic_ang[k] - partic_ang[k-1] )/(partic_acu[k] - partic_acu[k-1])
    
    return math.cos(angulo)

def generateRayleighCosAngle(randomVariable):
    r = randomVariable
    
    b = 0.835
    p = 1./b
    q = (b+3.)*(r-0.5)/b
    d = q*q + p*p*p
    u1 = -q+math.sqrt(d)
    u = (abs(u1))**(1./3.)
    if u1<0.: u = -u
    
    v1 = -q-math.sqrt(d)
    v = (abs(v1))**(1./3.)
    if v1 < 0.: v = -v
    
    coscorr = max(-1., min( 1., u+v))
    return coscorr
    
def generatep00075CosAngle(randomVariable, randomVariable2, partic_ang, partic_acu):
    if randomVariable2 < factor_p00075:
        return generateRayleighCosAngle(randomVariable)
    else:
        return generateParticCosAngle(randomVariable, partic_ang, partic_acu)


def generateHenyeyGreensteinCosAngle(randomVariable, g):
    r = randomVariable
    wT = g # g==<cos theta>
    
    cosa = -1.0 + 2.0 * ((1.0-wT+wT*r)*r*((1.0+wT)**2.) / ((1.0-wT+2.0*wT*r)**2.));

    return cosa


def generateOxfordCosAngle(randomVariable, randomVariable2):
    if randomVariable2 < factor_Oxford:
        return generateRayleighCosAngle(randomVariable)
    else:
        return generateHenyeyGreensteinCosAngle(randomVariable, 0.55)

def particMobley(cosAngle):
    angle = numpy.arccos(cosAngle)
    ang0 = data_ang[0]*math.pi/180.
    if angle < ang0:
        beta0 = data_val[0]
        return beta0 * ((angle/ang0)**(-1.346))
    else:
        return scipy.interpolate.interp1d(numpy.array(data_ang)*math.pi/180., numpy.array(data_val))(angle)
particMobley = numpy.vectorize(particMobley)
    

def rayleigh(cosAngle):
    a0 = 0.06225
    a1 = 0.835
    norm = 2.*(1.+a1/3.)
    return (1. + a1*cosAngle*cosAngle) / norm
    #return a0 * (x + a1*0.333*x*x*x) 

def henyeyGreenstein(cosAngle, g):
    if g != 0.0:
        return (1.0 - g**2.) / ((1.0 + g**2. - 2.0*g*cosAngle)**(3./2.) * 4.*math.pi);
    else:
        return cosAngle / (4.*math.pi);
henyeyGreenstein = numpy.vectorize(henyeyGreenstein)


def genMCHistograms(generator, samples=10240*1000, numBins=1000):
    angles = []
    for i in range(samples):
        angles.append(generator())
    angles = numpy.array(angles) # convert to numpy array
    numAng_orig, binsAng = scipy.histogram(numpy.arccos(angles)*180./math.pi, range=(0.,180.), bins=numBins)
    numCos_orig, binsCos = scipy.histogram(angles, range=(-1.,1.), bins=numBins)
    del angles # not needed anymore

    numAng=[]
    for i, number in enumerate(numAng_orig):
        binWidth = math.cos(binsAng[i]*math.pi/180.) - math.cos(binsAng[i+1]*math.pi/180.)
        numAng.append(float(number)/float(samples)/binWidth)
    numAng=numpy.array(numAng)

    numCos=[]
    for i, number in enumerate(numCos_orig):
        numCos.append(float(number)/float(samples)/float(2./float(numBins)))
    numCos=numpy.array(numCos)

    binsAng = numpy.array(binsAng[:-1])+(binsAng[1]-binsAng[0])/2.
    binsCos = numpy.array(binsCos[:-1])+(binsCos[1]-binsCos[0])/2.
    
    return dict(cos=dict(num=numCos, bins=binsCos), ang=dict(num=numAng, bins=binsAng))

#
from icecube import icetray, dataclasses, clsim, phys_services

def genMCHistogramsOpenCL(distribution, rng, iterations=1000, numBins=10000):
    # get OpenCL CPU devices
    openCLDevices = [device for device in clsim.I3CLSimOpenCLDevice.GetAllDevices() if device.cpu]
    if len(openCLDevices)==0:
        raise RuntimeError("No CPU OpenCL devices available!")
    openCLDevice = openCLDevices[0]

    openCLDevice.useNativeMath=False
    workgroupSize = 1
    workItemsPerIteration = 10240
    print("           using platform:", openCLDevice.platform)
    print("             using device:", openCLDevice.device)
    print("            workgroupSize:", workgroupSize)
    print("    workItemsPerIteration:", workItemsPerIteration)
    
    tester = clsim.I3CLSimRandomDistributionTester(device=openCLDevice,
                                                   workgroupSize=workgroupSize,
                                                   workItemsPerIteration=workItemsPerIteration,
                                                   randomService=rng,
                                                   randomDistribution=distribution)
    
    angles = tester.GenerateRandomNumbers(iterations)
    samples = len(angles)
    
    print("generated")
    
    angles = numpy.array(angles) # convert to numpy array
    print("converted")
    
    numAng_orig, binsAng = scipy.histogram(numpy.arccos(angles)*(180./math.pi), range=(0.,180.), bins=numBins)
    print("hist1 complete")

    numCos_orig, binsCos = scipy.histogram(angles, range=(-1.,1.), bins=numBins)
    print("hist2 complete")
    
    del angles # not needed anymore
    print("deleted")
    
    numAng=[]
    for i, number in enumerate(numAng_orig):
        binWidth = math.cos(binsAng[i]*math.pi/180.) - math.cos(binsAng[i+1]*math.pi/180.)
        numAng.append(float(number)/float(samples)/binWidth)
    numAng=numpy.array(numAng)
    
    numCos=[]
    for i, number in enumerate(numCos_orig):
        numCos.append(float(number)/float(samples)/float(2./float(numBins)))
    numCos=numpy.array(numCos)
    
    binsAng = numpy.array(binsAng[:-1])+(binsAng[1]-binsAng[0])/2.
    binsCos = numpy.array(binsCos[:-1])+(binsCos[1]-binsCos[0])/2.
    
    return dict(cos=dict(num=numCos, bins=binsCos), ang=dict(num=numAng, bins=binsAng))

def genMCHistogramsHost(distribution, rng, iterations=10000000, numBins=1000):
    print("generating (host)")

    angles = []
    for i in range(iterations):
        angles.append(distribution.SampleFromDistribution(rng, []))
    samples = len(angles)

    print("generated (host)")
    
    angles = numpy.array(angles) # convert to numpy array
    print("converted (host)")
    
    numAng_orig, binsAng = scipy.histogram(numpy.arccos(angles)*(180./math.pi), range=(0.,180.), bins=numBins)
    print("hist1 complete (host)")

    numCos_orig, binsCos = scipy.histogram(angles, range=(-1.,1.), bins=numBins)
    print("hist2 complete (host)")
    
    del angles # not needed anymore
    print("deleted (host)")
    
    numAng=[]
    for i, number in enumerate(numAng_orig):
        binWidth = math.cos(binsAng[i]*math.pi/180.) - math.cos(binsAng[i+1]*math.pi/180.)
        numAng.append(float(number)/float(samples)/binWidth)
    numAng=numpy.array(numAng)
    
    numCos=[]
    for i, number in enumerate(numCos_orig):
        numCos.append(float(number)/float(samples)/float(2./float(numBins)))
    numCos=numpy.array(numCos)
    
    binsAng = numpy.array(binsAng[:-1])+(binsAng[1]-binsAng[0])/2.
    binsCos = numpy.array(binsCos[:-1])+(binsCos[1]-binsCos[0])/2.
    
    return dict(cos=dict(num=numCos, bins=binsCos), ang=dict(num=numAng, bins=binsAng))


rng = phys_services.I3SPRNGRandomService(seed=3244, nstreams=2, streamnum=0)

#print clsim.GetPetzoldScatteringCosAngleDistribution().GetOpenCLFunction("func", "ARGS", "ARGSTOCALL", "CO", "OC")

#hist_p00075 = genMCHistograms(lambda: generatep00075CosAngle(prng.uniform(0.,1.), prng.uniform(0.,1.), partic_ang, partic_acu))
#hist_Oxford = genMCHistograms(lambda: generateOxfordCosAngle(prng.uniform(0.,1.), prng.uniform(0.,1.)))

hist_p00075 = genMCHistogramsOpenCL(clsim.GetAntaresScatteringCosAngleDistribution(), rng)
#hist_p00075 = genMCHistogramsOpenCL(clsim.GetPetzoldScatteringCosAngleDistribution(), rng)
#hist_p00075 = genMCHistogramsOpenCL(clsim.I3CLSimRandomValueRayleighScatteringCosAngle(), rng)

hist_p00075_host = genMCHistogramsHost(clsim.GetAntaresScatteringCosAngleDistribution(), rng)



fig = pylab.figure(1)
fig.subplots_adjust(left=0.06, bottom=0.06, top=0.98, right=0.98)

ax = fig.add_subplot(2, 2, 1)
bx = fig.add_subplot(2, 2, 2)
cx = fig.add_subplot(2, 2, 3)
dx = fig.add_subplot(2, 2, 4)

HG_cosTheta = 0.924
#HG_cosTheta = 0.8625


if True:
    fineBins = numpy.logspace(numpy.log10(0.1), numpy.log10(180.), num=1000, base=10.)

    ax.semilogy(hist_p00075["ang"]["bins"], hist_p00075["ang"]["num"], linewidth=2, color='r', label="MC generated (OpenCL)")
    ax.semilogy(hist_p00075_host["ang"]["bins"], hist_p00075_host["ang"]["num"], linewidth=2, color='y', label="MC generated (C++/CPU)")

    ax.semilogy(fineBins, particMobley(numpy.cos(numpy.array(fineBins)*math.pi/180.))*2.*math.pi, linewidth=2, color='k', label=r"\textbf{Petzold} (``avg. part.'') (c.f. Mobley et al., 1993) (from km3)")
    ax.semilogy(fineBins, rayleigh(numpy.cos(numpy.array(fineBins)*math.pi/180.)), linewidth=2, color='g', label=r"\textbf{``Rayleigh''} (c.f. Morel et al., 1974) $(\propto 1+0.835 \cos^2 \theta)$")

    ax.semilogy(fineBins, factor_p00075 * rayleigh(numpy.cos(numpy.array(fineBins)*math.pi/180.))  +  (1.-factor_p00075)*particMobley(numpy.cos(numpy.array(fineBins)*math.pi/180.))*2.*math.pi, linewidth=2, color='b', label=r"\textbf{p0.0075 model} ($\eta = %s$)" % (factor_p00075))

    integral = scipy.integrate.romberg(lambda x: interpolatedPhaseFunctionKopelevic(numpy.arccos(numpy.array(x))*180./math.pi, 374.5), -1., 0.9999)
    ax.semilogy(fineBins, interpolatedPhaseFunctionKopelevic(numpy.array(fineBins), 374.5)/integral, linewidth=2, color='0.5', label=r"\textbf{Kopelevic} $\nu_s=\nu_l=0.0075$")

    ax.loglog(fineBins, henyeyGreenstein(numpy.cos(numpy.array(fineBins)*math.pi/180.), HG_cosTheta)*2.*math.pi, linewidth=2, color='m', label=r"\textbf{Henyey-Greenstein} $\left<\cos\theta\right>=%5.3f$" % (HG_cosTheta))

    ax.set_xlabel("scattering angle $\\theta [^\\circ]$")
    ax.set_ylabel("$\\beta(\\theta)$")
    ax.grid(True)
    #ax.legend(loc="upper right")
    ax.set_xlim(0.1,180.)

if True:
    fineBins = numpy.linspace(-1., 1., num=1000)

    bx.semilogy(hist_p00075["cos"]["bins"], hist_p00075["cos"]["num"], linewidth=2, color='r', label="MC generated (OpenCL)")
    bx.semilogy(hist_p00075_host["cos"]["bins"], hist_p00075_host["cos"]["num"], linewidth=2, color='y', label="MC generated (C++/CPU)")

    bx.semilogy(fineBins, particMobley(numpy.array(fineBins))*2.*math.pi, linewidth=2, color='k', label=r"\textbf{Petzold} (``avg. part.'') (c.f. Mobley et al., 1993) (from km3)")
    
    
    bx.semilogy(fineBins, rayleigh(numpy.array(fineBins)), linewidth=2, color='g', label=r"\textbf{``Rayleigh''} (c.f. Morel et al., 1974) $(\propto 1+0.835 \cos^2 \theta)$")

    bx.semilogy(fineBins, factor_p00075 * rayleigh(numpy.array(fineBins))  +  (1.-factor_p00075)*particMobley(numpy.array(fineBins))*2.*math.pi, linewidth=2, color='b', label=r"\textbf{p0.0075 model} ($\eta = %s$)" % (factor_p00075))

    integral = scipy.integrate.romberg(lambda x: interpolatedPhaseFunctionKopelevic(numpy.arccos(numpy.array(x))*180./math.pi, 374.5), -1., 0.9999999)
    kopelevicMeanCosTheta = scipy.integrate.romberg(lambda x: x*interpolatedPhaseFunctionKopelevic(numpy.arccos(x)*180./math.pi, 374.5)/integral, -1.0, 0.9999999)
    bx.semilogy(fineBins, interpolatedPhaseFunctionKopelevic(numpy.arccos(numpy.array(fineBins))*180./math.pi, 374.5)/integral, linewidth=2, color='0.5', label=r"\textbf{Kopelevic} (particles only) $\nu_s=\nu_l=0.0075$ $\Rightarrow$ $\left<\cos\theta\right>=%5.3f$" % (kopelevicMeanCosTheta))

    bx.semilogy(fineBins, henyeyGreenstein(numpy.array(fineBins), HG_cosTheta)*2.*math.pi, linewidth=2, color='m', label=r"\textbf{Henyey-Greenstein} $\left<\cos\theta\right>=%5.3f$" % (HG_cosTheta))
    
    bx.set_xlabel(r"scattering angle $\cos\theta [^\circ]$")
    bx.set_ylabel("$\\beta(\\theta)$")
    bx.grid(True)
    bx.legend(loc="upper left")

if True:
    fineBins = numpy.logspace(numpy.log10(0.1), numpy.log10(180.), num=1000, base=10.)

    #cx.semilogy(hist_Oxford["ang"]["bins"], hist_Oxford["ang"]["num"], linewidth=2, color='r', label="MC generated")

    cx.loglog(fineBins, henyeyGreenstein(numpy.cos(numpy.array(fineBins)*math.pi/180.), 0.54)*2.*math.pi, linewidth=2, color='k', label=r"\textbf{Henyey-Greenstein} (``LC'') $\left<\cos\theta\right>=0.54$")
    cx.loglog(fineBins, rayleigh(numpy.cos(numpy.array(fineBins)*math.pi/180.)), linewidth=2, color='g', label=r"\textbf{``Rayleigh''} (``SC'') $(\propto 1+0.835 \cos^2 \theta)$")

    cx.loglog(fineBins, factor_Oxford * rayleigh(numpy.cos(numpy.array(fineBins)*math.pi/180.))  +  (1.-factor_Oxford)*henyeyGreenstein(numpy.cos(numpy.array(fineBins)*math.pi/180.), 0.54)*2.*math.pi, linewidth=2, color='b', label=r"\textbf{Oxford model} ($\eta = %s$)" % (factor_Oxford))
    cx.loglog(fineBins, factor_p00075 * rayleigh(numpy.cos(numpy.array(fineBins)*math.pi/180.))  +  (1.-factor_p00075)*particMobley(numpy.cos(numpy.array(fineBins)*math.pi/180.))*2.*math.pi,           linewidth=1, color='b', linestyle='--', label=r"\textbf{p0.0075 model} ($\eta = %s$)" % (factor_p00075))


    cx.set_xlabel("scattering angle $\\theta [^\\circ]$")
    cx.set_ylabel("$\\beta(\\theta)$")
    cx.grid(True)
    cx.legend(loc="upper right")
    cx.set_xlim(0.1,180.)
    cx.set_ylim(1e-2,1e5)

if True:
    fineBins = numpy.linspace(-1., 1.0, num=1000)

    #dx.semilogy(hist_Oxford["cos"]["bins"], hist_Oxford["cos"]["num"], linewidth=2, color='r', label="MC generated")

    dx.semilogy(fineBins, henyeyGreenstein(numpy.array(fineBins), 0.54)*2.*math.pi, linewidth=2, color='k', label=r"\textbf{Henyey-Greenstein} (``LC'') $\left<\cos\theta\right>=0.54$")
    dx.semilogy(fineBins, rayleigh(numpy.array(fineBins)), linewidth=2, color='g', label=r"\textbf{``Rayleigh''} (``SC'') $(\propto 1+0.835 \cos^2 \theta)$")

    dx.semilogy(fineBins, factor_Oxford * rayleigh(numpy.array(fineBins))  +  (1.-factor_Oxford)*henyeyGreenstein(numpy.array(fineBins), 0.54)*2.*math.pi, linewidth=2, color='b', label=r"\textbf{Oxford model} ($\eta = %s$)" % (factor_Oxford))
    dx.semilogy(fineBins, factor_p00075 * rayleigh(numpy.array(fineBins))  +  (1.-factor_p00075)*particMobley(numpy.array(fineBins))*2.*math.pi,           linewidth=1, color='b', linestyle='--', label=r"\textbf{p0.0075 model} ($\eta = %s$)" % (factor_p00075))


    dx.set_xlabel(r"scattering angle $\cos\theta [^\circ]$")
    dx.set_ylabel(r"$\beta(\theta)$")
    dx.grid(True)
    dx.legend(loc="upper left")
    dx.set_xlim(-1.,1.)
    dx.set_ylim(1e-2,1e3)



pylab.savefig("scattering_angle_distribution_Antares.pdf", transparent=False)
