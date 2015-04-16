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


def _henyeyGreenstein(cosAngle, g):
    if g != 0.0:
        return 0.5 * ( (1.0 - g**2.) / ((1.0 + g**2. - 2.0*g*cosAngle)**(3./2.)) )
    else:
        return 0.5 * cosAngle
henyeyGreenstein = numpy.vectorize(_henyeyGreenstein)

def _simplifiedLiu(cosAngle, g):
    alpha = 2.*g/(1.-g)

    #print "alpha =", alpha
    norm = (alpha+1.)/(2.**(alpha+1.))
    #print "norm =", norm
    
    return norm*((1.+cosAngle)**alpha)
simplifiedLiu = numpy.vectorize(_simplifiedLiu)

#integral = scipy.integrate.quadrature(lambda x: simplifiedLiu(x, 0.943), -1., 1., maxiter=200)
#print "integral =", integral[0]


from icecube import icetray, dataclasses, clsim, phys_services

rng = phys_services.I3SPRNGRandomService(seed=3245, nstreams=2, streamnum=0)

def genMCHistogramsOpenCL(distribution, rng, iterations=100, numBins=1000):
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
    
    print("maxWorkgroupSizeForKernel:", tester.maxWorkgroupSize)
    
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

#print clsim.GetIceCubeScatteringCosAngleDistribution().GetOpenCLFunction("func", "ARGS", "ARGSTOCALL", "CO", "OC")

mediumProps = clsim.MakeIceCubeMediumProperties()
hist_cosangles = genMCHistogramsOpenCL(mediumProps.ScatteringCosAngleDistribution, rng)
hist_cosangles_host = genMCHistogramsHost(mediumProps.ScatteringCosAngleDistribution, rng)






fig = pylab.figure(1)
fig.subplots_adjust(left=0.06, bottom=0.06, top=0.98, right=0.98)

ax = fig.add_subplot(2, 2, 1)
bx = fig.add_subplot(2, 2, 2)
cx = fig.add_subplot(2, 2, 3)
dx = fig.add_subplot(2, 2, 4)



mean_cosTheta = 0.943       # closer to reality
mean_cosTheta_inSim = 0.9   # used in simulation

factor_SPICE_Mie = 0.45


if True:
    fineBins = numpy.logspace(numpy.log10(0.1), numpy.log10(180.), num=1000, base=10.)

    ax.loglog(fineBins, henyeyGreenstein(numpy.cos(numpy.array(fineBins)*math.pi/180.), mean_cosTheta), linewidth=2, color='m', label=r"\textbf{Henyey-Greenstein (HG)} $\left<\cos\theta\right>=%5.3f$" % (mean_cosTheta))

    ax.loglog(fineBins, simplifiedLiu(numpy.cos(numpy.array(fineBins)*math.pi/180.), mean_cosTheta), linewidth=2, color='g', label=r"\textbf{simplified Liu (sL)} $\left<\cos\theta\right>=%5.3f$" % (mean_cosTheta))

    ax.loglog(fineBins, factor_SPICE_Mie * simplifiedLiu(numpy.cos(fineBins*math.pi/180.), mean_cosTheta)  +  (1.-factor_SPICE_Mie)*henyeyGreenstein(numpy.cos(fineBins*math.pi/180.), mean_cosTheta), linewidth=2, color='b', label=r"$%4.2f \times \textbf{HG} + %4.2f \times \textbf{sL}$ $\left<\cos\theta\right>=%5.3f$" % (1.-factor_SPICE_Mie, factor_SPICE_Mie, mean_cosTheta))

    ax.loglog(fineBins, factor_SPICE_Mie * simplifiedLiu(numpy.cos(fineBins*math.pi/180.), mean_cosTheta_inSim)  +  (1.-factor_SPICE_Mie)*henyeyGreenstein(numpy.cos(fineBins*math.pi/180.), mean_cosTheta_inSim), linewidth=1, linestyle='--', color='b', label=r"$%4.2f \times \textbf{HG} + %4.2f \times \textbf{sL}$ $\left<\cos\theta\right>=%5.3f$" % (1.-factor_SPICE_Mie, factor_SPICE_Mie, mean_cosTheta_inSim))

    ax.set_xlabel("scattering angle $\\theta [^\\circ]$")
    ax.set_ylabel("$\\beta(\\theta)$")
    ax.grid(True)
    #ax.legend(loc="upper right")
    ax.set_xlim(0.1,180.)
    ax.set_ylim(1e-3,1e3)

if True:
    fineBins = numpy.linspace(-1., 1., num=1000)

    bx.semilogy(fineBins, henyeyGreenstein(fineBins, mean_cosTheta), linewidth=2, color='m', label=r"\textbf{Henyey-Greenstein (HG)} $\left<\cos\theta\right>=%5.3f$" % (mean_cosTheta))

    bx.semilogy(fineBins, simplifiedLiu(fineBins, mean_cosTheta), linewidth=2, color='g', label=r"\textbf{simplified Liu (sL)} $\left<\cos\theta\right>=%5.3f$" % (mean_cosTheta))

    bx.semilogy(fineBins, factor_SPICE_Mie * simplifiedLiu(fineBins, mean_cosTheta)  +  (1.-factor_SPICE_Mie)*henyeyGreenstein(fineBins, mean_cosTheta), linewidth=2, color='b', label=r"$%4.2f \times \textbf{HG} + %4.2f \times \textbf{sL}$ $\left<\cos\theta\right>=%5.3f$" % (1.-factor_SPICE_Mie, factor_SPICE_Mie, mean_cosTheta))

    bx.semilogy(fineBins, factor_SPICE_Mie * simplifiedLiu(fineBins, mean_cosTheta_inSim)  +  (1.-factor_SPICE_Mie)*henyeyGreenstein(fineBins, mean_cosTheta_inSim), linewidth=1, linestyle='--', color='b', label=r"$%4.2f \times \textbf{HG} + %4.2f \times \textbf{sL}$ $\left<\cos\theta\right>=%5.3f$" % (1.-factor_SPICE_Mie, factor_SPICE_Mie, mean_cosTheta_inSim))

    bx.set_xlabel(r"scattering angle $\cos\theta [^\circ]$")
    bx.set_ylabel("$\\beta(\\theta)$")
    bx.grid(True)
    bx.legend(loc="upper left")
    bx.set_ylim(1e-3,1e3)

if True:
    fineBins = numpy.logspace(numpy.log10(0.1), numpy.log10(180.), num=1000, base=10.)
    
    cx.loglog(hist_cosangles["ang"]["bins"], hist_cosangles["ang"]["num"], linewidth=2, color='r', label="MC generated (OpenCL)")
    cx.loglog(hist_cosangles_host["ang"]["bins"], hist_cosangles_host["ang"]["num"], linewidth=2, color='y', label="MC generated (C++/CPU)")
    
    cx.loglog(fineBins, henyeyGreenstein(numpy.cos(numpy.array(fineBins)*math.pi/180.), mean_cosTheta_inSim), linewidth=2, color='m', label=r"\textbf{Henyey-Greenstein (HG)} $\left<\cos\theta\right>=%5.3f$" % (mean_cosTheta_inSim))
    
    cx.loglog(fineBins, simplifiedLiu(numpy.cos(numpy.array(fineBins)*math.pi/180.), mean_cosTheta_inSim), linewidth=2, color='g', label=r"\textbf{simplified Liu (sL)} $\left<\cos\theta\right>=%5.3f$" % (mean_cosTheta_inSim))
    
    cx.loglog(fineBins, factor_SPICE_Mie * simplifiedLiu(numpy.cos(fineBins*math.pi/180.), mean_cosTheta_inSim)  +  (1.-factor_SPICE_Mie)*henyeyGreenstein(numpy.cos(fineBins*math.pi/180.), mean_cosTheta_inSim), linewidth=2, color='b', label=r"$%4.2f \times \textbf{HG} + %4.2f \times \textbf{sL}$ $\left<\cos\theta\right>=%5.3f$" % (1.-factor_SPICE_Mie, factor_SPICE_Mie, mean_cosTheta_inSim))

    cx.loglog(fineBins, factor_SPICE_Mie * simplifiedLiu(numpy.cos(fineBins*math.pi/180.), mean_cosTheta)  +  (1.-factor_SPICE_Mie)*henyeyGreenstein(numpy.cos(fineBins*math.pi/180.), mean_cosTheta), linewidth=1, linestyle='--', color='b', label=r"$%4.2f \times \textbf{HG} + %4.2f \times \textbf{sL}$ $\left<\cos\theta\right>=%5.3f$" % (1.-factor_SPICE_Mie, factor_SPICE_Mie, mean_cosTheta))
    
    cx.set_xlabel("scattering angle $\\theta [^\\circ]$")
    cx.set_ylabel("$\\beta(\\theta)$")
    cx.grid(True)
    #cx.legend(loc="upper right")
    cx.set_xlim(0.1,180.)
    cx.set_ylim(1e-3,1e3)


if True:
    fineBins = numpy.linspace(-1., 1., num=1000)
    
    dx.semilogy(hist_cosangles["cos"]["bins"], hist_cosangles["cos"]["num"], linewidth=2, color='r', label="MC generated (OpenCL)")
    dx.semilogy(hist_cosangles_host["cos"]["bins"], hist_cosangles_host["cos"]["num"], linewidth=2, color='y', label="MC generated (C++/CPU)")
    
    dx.semilogy(fineBins, henyeyGreenstein(fineBins, mean_cosTheta_inSim), linewidth=2, color='m', label=r"\textbf{Henyey-Greenstein (HG)} $\left<\cos\theta\right>=%5.3f$" % (mean_cosTheta_inSim))
    
    dx.semilogy(fineBins, simplifiedLiu(fineBins, mean_cosTheta_inSim), linewidth=2, color='g', label=r"\textbf{simplified Liu (sL)} $\left<\cos\theta\right>=%5.3f$" % (mean_cosTheta_inSim))
    
    dx.semilogy(fineBins, factor_SPICE_Mie * simplifiedLiu(fineBins, mean_cosTheta_inSim)  +  (1.-factor_SPICE_Mie)*henyeyGreenstein(fineBins, mean_cosTheta_inSim), linewidth=2, color='b', label=r"$%4.2f \times \textbf{HG} + %4.2f \times \textbf{sL}$ $\left<\cos\theta\right>=%5.3f$" % (1.-factor_SPICE_Mie, factor_SPICE_Mie, mean_cosTheta_inSim))
    
    dx.semilogy(fineBins, factor_SPICE_Mie * simplifiedLiu(fineBins, mean_cosTheta)  +  (1.-factor_SPICE_Mie)*henyeyGreenstein(fineBins, mean_cosTheta), linewidth=1, linestyle='--', color='b', label=r"$%4.2f \times \textbf{HG} + %4.2f \times \textbf{sL}$ $\left<\cos\theta\right>=%5.3f$" % (1.-factor_SPICE_Mie, factor_SPICE_Mie, mean_cosTheta))

    dx.set_xlabel(r"scattering angle $\cos\theta [^\circ]$")
    dx.set_ylabel("$\\beta(\\theta)$")
    dx.grid(True)
    dx.legend(loc="upper left")
    dx.set_ylim(1e-3,1e3)



pylab.savefig("scattering_angle_distribution_IceCube.pdf", transparent=False)
