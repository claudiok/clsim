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
        'legend.fontsize': 6,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'text.usetex': True,
        'figure.figsize': fig_size}
matplotlib.rcParams.update(params)
matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern']})

def addAnnotationToPlot(plot, text, loc=1, size=8., rotation=0.):
    from mpl_toolkits.axes_grid. anchored_artists import AnchoredText
    at = AnchoredText(text,
                      prop=dict(size=size, rotation=rotation),
                      frameon=True,
                      loc=loc, # 1=='upper right'
                      )
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    plot.add_artist(at)

import pylab
import scipy


from icecube import icetray, dataclasses, clsim, phys_services
from I3Tray import I3Units


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



rng = phys_services.I3SPRNGRandomService(seed=3244, nstreams=2, streamnum=0)

def genMCHistogramsOpenCL(distribution, hist_range, distribution_params=[], iterations=1000, numBins=1000):
    tester = clsim.I3CLSimRandomDistributionTester(device=openCLDevice,
                                                   workgroupSize=workgroupSize,
                                                   workItemsPerIteration=workItemsPerIteration,
                                                   randomService=rng,
                                                   randomDistribution=distribution,
                                                   runtimeParameters=distribution_params)
    
    values = tester.GenerateRandomNumbers(iterations)
    samples = len(values)
    print("generated")
    
    range_width=hist_range[1]-hist_range[0]
    
    num_orig, bins = scipy.histogram(values, range=hist_range, bins=numBins)
    print("hist1 complete")
    
    del values # not needed anymore
    print("deleted")
    
    num=[]
    for number in num_orig:
        num.append(float(number)/float(samples)/float(range_width/float(numBins)))
    num=numpy.array(num)
    
    bins = numpy.array(bins[:-1])+(bins[1]-bins[0])/2.
    
    return dict(num=num, bins=bins)


def genMCHistogramsHost(distribution, hist_range, distribution_params=[], iterations=100000, numBins=1000):
    print("generating (host)")

    values = []
    for i in range(iterations):
        values.append(distribution.SampleFromDistribution(rng, distribution_params))
    samples = len(values)
    print("generated (host)")
    
    range_width=hist_range[1]-hist_range[0]
    
    num_orig, bins = scipy.histogram(values, range=hist_range, bins=numBins)
    print("hist1 complete (host)")
    
    del values # not needed anymore
    print("deleted (host)")
    
    num=[]
    for number in num_orig:
        num.append(float(number)/float(samples)/float(range_width/float(numBins)))
    num=numpy.array(num)
    
    bins = numpy.array(bins[:-1])+(bins[1]-bins[0])/2.
    
    return dict(num=num, bins=bins)

def gauss(mu, sigma, x):
    return (1.0/(sigma*numpy.sqrt(2*numpy.pi)))*numpy.exp(-(x-mu)**2/(2.0*sigma**2))

def plotProfileAndMC(ax, FB_WIDTH, color1, color2, label=None, **kwargs):
    profile_dist = clsim.I3CLSimRandomValueIceCubeFlasherTimeProfile()
    
    xVals = numpy.linspace(-100.,200.,2000)
    area = scipy.integrate.trapz(profile_dist._the_pulse(xVals, FB_WIDTH), xVals) 

    hist = genMCHistogramsHost(profile_dist, hist_range=(0., 200.), distribution_params=[FB_WIDTH*0.5*I3Units.ns])
    ax.plot(hist["bins"], hist["num"]*area, color=color2, linewidth=2., label=label+" (MC)", **kwargs)

    ax.plot(xVals, profile_dist._the_pulse(xVals, FB_WIDTH), color=color1, label=label+" (func)", **kwargs)
    
    ax.plot(xVals, gauss(0., FB_WIDTH*0.5 / 2.3548, xVals)*area, linestyle='--', color='g', label=r"previous versions")
    
####



normal_dist = clsim.I3CLSimRandomValueNormalDistribution()
normal_dist_fixed_mean_0 = clsim.I3CLSimRandomValueFixParameter(normal_dist, 0, 0.)
normal_dist_fixed_mean_5 = clsim.I3CLSimRandomValueFixParameter(normal_dist, 0, 50.)
normal_dist_fixed_sigma_10 = clsim.I3CLSimRandomValueFixParameter(normal_dist, 1, 10.)




####

fig = pylab.figure(3)
fig.subplots_adjust(left=0.09, bottom=0.05, top=0.95, right=0.98)

ax = fig.add_subplot(3, 2, 1)
bx = fig.add_subplot(3, 2, 2)
cx = fig.add_subplot(3, 2, 3)
dx = fig.add_subplot(3, 2, 4)
ex = fig.add_subplot(3, 2, 5)
fx = fig.add_subplot(3, 2, 6)


xVals = numpy.linspace(-100., 100., 1000)

histHost = genMCHistogramsHost(normal_dist, hist_range=(-100., 100.), distribution_params=[0., 10.])
histOpenCL = genMCHistogramsOpenCL(normal_dist, hist_range=(-100., 100.), distribution_params=[0., 10.])
ax.plot(histHost["bins"], histHost["num"], color='r', linewidth=2., label=r"MC (Host)")
ax.plot(histOpenCL["bins"], histOpenCL["num"], color='g', linewidth=2., label=r"MC (OpenCL)")
ax.semilogy(xVals, gauss(0., 10., xVals), color='k', label=r"func")


histHost = genMCHistogramsHost(normal_dist, hist_range=(-100., 100.), distribution_params=[50., 10.])
histOpenCL = genMCHistogramsOpenCL(normal_dist, hist_range=(-100., 100.), distribution_params=[50., 10.])
bx.plot(histHost["bins"], histHost["num"], color='r', linewidth=2., label=r"MC (Host)")
bx.plot(histOpenCL["bins"], histOpenCL["num"], color='g', linewidth=2., label=r"MC (OpenCL)")
bx.semilogy(xVals, gauss(50., 10., xVals), color='k', label=r"func")




histHost = genMCHistogramsHost(normal_dist_fixed_mean_0, hist_range=(-100., 100.), distribution_params=[10.])
histOpenCL = genMCHistogramsOpenCL(normal_dist_fixed_mean_0, hist_range=(-100., 100.), distribution_params=[10.])
cx.plot(histHost["bins"], histHost["num"], color='r', linewidth=2., label=r"MC (Host)")
cx.plot(histOpenCL["bins"], histOpenCL["num"], color='g', linewidth=2., label=r"MC (OpenCL)")
cx.semilogy(xVals, gauss(0., 10., xVals), color='k', label=r"func")


histHost = genMCHistogramsHost(normal_dist_fixed_mean_5, hist_range=(-100., 100.), distribution_params=[10.])
histOpenCL = genMCHistogramsOpenCL(normal_dist_fixed_mean_5, hist_range=(-100., 100.), distribution_params=[10.])
dx.plot(histHost["bins"], histHost["num"], color='r', linewidth=2., label=r"MC (Host)")
dx.plot(histOpenCL["bins"], histOpenCL["num"], color='g', linewidth=2., label=r"MC (OpenCL)")
dx.semilogy(xVals, gauss(50., 10., xVals), color='k', label=r"func")




histHost = genMCHistogramsHost(normal_dist_fixed_sigma_10, hist_range=(-100., 100.), distribution_params=[0.])
histOpenCL = genMCHistogramsOpenCL(normal_dist_fixed_sigma_10, hist_range=(-100., 100.), distribution_params=[0.])
ex.plot(histHost["bins"], histHost["num"], color='r', linewidth=2., label=r"MC (Host)")
ex.plot(histOpenCL["bins"], histOpenCL["num"], color='g', linewidth=2., label=r"MC (OpenCL)")
ex.semilogy(xVals, gauss(0., 10., xVals), color='k', label=r"func")


histHost = genMCHistogramsHost(normal_dist_fixed_sigma_10, hist_range=(-100., 100.), distribution_params=[50.])
histOpenCL = genMCHistogramsOpenCL(normal_dist_fixed_sigma_10, hist_range=(-100., 100.), distribution_params=[50.])
fx.plot(histHost["bins"], histHost["num"], color='r', linewidth=2., label=r"MC (Host)")
fx.plot(histOpenCL["bins"], histOpenCL["num"], color='g', linewidth=2., label=r"MC (OpenCL)")
fx.semilogy(xVals, gauss(50., 10., xVals), color='k', label=r"func")




for x in [ax, bx, cx, dx, ex, fx]:
    x.set_xlim(-100.,100.)
    #x.set_ylim(0.,1.1)
    x.set_ylim(1e-8,1e-1)
    x.legend(loc='upper left')
    x.grid(True)
    x.set_xlabel("x")
    x.set_ylabel("y")

pylab.savefig("random_distributions_test.pdf", transparent=False)



