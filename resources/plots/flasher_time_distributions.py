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
import scipy.interpolate
import scipy.integrate

from os.path import expandvars

from icecube import icetray, dataclasses, clsim, phys_services
from I3Tray import I3Units


scanned_FB_WIDTH15  = numpy.loadtxt(expandvars("$I3_BUILD/clsim/resources/flasher_data/optical_pulse_shape_FB_WIDTH15.txt"),  unpack=True)
scanned_FB_WIDTH20  = numpy.loadtxt(expandvars("$I3_BUILD/clsim/resources/flasher_data/optical_pulse_shape_FB_WIDTH20.txt"),  unpack=True)
scanned_FB_WIDTH124 = numpy.loadtxt(expandvars("$I3_BUILD/clsim/resources/flasher_data/optical_pulse_shape_FB_WIDTH124.txt"), unpack=True)

scanned_FB_WIDTH15[1] = scanned_FB_WIDTH15[1] / numpy.max(scanned_FB_WIDTH15[1])    # this one needs some re-scaling

scanned_FB_WIDTH20[1] = scanned_FB_WIDTH20[1] / numpy.max(scanned_FB_WIDTH20[1])    # this one also needs re-scaling
scanned_FB_WIDTH20[0] = scanned_FB_WIDTH20[0] - 22.88473                            # and has an offset, too



rng = phys_services.I3SPRNGRandomService(seed=3244, nstreams=2, streamnum=0)


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


fig = pylab.figure(3)
fig.subplots_adjust(left=0.09, bottom=0.05, top=0.95, right=0.98)

ax = fig.add_subplot(3, 2, 1)
bx = fig.add_subplot(3, 2, 2)
cx = fig.add_subplot(3, 2, 3)
dx = fig.add_subplot(3, 2, 4)
ex = fig.add_subplot(3, 2, 5)
fx = fig.add_subplot(3, 2, 6)


plotProfileAndMC(ax, FB_WIDTH=15.,  color1='k', color2='r', label="width: 15")
plotProfileAndMC(bx, FB_WIDTH=20.,  color1='k', color2='r', label="width: 20")
plotProfileAndMC(cx, FB_WIDTH=40.,  color1='k', color2='r', label="width: 40")
plotProfileAndMC(dx, FB_WIDTH=60.,  color1='k', color2='r', label="width: 60")
plotProfileAndMC(ex, FB_WIDTH=80.,  color1='k', color2='r', label="width: 80")
plotProfileAndMC(fx, FB_WIDTH=124., color1='k', color2='r', label="width: 124")

ax.plot(scanned_FB_WIDTH15[0],  scanned_FB_WIDTH15[1],  linestyle='--', color='b', label="scanned from wiki")
bx.plot(scanned_FB_WIDTH20[0],  scanned_FB_WIDTH20[1],  linestyle='--', color='b', label="scanned from wiki")
fx.plot(scanned_FB_WIDTH124[0], scanned_FB_WIDTH124[1], linestyle='--', color='b', label="scanned from wiki")


for x in [ax, bx, cx, dx, ex, fx]:
    x.set_xlim(-100.,120.)
    x.set_ylim(0.,1.1)
    #x.set_ylim(1e-4,1.1)
    x.legend(loc='upper left')
    x.grid(True)
    x.set_xlabel("time delay $[\\mathrm{ns}]$")
    x.set_ylabel("a.u.")

pylab.savefig("flasher_time_distributions.pdf", transparent=False)



