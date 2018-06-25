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

#rng = phys_services.I3SPRNGRandomService(seed=3244, nstreams=2, streamnum=0)

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



def applyOpenCLWlenDependentFunction(xValues, functionOpenCL, useReferenceFunction=False):
    #print "         number of values:", len(xValues)
    
    tester = clsim.I3CLSimFunctionTester(device=openCLDevice,
                                                   workgroupSize=workgroupSize,
                                                   workItemsPerIteration=workItemsPerIteration,
                                                   wlenDependentValue=functionOpenCL)
    
    #print "maxWorkgroupSizeForKernel:", tester.maxWorkgroupSize
    
    # the function currently only accepts I3VectorFloat as its input type
    vector = dataclasses.I3VectorFloat(numpy.array(xValues))
    
    if useReferenceFunction:
        yValues = numpy.array(tester.EvaluateReferenceFunction(vector))
    else:
        yValues = numpy.array(tester.EvaluateFunction(vector))
    
    return yValues


domAngularAcceptance_holeIce = clsim.GetIceCubeDOMAngularSensitivity(holeIce=expandvars("$I3_BUILD/ice-models/resources/models/angsens/as.h2-50cm"))
domAngularAcceptance = clsim.GetIceCubeDOMAngularSensitivity(holeIce=expandvars("$I3_BUILD/ice-models/resources/models/angsens/as.nominal"))

####


fig = pylab.figure(3)
fig.subplots_adjust(left=0.09, bottom=0.05, top=0.95, right=0.98)

ax = fig.add_subplot(3, 1, 1)
bx = fig.add_subplot(3, 1, 2)
cx = fig.add_subplot(3, 1, 3)


cosines=numpy.linspace(-1.,1.,num=10000)


acceptance_reference = applyOpenCLWlenDependentFunction(cosines, domAngularAcceptance, useReferenceFunction=True)
acceptance_OpenCL = applyOpenCLWlenDependentFunction(cosines, domAngularAcceptance, useReferenceFunction=False)

acceptance_reference_holeIce = applyOpenCLWlenDependentFunction(cosines, domAngularAcceptance_holeIce, useReferenceFunction=True)
acceptance_OpenCL_holeIce = applyOpenCLWlenDependentFunction(cosines, domAngularAcceptance_holeIce, useReferenceFunction=False)


ax.semilogy(cosines, acceptance_reference_holeIce, linewidth=3., color='k', linestyle='--', label=r"pure acceptance (hole ice) (C++)")
ax.semilogy(cosines, acceptance_OpenCL_holeIce, linewidth=1., color=(0.3,0.3,1.0), linestyle='--', label=r"pure acceptance (hole ice) (OpenCL)")

ax.semilogy(cosines, acceptance_reference, linewidth=3., color='k', linestyle='-', label=r"pure acceptance (nominal) (C++)")
ax.semilogy(cosines, acceptance_OpenCL, linewidth=1., color='r', linestyle='-', label=r"pure acceptance (nominal) (OpenCL)")





ax.set_xlim(-1.,1.)
ax.legend(loc='lower right')
ax.grid(True)
ax.set_xlabel("wavelength $\\lambda [\\mathrm{nm}]$")
ax.set_ylabel("DOM acceptance")



pylab.savefig("dom_angular_sensitivity.pdf", transparent=True)



