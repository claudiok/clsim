#!/usr/bin/env python

from __future__ import print_function

import os

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
import scipy.misc

nm=1.
m=1.

def absLenIceCube(wavelen, layer):
    if layer != 146: raise RuntimeError("Python implementation can only return values for layer 146.")
    
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



def getPhaseRefIndex(wavelength):
    x = wavelength/1000.# wavelength in micrometer
    return 1.55749 - 1.57988*x + 3.99993*x**2. - 4.68271*x**3. + 2.09354*x**4.

def _getDispersionPhase(wavelength):
    return scipy.misc.derivative(getPhaseRefIndex, wavelength)
getDispersionPhase = numpy.vectorize(_getDispersionPhase)

#def getDispersionPhase(wavelength):
#    x = wavelength/1000.# wavelength in micrometer
#    return (-1.57988 + 2.*3.99993*x - 3.*4.68271*x**2. + 4.*2.09354*x**3.)/1000.

def getGroupRefIndex_derivative(wavelength):
    n_inv = 1./getPhaseRefIndex(wavelength)
    y = getDispersionPhase(wavelength);
    return 1./((1.0 + y*wavelength*n_inv) * n_inv)

def getGroupRefIndex_parameterized(wavelength):
    np = getPhaseRefIndex(wavelength)
    x = wavelength/1000.# wavelength in micrometer
    return np * (1. + 0.227106 - 0.954648*x + 1.42568*x**2. - 0.711832*x**3.)
    
def Cherenkov_dN_dXdwlen(wlen, beta=1.):
    return (2.*math.pi/(137.*(wlen**2.)))*(1. - 1./((beta*getPhaseRefIndex(wlen))**2.))

#print Cherenkov_dN_dXdwlen(470.)

numberOfPhotonsPerNanometer, err = scipy.integrate.quadrature(Cherenkov_dN_dXdwlen, 265., 675.)
#print err
numberOfPhotonsPerMeter = numberOfPhotonsPerNanometer*1e9

print("photons per meter between [290..610]nm =", numberOfPhotonsPerMeter)



####



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



def applyOpenCLWlenDependentFunction(xValues, functionOpenCL, getDerivative=False, useReferenceFunction=False):
    print("         number of values:", len(xValues))
    
    tester = clsim.I3CLSimFunctionTester(device=openCLDevice,
                                                   workgroupSize=workgroupSize,
                                                   workItemsPerIteration=workItemsPerIteration,
                                                   wlenDependentValue=functionOpenCL)
    
    print("maxWorkgroupSizeForKernel:", tester.maxWorkgroupSize)
    
    # the function currently only accepts I3VectorFloat as its input type
    vector = dataclasses.I3VectorFloat(numpy.array(xValues)*I3Units.nanometer)

    if useReferenceFunction:
        if getDerivative:
            yValues = numpy.array(tester.EvaluateReferenceDerivative(vector))/(1./I3Units.nanometer)
        else:
            yValues = numpy.array(tester.EvaluateReferenceFunction(vector))
    else:
        if getDerivative:
            yValues = numpy.array(tester.EvaluateDerivative(vector))/(1./I3Units.nanometer)
        else:
            yValues = numpy.array(tester.EvaluateFunction(vector))

    return yValues

def applyOpenCLMediumPropertyFunction(xValues, layer, mediumProps, mode):
    tester = clsim.I3CLSimMediumPropertiesTester(device=openCLDevice,
                                                 workgroupSize=workgroupSize,
                                                 workItemsPerIteration=workItemsPerIteration,
                                                 mediumProperties=mediumProps,
                                                 randomService=None)
    
    vector = dataclasses.I3VectorFloat(numpy.array(xValues)*I3Units.nanometer)
    
    if mode=="phaseRefIndex":
        yValues = numpy.array(tester.EvaluatePhaseRefIndex(vector, layer))
    elif mode=="dispersion":
        yValues = numpy.array(tester.EvaluateDispersion(vector, layer))
    elif mode=="groupVelocity":
        yValues = numpy.array(tester.EvaluateGroupVelocity(vector, layer))
    elif mode=="absorptionLength":
        yValues = numpy.array(tester.EvaluateAbsorptionLength(vector, layer))
    elif mode=="scatteringLength":
        yValues = numpy.array(tester.EvaluateScatteringLength(vector, layer))
    else:
        raise RuntimeError("Mode \"%s\" is not valid." % mode)
    
    return yValues


#mediumProps = clsim.MakeIceCubeMediumPropertiesPhotonics(tableFile=os.path.expandvars("$I3_BUILD/clsim/resources/ice/photonics_wham/Ice_table.wham.i3coords.cos090.11jul2011.txt"))
mediumProps = clsim.MakeIceCubeMediumProperties()
print("numer of layers =", mediumProps.LayersNum)
currentLayer = mediumProps.LayersNum/2
print("  current layer =", currentLayer)
phaseRefIndex = mediumProps.GetPhaseRefractiveIndex(currentLayer)

####

print("a")

fig = pylab.figure(3)
fig.subplots_adjust(left=0.09, bottom=0.05, top=0.95, right=0.98)

ax = fig.add_subplot(3, 1, 1)
bx = fig.add_subplot(3, 1, 2)
cx = fig.add_subplot(3, 1, 3)

wlens=numpy.linspace(260.,690.,num=100)


ax.plot(wlens, getPhaseRefIndex(wlens), linewidth=6., color='g', linestyle='solid', label=r"$n_\mathrm{p}$ (python)")
l, = ax.plot(wlens, getGroupRefIndex_parameterized(wlens), linewidth=6., color='g', label=r"$n_\mathrm{g}$ (python)")
l.set_dashes([5,5])
l, = ax.plot(wlens, getGroupRefIndex_derivative(wlens), linewidth=3., color='0.5', label=r"$n_\mathrm{g}$ (calculated from $n_\mathrm{p}$)")
l.set_dashes([10,10])

print("b")


nphase_reference = applyOpenCLWlenDependentFunction(wlens, mediumProps.GetPhaseRefractiveIndex(currentLayer), useReferenceFunction=True)
if mediumProps.GetGroupRefractiveIndexOverride(currentLayer) is None:
    ngroup_reference = None
else:
    ngroup_reference = applyOpenCLWlenDependentFunction(wlens, mediumProps.GetGroupRefractiveIndexOverride(currentLayer), useReferenceFunction=True)

print("c")

ax.plot(wlens, nphase_reference, linewidth=3., color='k', linestyle='solid', label=r"$n_\mathrm{p}$ (C++)")
if ngroup_reference is not None:
    l, = ax.plot(wlens, ngroup_reference, linewidth=3., color='k', label=r"$n_\mathrm{g}$ (C++)")
    l.set_dashes([5,5])


print("d")


nphase = applyOpenCLMediumPropertyFunction(wlens, currentLayer, mediumProps, mode="phaseRefIndex")

# get the ngroup parameterization from the ice properties
ngroup = dataclasses.I3Constants.c/applyOpenCLMediumPropertyFunction(wlens, currentLayer, mediumProps, mode="groupVelocity")

# calculate from ngroup from nphase
dispersion = applyOpenCLMediumPropertyFunction(wlens, currentLayer, mediumProps, mode="dispersion")
ngroup_calc = nphase/(1. + dispersion*wlens*I3Units.nanometer/nphase)


ax.plot(wlens, nphase, linewidth=1., color='r', linestyle='solid', label=r"$n_\mathrm{p}$ (OpenCL)")
l, = ax.plot(wlens, ngroup, linewidth=1., color='r', linestyle='solid', label=r"$n_\mathrm{g}$ (OpenCL)")
l.set_dashes([5,5])
l, = ax.plot(wlens, ngroup_calc, linewidth=1., color=(0.5,0.0,0.0), linestyle='solid', label=r"$n_\mathrm{g}$ (from $n_\mathrm{p}$) (OpenCL)")
l.set_dashes([10,10])



ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(20))
ax.set_xlim(260.,690.)
ax.legend()
ax.grid(True)
ax.set_xlabel("wavelength $\\lambda [\\mathrm{nm}]$")
ax.set_ylabel("refractive index $n$")


for layerNum in range(0, mediumProps.LayersNum, 10):
    absCoeff = 1./applyOpenCLMediumPropertyFunction(wlens, layerNum, mediumProps, mode="absorptionLength")
    bx.plot(wlens, 1./absCoeff, linewidth=1., linestyle='-', color='b', alpha=0.1) #, label=r"layer %u" % (layerNum))


exampleLayerNum=114
#bx.plot(wlens, absLenIceCube(wlens, layer=exampleLayerNum), linewidth=6., linestyle='-', color='g', label=r"python")

absCoeff = 1./applyOpenCLWlenDependentFunction(wlens, mediumProps.GetAbsorptionLength(exampleLayerNum), useReferenceFunction=True)
bx.plot(wlens, 1./absCoeff, linewidth=3., linestyle='-', color='k', label=r"C++")

absCoeff = 1./applyOpenCLMediumPropertyFunction(wlens, exampleLayerNum, mediumProps, mode="absorptionLength")
bx.plot(wlens, 1./absCoeff, linewidth=1., linestyle='-', color='r', label=r"OpenCL")



bx.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(20))
bx.set_xlim(260.,690.)
#bx.set_ylim(0.,100.)
bx.legend(loc="upper left")
bx.grid(True)
bx.set_xlabel("wavelength $\\lambda [\\mathrm{nm}]$")
bx.set_ylabel("absorption length $[\\mathrm{m}]$")
#bx.set_ylim(0., 0.2)

for layerNum in range(0, mediumProps.LayersNum, 10):
    scatCoeff = 1./applyOpenCLMediumPropertyFunction(wlens, layerNum, mediumProps, mode="scatteringLength")
    cx.plot(wlens, 1./scatCoeff, linewidth=1., linestyle='-', color='b', alpha=0.1) #, label=r"layer %u" % (layerNum))

scatCoeff = 1./applyOpenCLWlenDependentFunction(wlens, mediumProps.GetScatteringLength(exampleLayerNum), useReferenceFunction=True)
cx.plot(wlens, 1./scatCoeff, linewidth=3., linestyle='-', color='k', label=r"C++")

scatCoeff = 1./applyOpenCLMediumPropertyFunction(wlens, exampleLayerNum, mediumProps, mode="scatteringLength")
cx.plot(wlens, 1./scatCoeff, linewidth=1., linestyle='-', color='r', label=r"OpenCL")


cx.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(20))
cx.set_xlim(260.,690.)
cx.grid(True)
cx.legend(loc='upper right')
cx.set_xlabel(r"wavelength $\lambda [\mathrm{nm}]$")
cx.set_ylabel(r"scattering length $\lambda_\mathrm{scat;geom}$ $[\mathrm{m}]$")

pylab.savefig("medium_properties_IceCube.pdf", transparent=False)

