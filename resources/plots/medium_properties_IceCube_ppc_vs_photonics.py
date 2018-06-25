#!/usr/bin/env python

from __future__ import print_function

import math
import numpy

from os.path import expandvars

import matplotlib
matplotlib.use("PDF")

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
import scipy.misc

nm=1.
m=1.


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


def addAnnotationToPlot(plot, text, loc=1, size=6.):
    from mpl_toolkits.axes_grid. anchored_artists import AnchoredText
    at = AnchoredText(text,
                      prop=dict(size=size), frameon=True,
                      loc=loc,
                      )
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    plot.add_artist(at)


#modelName="AHA"
#meanCos=0.94

modelName="SPICE-Mie"


if modelName=="AHA":
    meanCosString = "%4.2f" % meanCos
    meanCosString2 = "%03.0f" % (meanCos*100)
    mediumPropsPPC = clsim.MakeIceCubeMediumProperties(iceDataDirectory=expandvars("$I3_BUILD/clsim/resources/ice/ppc_aha_%s/" % meanCosString))

    photonicsFilename = "Ice_table.aha.i3coords.cos%s.17may2007.txt" % meanCosString2
    photonicsFilenameLatex = "Ice\_table.aha.i3coords.cos%s.17may2007.txt" % meanCosString2
    mediumPropsPhotonics = clsim.MakeIceCubeMediumPropertiesPhotonics(tableFile=expandvars("$I3_BUILD/clsim/resources/ice/photonics_aha/" + photonicsFilename))
elif modelName=="SPICE-Mie":
    meanCos=0.9 # override the cosine theta
    mediumPropsPPC = clsim.MakeIceCubeMediumProperties(iceDataDirectory=expandvars("$I3_BUILD/ice-models/resources/models/spice_mie/"))

    photonicsFilename = "Ice_table.mie.i3coords.cos090.08Apr2011.txt"
    photonicsFilenameLatex = "Ice\_table.mie.i3coords.cos090.08Apr2011.txt"
    mediumPropsPhotonics = clsim.MakeIceCubeMediumPropertiesPhotonics(tableFile=expandvars("$I3_BUILD/clsim/resources/ice/photonics_spice_mie/" + photonicsFilename))
else:
    raise RuntimeError("unknown model:", modelName)

print("numer of layers (PPC      ) is {}, starting at z={}m, height={}m".format(mediumPropsPPC.LayersNum, mediumPropsPPC.LayersZStart/I3Units.m, mediumPropsPPC.LayersHeight/I3Units.m))
print("numer of layers (photonics) is {}, starting at z={}m, height={}m".format(mediumPropsPhotonics.LayersNum, mediumPropsPhotonics.LayersZStart/I3Units.m, mediumPropsPhotonics.LayersHeight/I3Units.m))

currentLayer = 0


####

print("a")

fig = pylab.figure(3)
fig.subplots_adjust(left=0.05, bottom=0.05, top=0.95, right=0.98)

ax =  fig.add_subplot(4, 2, 1)
ax2 = fig.add_subplot(4, 2, 3)
bx =  fig.add_subplot(4, 2, 2)
bx2 = fig.add_subplot(4, 2, 4)
cx =  fig.add_subplot(4, 2, 6)
cx2 = fig.add_subplot(4, 2, 8)

cx2_prime = fig.add_subplot(2,2,3)

wlens=numpy.linspace(300.,600.,num=100)



nphase_reference_PPC = applyOpenCLWlenDependentFunction(wlens, mediumPropsPPC.GetPhaseRefractiveIndex(0), useReferenceFunction=True)
if mediumPropsPPC.GetGroupRefractiveIndexOverride(0) is None:
    ngroup_reference_PPC = None
else:
    ngroup_reference_PPC = applyOpenCLWlenDependentFunction(wlens, mediumPropsPPC.GetGroupRefractiveIndexOverride(0), useReferenceFunction=True)

nphase_reference_photonics = applyOpenCLWlenDependentFunction(wlens, mediumPropsPhotonics.GetPhaseRefractiveIndex(0), useReferenceFunction=True)
if mediumPropsPhotonics.GetGroupRefractiveIndexOverride(0) is None:
    ngroup_reference_photonics = None
else:
    ngroup_reference_photonics = applyOpenCLWlenDependentFunction(wlens, mediumPropsPhotonics.GetGroupRefractiveIndexOverride(0), useReferenceFunction=True)


print("c")

ax.plot(wlens, nphase_reference_PPC,       linewidth=3., color='k',   linestyle='solid', label=r"$n_\mathrm{p}$ (PPC)")
ax.plot(wlens, nphase_reference_photonics, linewidth=3., color='0.5', linestyle='solid', label=r"$n_\mathrm{p}$ (photonics)")
if ngroup_reference_PPC is not None:
    l, = ax.plot(wlens, ngroup_reference_PPC,       linewidth=3., color='k',   label=r"$n_\mathrm{g}$ (PPC)")
    l.set_dashes([5,5])
    
if ngroup_reference_photonics is not None:
    l, = ax.plot(wlens, ngroup_reference_photonics, linewidth=3., color='0.5', label=r"$n_\mathrm{g}$ (photonics)")
    l.set_dashes([5,5])


ax2.plot(wlens, nphase_reference_PPC/nphase_reference_photonics, linewidth=3., color='k',   linestyle='solid', label=r"$n_\mathrm{p}$ (PPC/photonics)")
l, = ax2.plot(wlens, ngroup_reference_PPC/ngroup_reference_photonics, linewidth=3., color='k',   linestyle='solid', label=r"$n_\mathrm{g}$ (PPC/photonics)")
l.set_dashes([5,5])

ax2.plot([wlens[0], wlens[-1]], [1.,1.], linewidth=1., color='k', linestyle='-')


ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(19))
ax.set_xlim(305.,595.)
ax.legend()
ax.grid(True)
ax.set_xlabel("wavelength $\\lambda [\\mathrm{nm}]$")
ax.set_ylabel("refractive index $n$")

ax2.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(19))
ax2.set_xlim(305.,595.)
ax2.set_xlabel("wavelength $\\lambda [\\mathrm{nm}]$")
ax2.set_ylabel("ratio PPC/photonics")
ax2.legend()
ax2.grid(True)
ax2.set_ylim(0.995,1.005)


for layerNum in range(1, mediumPropsPhotonics.GetLayersNum()): # skip bottom layer (number 0)
    zStart_PPC = mediumPropsPPC.GetLayersZStart() + float(layerNum) * mediumPropsPPC.GetLayersHeight()
    zStart_photonics = mediumPropsPhotonics.GetLayersZStart() + float(layerNum) * mediumPropsPhotonics.GetLayersHeight()
    
    print("layer {} starts at z={}m (PPC) / z={}m (photonics)".format(layerNum, zStart_PPC/I3Units.m, zStart_photonics/I3Units.m))
    
    absCoeff_PPC = 1./applyOpenCLMediumPropertyFunction(wlens, layerNum, mediumPropsPPC, mode="absorptionLength")
    bx.plot(wlens, absCoeff_PPC, linewidth=1., linestyle='-', color='b', alpha=0.1, label=r"(PPC) layer %u" % (layerNum))
    #print absCoeff_PPC

    absCoeff_photonics = 1./applyOpenCLMediumPropertyFunction(wlens, layerNum, mediumPropsPhotonics, mode="absorptionLength")
    bx.plot(wlens, absCoeff_photonics, linewidth=1., linestyle='-', color='g', alpha=0.1, label=r"(photonics) layer %u" % (layerNum))
    #print absCoeff_photonics

    bx2.plot(wlens, absCoeff_PPC/absCoeff_photonics, linewidth=1., linestyle='-', color='k', alpha=0.1, label=r"PPC / photonics layer %u" % (layerNum))
 
    
bx2.plot([wlens[0], wlens[-1]], [1.,1.], linewidth=1., color='k', linestyle='-')


bx.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(19))
bx.set_xlim(305.,595.)
#bx.set_ylim(0.,100.)
#bx.legend(loc="upper left")
bx.grid(True)
bx.set_xlabel("wavelength $\\lambda [\\mathrm{nm}]$")
bx.set_ylabel("absorption coeff $[\\mathrm{m}^{-1}]$")
#bx.set_ylim(0., 0.2)

bx2.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(19))
bx2.set_xlim(305.,595.)
bx2.set_xlabel("wavelength $\\lambda [\\mathrm{nm}]$")
bx2.set_ylabel("ratio PPC/photonics")
bx2.grid(True)
#bx2.legend(loc="upper left")
bx2.set_ylim(0.8, 1.2)


detectorCenterDepth = 1948.07*I3Units.m
#stopLayer = 146
stopLayer = 171

#for layerNum in range(1,stopLayer):
for layerNum in range(1, mediumPropsPhotonics.GetLayersNum()): # skip bottom layer (number 0)
    zStart_PPC = mediumPropsPPC.GetLayersZStart() + float(layerNum) * mediumPropsPPC.GetLayersHeight()
    zStart_photonics = mediumPropsPhotonics.GetLayersZStart() + float(layerNum) * mediumPropsPhotonics.GetLayersHeight()
    
    print("layer {} starts at z={}m (PPC) / z={}m (photonics)".format(layerNum, zStart_PPC/I3Units.m, zStart_photonics/I3Units.m))
    
    scatCoeff_PPC = 1./applyOpenCLMediumPropertyFunction(wlens, layerNum, mediumPropsPPC, mode="scatteringLength")
    cx.plot(wlens, scatCoeff_PPC, linewidth=1., linestyle='-', color='b', alpha=0.1, label=r"(PPC) layer %u" % (layerNum))
    #print scatCoeff_PPC

    scatCoeff_photonics = 1./applyOpenCLMediumPropertyFunction(wlens, layerNum, mediumPropsPhotonics, mode="scatteringLength")
    cx.plot(wlens, scatCoeff_photonics, linewidth=1., linestyle='-', color='g', alpha=0.1, label=r"(photonics) layer %u" % (layerNum))
    #print scatCoeff_photonics

    cx2.plot(wlens, scatCoeff_PPC/scatCoeff_photonics, linewidth=1., linestyle='-', color='k', alpha=0.2, label=r"PPC / photonics layer %u" % (layerNum))
    cx2_prime.plot(wlens, scatCoeff_PPC/scatCoeff_photonics, linewidth=1., linestyle='-', color='k', alpha=0.2, label=r"PPC / photonics layer %u" % (layerNum))
    

cx2.plot([wlens[0], wlens[-1]], [1.,1.], linewidth=1., color='k', linestyle='-')

cx.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(19))
cx.set_xlim(305.,595.)
#cx.set_ylim(0.,100.)
#cx.legend(loc="upper left")
cx.grid(True)
cx.set_xlabel("wavelength $\\lambda [\\mathrm{nm}]$")
cx.set_ylabel("geom. scattering coeff $[\\mathrm{m}^{-1}]$")
cx.set_ylim(0., 4.0)

cx2.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(19))
cx2.set_xlim(305.,595.)
cx2.set_xlabel("wavelength $\\lambda [\\mathrm{nm}]$")
cx2.set_ylabel("ratio PPC/photonics")
cx2.grid(True)
#cx2.legend(loc="upper left")
cx2.set_ylim(0.95, 1.05)


cx2_prime.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(19))
cx2_prime.set_xlim(305.,595.)
cx2_prime.set_xlabel("wavelength $\\lambda [\\mathrm{nm}]$")
cx2_prime.set_ylabel("ratio PPC/photonics")
cx2_prime.grid(True)
#cx2_prime.legend(loc="upper left")
cx2_prime.set_ylim(0.7, 1.3)



depthStart = detectorCenterDepth-mediumPropsPPC.GetLayersZStart() - mediumPropsPPC.GetLayersHeight() # one layer up
depthStop = detectorCenterDepth-mediumPropsPPC.GetLayersZStart() - mediumPropsPPC.GetLayersHeight() * float(stopLayer+1) # stopLayers x up
addAnnotationToPlot(cx, r"layers 1-{} (depth from {}m to {}m)".format(stopLayer, depthStart/I3Units.m, depthStop/I3Units.m))
addAnnotationToPlot(cx2, r"layers 1-{} (depth from {}m to {}m)".format(stopLayer, depthStart/I3Units.m, depthStop/I3Units.m))
addAnnotationToPlot(cx2_prime, r"layers 1-{} (depth from {}m to {}m)".format(stopLayer, depthStart/I3Units.m, depthStop/I3Units.m))


#for layerNum in range(0, mediumProps.LayersNum, 10):
#    scatCoeff = 1./applyOpenCLMediumPropertyFunction(wlens, layerNum, mediumProps, mode="scatteringLength")
#    cx.plot(wlens, 1./scatCoeff, linewidth=1., linestyle='-', color='b', alpha=0.1) #, label=r"layer %u" % (layerNum))
#
#scatCoeff = 1./applyOpenCLWlenDependentFunction(wlens, mediumProps.GetScatteringLength(exampleLayerNum), useReferenceFunction=True)
#cx.plot(wlens, 1./scatCoeff, linewidth=3., linestyle='-', color='k', label=r"C++")
#
#scatCoeff = 1./applyOpenCLMediumPropertyFunction(wlens, exampleLayerNum, mediumProps, mode="scatteringLength")
#cx.plot(wlens, 1./scatCoeff, linewidth=1., linestyle='-', color='r', label=r"OpenCL")
#
#
#cx.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(20))
#cx.set_xlim(300.,600.)
#cx.grid(True)
#cx.legend(loc='upper right')
#cx.set_xlabel(r"wavelength $\lambda [\mathrm{nm}]$")
#cx.set_ylabel(r"scattering length $\lambda_\mathrm{scat;geom}$ $[\mathrm{m}]$")

pylab.suptitle(r"model: %s, $\cos \theta = %4.2f$, photonics file: %s" % (modelName, meanCos, photonicsFilenameLatex))

pylab.savefig("medium_properties_IceCube_%s_ppc_vs_photonics_cos%4.2f.pdf" % (modelName, meanCos), transparent=False)

