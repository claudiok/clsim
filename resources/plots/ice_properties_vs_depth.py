#!/usr/bin/env python

from __future__ import print_function

import math
import numpy

from os.path import expandvars

import matplotlib
matplotlib.use("PDF")

fig_size = [11.7,8.3] # din A4
params = {'backend': 'pdf',
        'axes.labelsize': 12,
        'text.fontsize': 12,
        'legend.fontsize': 8,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
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

wlens=numpy.linspace(300.,600.,num=100)



detectorCenterDepth = 1948.07*I3Units.m
plotAtWavelength = 400.



color_SPICE='g'
color_WHAM='b'





#mediumPropsSPICE = clsim.MakeIceCubeMediumProperties(iceDataDirectory=expandvars("$I3_BUILD/ice-models/resources/models/spice_mie/"))
mediumPropsSPICE = clsim.MakeIceCubeMediumPropertiesPhotonics(tableFile=expandvars("$I3_BUILD/clsim/resources/ice/photonics_spice_mie/Ice_table.mie.i3coords.cos090.08Apr2011.txt"))
name_SPICE="SPICE-Mie"
meanCos_SPICE = 0.9

#mediumPropsSPICE  = clsim.MakeIceCubeMediumPropertiesPhotonics(tableFile=expandvars("$I3_BUILD/clsim/resources/ice/photonics_aha/Ice_table.aha.i3coords.cos080.17may2007.txt"))
#name_SPICE="AHA080"
#meanCos_SPICE = 0.80





mediumPropsWHAM  = clsim.MakeIceCubeMediumPropertiesPhotonics(tableFile=expandvars("$I3_BUILD/clsim/resources/ice/photonics_wham/Ice_table.wham.i3coords.cos094.11jul2011.txt"))
name_WHAM="WHAM094"
meanCos_WHAM = 0.94

#mediumPropsWHAM  = clsim.MakeIceCubeMediumPropertiesPhotonics(tableFile=expandvars("$I3_BUILD/clsim/resources/ice/photonics_wham/Ice_table.wham.i3coords.cos090.11jul2011.txt"))
#name_WHAM="WHAM090"
#meanCos_WHAM = 0.90

#mediumPropsWHAM  = clsim.MakeIceCubeMediumPropertiesPhotonics(tableFile=expandvars("$I3_BUILD/clsim/resources/ice/photonics_wham/Ice_table.wham.i3coords.cos080.11jul2011.txt"))
#name_WHAM="WHAM080"
#meanCos_WHAM = 0.80


#mediumPropsWHAM  = clsim.MakeIceCubeMediumPropertiesPhotonics(tableFile=expandvars("$I3_BUILD/clsim/resources/ice/photonics_aha/Ice_table.aha.i3coords.cos080.17may2007.txt"))
#name_WHAM="AHA080"
#meanCos_WHAM = 0.80


#mediumPropsWHAM  = clsim.MakeIceCubeMediumPropertiesPhotonics(tableFile=expandvars("$I3_BUILD/clsim/resources/ice/photonics_aha/Ice_table.aha.i3coords.cos094.17may2007.txt"))
#name_WHAM="AHA094"
#meanCos_WHAM = 0.94



print("numer of layers ({}) is {}, starting at z={}m, height={}m".format(name_SPICE, mediumPropsSPICE.LayersNum, mediumPropsSPICE.LayersZStart/I3Units.m, mediumPropsSPICE.LayersHeight/I3Units.m))
print("numer of layers ({}) is {}, starting at z={}m, height={}m".format(name_WHAM,  mediumPropsWHAM.LayersNum, mediumPropsWHAM.LayersZStart/I3Units.m, mediumPropsWHAM.LayersHeight/I3Units.m))

currentLayer = 0


####

print("a")

fig = pylab.figure(3)
fig.subplots_adjust(left=0.06, bottom=0.05, top=0.96, right=0.98)

if False:
    ax =  fig.add_subplot(4, 2, 1)
    ax2 = fig.add_subplot(4, 2, 3)
    bx =  fig.add_subplot(4, 2, 2)
    bx2 = fig.add_subplot(4, 2, 4)
    cx =  fig.add_subplot(4, 2, 6)
    cx2 = fig.add_subplot(4, 2, 8)

if True:
    bx =  fig.add_subplot(2, 2, 1)
    bx2 = fig.add_subplot(2, 2, 3)
    cx =  fig.add_subplot(2, 2, 2)
    cx2 = fig.add_subplot(2, 2, 4)


if False:
    nphase_reference_SPICE = applyOpenCLWlenDependentFunction(wlens, mediumPropsSPICE.GetPhaseRefractiveIndex(0), useReferenceFunction=True)
    if mediumPropsSPICE.GetGroupRefractiveIndexOverride(0) is None:
        ngroup_reference_SPICE = None
    else:
        ngroup_reference_SPICE = applyOpenCLWlenDependentFunction(wlens, mediumPropsSPICE.GetGroupRefractiveIndexOverride(0), useReferenceFunction=True)

    nphase_reference_WHAM = applyOpenCLWlenDependentFunction(wlens, mediumPropsWHAM.GetPhaseRefractiveIndex(0), useReferenceFunction=True)
    if mediumPropsWHAM.GetGroupRefractiveIndexOverride(0) is None:
        ngroup_reference_WHAM = None
    else:
        ngroup_reference_WHAM = applyOpenCLWlenDependentFunction(wlens, mediumPropsWHAM.GetGroupRefractiveIndexOverride(0), useReferenceFunction=True)


    print("c")

    ax.plot(wlens, nphase_reference_SPICE, linewidth=2., color=color_SPICE, linestyle='solid', label=r"$n_\mathrm{p}$ " + "({})".format(name_SPICE))
    ax.plot(wlens, nphase_reference_WHAM,  linewidth=2., color=color_WHAM,  linestyle='solid', label=r"$n_\mathrm{p}$ " + "({})".format(name_WHAM))
    if ngroup_reference_SPICE is not None:
        l, = ax.plot(wlens, ngroup_reference_SPICE,       linewidth=3., color=color_SPICE,   label=r"$n_\mathrm{g}$ " + "({})".format(name_SPICE))
        l.set_dashes([5,5])
    
    if ngroup_reference_WHAM is not None:
        l, = ax.plot(wlens, ngroup_reference_WHAM, linewidth=2., color=color_WHAM, label=r"$n_\mathrm{g}$ " + "({})".format(name_WHAM))
        l.set_dashes([5,5])


    ax2.plot(wlens, nphase_reference_SPICE/nphase_reference_WHAM, linewidth=2., color='k',   linestyle='solid', label=r"$n_\mathrm{p}$ " + "({}/{})".format(name_SPICE,name_WHAM))
    l, = ax2.plot(wlens, ngroup_reference_SPICE/ngroup_reference_WHAM, linewidth=2., color='k',   linestyle='solid', label=r"$n_\mathrm{g}$ " + "({}/{})".format(name_SPICE,name_WHAM))
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
    ax2.set_ylabel("ratio {}/{}".format(name_SPICE,name_WHAM))
    ax2.legend()
    ax2.grid(True)
    ax2.set_ylim(0.995,1.005)


depths_SPICE = []
depths_WHAM = []
absorptionCoeffs_vals_SPICE = []
absorptionCoeffs_vals_WHAM = []
scatteringCoeffs_vals_SPICE = []
scatteringCoeffs_vals_WHAM = []

for layerNum in range(0, mediumPropsSPICE.GetLayersNum()):
    zmean_SPICE = mediumPropsSPICE.GetLayersZStart() + float(layerNum) * mediumPropsSPICE.GetLayersHeight() + mediumPropsSPICE.GetLayersHeight()/2.
    
    
    absCoeff_SPICE  = 1./applyOpenCLMediumPropertyFunction([plotAtWavelength-10.,plotAtWavelength,plotAtWavelength+10.], layerNum, mediumPropsSPICE, mode="absorptionLength")[1]
    scatCoeff_SPICE = 1./applyOpenCLMediumPropertyFunction([plotAtWavelength-10.,plotAtWavelength,plotAtWavelength+10.], layerNum, mediumPropsSPICE, mode="scatteringLength")[1]

    absorptionCoeffs_vals_SPICE.append(absCoeff_SPICE)
    scatteringCoeffs_vals_SPICE.append(scatCoeff_SPICE)
    
    depths_SPICE.append(detectorCenterDepth-zmean_SPICE)


for layerNum in range(0, mediumPropsWHAM.GetLayersNum()):
    zmean_WHAM  = mediumPropsWHAM.GetLayersZStart() + float(layerNum) * mediumPropsWHAM.GetLayersHeight() + mediumPropsWHAM.GetLayersHeight()/2.


    absCoeff_WHAM   = 1./applyOpenCLMediumPropertyFunction([plotAtWavelength-10.,plotAtWavelength,plotAtWavelength+10.], layerNum, mediumPropsWHAM, mode="absorptionLength")[1]
    scatCoeff_WHAM  = 1./applyOpenCLMediumPropertyFunction([plotAtWavelength-10.,plotAtWavelength,plotAtWavelength+10.], layerNum, mediumPropsWHAM, mode="scatteringLength")[1]

    absorptionCoeffs_vals_WHAM.append(absCoeff_WHAM)
    scatteringCoeffs_vals_WHAM.append(scatCoeff_WHAM)

    depths_WHAM.append(detectorCenterDepth-zmean_WHAM)

# re-order
depths_SPICE = depths_SPICE[::-1]
depths_WHAM = depths_WHAM[::-1]
absorptionCoeffs_vals_SPICE = absorptionCoeffs_vals_SPICE[::-1]
absorptionCoeffs_vals_WHAM = absorptionCoeffs_vals_WHAM[::-1]
scatteringCoeffs_vals_SPICE = scatteringCoeffs_vals_SPICE[::-1]
scatteringCoeffs_vals_WHAM = scatteringCoeffs_vals_WHAM[::-1]


depths_SPICE = numpy.array(depths_SPICE)
depths_WHAM = numpy.array(depths_WHAM)

absorptionCoeffs_vals_SPICE_interpolated = scipy.interpolate.interp1d(x=depths_SPICE, y=absorptionCoeffs_vals_SPICE, bounds_error=False)
absorptionCoeffs_vals_WHAM_interpolated  = scipy.interpolate.interp1d(x=depths_WHAM,  y=absorptionCoeffs_vals_WHAM,  bounds_error=False)
scatteringCoeffs_vals_SPICE_interpolated = scipy.interpolate.interp1d(x=depths_SPICE, y=scatteringCoeffs_vals_SPICE, bounds_error=False)
scatteringCoeffs_vals_WHAM_interpolated  = scipy.interpolate.interp1d(x=depths_WHAM,  y=scatteringCoeffs_vals_WHAM,  bounds_error=False)


manyDepths = numpy.linspace(min(depths_SPICE[0], depths_WHAM[0]), max(depths_SPICE[-1], depths_WHAM[-1]), 5000)


bx.semilogy(manyDepths/I3Units.m, absorptionCoeffs_vals_SPICE_interpolated(manyDepths), linewidth=2., linestyle='-', color=color_SPICE, alpha=1.0, label=r"{} @ {}nm".format(name_SPICE, plotAtWavelength))
bx.semilogy(manyDepths/I3Units.m, absorptionCoeffs_vals_WHAM_interpolated(manyDepths),  linewidth=2., linestyle='-', color=color_WHAM,  alpha=1.0, label=r"{} @ {}nm".format(name_WHAM, plotAtWavelength))

bx2.plot(manyDepths/I3Units.m, absorptionCoeffs_vals_SPICE_interpolated(manyDepths)/absorptionCoeffs_vals_WHAM_interpolated(manyDepths), linewidth=2., linestyle='-', color='k',   alpha=1.0)
bx2.plot([1000., 3000.], [1.,1.], linewidth=1., color='k', linestyle='--')


bx.grid(True)
bx.set_xlabel("depth $[\\mathrm{m}]$")
bx.set_ylabel("absorption coeff $[\\mathrm{m}^{-1}]$")
bx.set_ylim(0.003,0.1)
bx.legend(loc="upper center", ncol=2)
bx.set_xlim(1100.,2800.)

bx2.set_xlabel("depth $[\\mathrm{m}]$")
bx2.set_ylabel("ratio {}/{}".format(name_SPICE,name_WHAM))
bx2.grid(True)
#bx2.legend(loc="upper left")
bx2.set_ylim(0.4, 3.0)
bx2.set_xlim(1100.,2800.)




cx.semilogy(manyDepths/I3Units.m, scatteringCoeffs_vals_SPICE_interpolated(manyDepths)*(1.-meanCos_SPICE), linewidth=2., linestyle='-', color=color_SPICE, alpha=1.0, label=r"{} @ {}nm, $\left< \cos\theta \right>={}$".format(name_SPICE, plotAtWavelength, meanCos_SPICE))
cx.semilogy(manyDepths/I3Units.m, scatteringCoeffs_vals_WHAM_interpolated(manyDepths) *(1.-meanCos_WHAM),  linewidth=2., linestyle='-', color=color_WHAM,  alpha=1.0, label=r"{} @ {}nm, $\left< \cos\theta \right>={}$".format(name_WHAM, plotAtWavelength, meanCos_WHAM))

cx2.plot(manyDepths/I3Units.m, scatteringCoeffs_vals_SPICE_interpolated(manyDepths)*(1.-meanCos_SPICE)/(scatteringCoeffs_vals_WHAM_interpolated(manyDepths)*(1.-meanCos_WHAM)), linewidth=2., linestyle='-', color='k',   alpha=1.0)
cx2.plot([1000., 3000.], [1.,1.], linewidth=1., color='k', linestyle='--')


cx.grid(True)
cx.set_xlabel("depth $[\\mathrm{m}]$")
cx.set_ylabel("eff. scattering coeff $[\\mathrm{m}^{-1}]$")
cx.set_ylim(0.01,0.6)
cx.legend(loc="upper center", ncol=2)
cx.set_xlim(1100.,2800.)

cx2.set_xlabel("depth $[\\mathrm{m}]$")
cx2.set_ylabel("ratio {}/{}".format(name_SPICE,name_WHAM))
cx2.grid(True)
#cx2.legend(loc="upper left")
cx2.set_ylim(0.4, 3.0)
cx2.set_xlim(1100.,2800.)



#pylab.suptitle(r"model: %s, $\cos \theta = %4.2f$, photonics file: %s" % (modelName, meanCos, photonicsFilenameLatex))

outfileName = "ice_properties_vs_depth___{}_vs_{}.pdf".format(name_SPICE,name_WHAM)
pylab.savefig(outfileName, transparent=False)
print("wrote", outfileName)

