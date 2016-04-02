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

# generated on photocathode per (injected on 1m^2)
qe_dom2007a_with_area = [
    #0.0002027029,
    0.0000000000,
    0.0000064522,
    0.0000064522,
    0.0000064522,
    0.0000021980,
    0.0001339040,
    0.0005556810,
    0.0016953000,
    0.0035997000,
    0.0061340900,
    0.0074592700,
    0.0090579800,
    0.0099246700,
    0.0105769000,
    0.0110961000,
    0.0114214000,
    0.0114425000,
    0.0111527000,
    0.0108086000,
    0.0104458000,
    0.0099763100,
    0.0093102500,
    0.0087516600,
    0.0083225800,
    0.0079767200,
    0.0075625100,
    0.0066377000,
    0.0053335800,
    0.0043789400,
    0.0037583500,
    0.0033279800,
    0.0029212500,
    0.0025334900,
    0.0021115400,
    0.0017363300,
    0.0013552700,
    0.0010546600,
    0.0007201020,
    0.0004843820,
    0.0002911110,
    0.0001782310,
    0.0001144300,
    0.0000509155
]
qe_dom2007a_wlen = numpy.linspace(260.,260.+10.*42, 43)

#dom_radius = 0.178 # dom radius [m] [i3mcml]
dom_radius = 0.16510 # dom radius [m] [ppc]
dom_area = math.pi*dom_radius**2.

qe_dom2007a_unit_area = numpy.array(qe_dom2007a_with_area)/dom_area
qe_dom2007a = scipy.interpolate.interp1d(qe_dom2007a_wlen, qe_dom2007a_unit_area) #, kind='cubic')


nm=1.
m=1.
h_times_c = 1239.84172 # in eV*nm

def getPhaseRefIndex(wavelength):
    x = wavelength/1000.# wavelength in micrometer
    return 1.55749 - 1.57988*x + 3.99993*x**2. - 4.68271*x**3. + 2.09354*x**4.

def getGroupRefIndex(wavelength):
    x = wavelength/1000.# wavelength in micrometer
    np = getPhaseRefIndex(wavelength)
    return np * (1. + 0.227106 - 0.954648*x + 1.42568*x**2. - 0.711832*x**3.)
    
def Cherenkov_dN_dXdwlen(wlen, beta=1.):
    value = (2.*math.pi/(137.*(wlen**2.)))*(1. - 1./((beta*getPhaseRefIndex(wlen))**2.))
    return numpy.where(value>0.,value,0.)

def Cherenkov_dN_dXdE(energy, beta=1.):
    value = (2.*math.pi/(137.*h_times_c))*(1. - 1./((beta*getPhaseRefIndex(h_times_c/energy))**2.))
    return numpy.where(value>0.,value,0.)


from icecube import icetray, dataclasses, clsim, phys_services
from I3Tray import I3Units

rng = phys_services.I3SPRNGRandomService(seed=3244, nstreams=2, streamnum=0)

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
    #print "         number of values:", len(xValues)
    
    tester = clsim.I3CLSimFunctionTester(device=openCLDevice,
                                         workgroupSize=workgroupSize,
                                         workItemsPerIteration=workItemsPerIteration,
                                         wlenDependentValue=functionOpenCL)
    
    #print "maxWorkgroupSizeForKernel:", tester.maxWorkgroupSize
    
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

def genMCHistogramsOpenCL(distribution, hist_range, iterations=1000, numBins=1000):
    tester = clsim.I3CLSimRandomDistributionTester(device=openCLDevice,
                                                   workgroupSize=workgroupSize,
                                                   workItemsPerIteration=workItemsPerIteration,
                                                   randomService=rng,
                                                   randomDistribution=distribution)
    
    values = tester.GenerateRandomNumbers(iterations)
    samples = len(values)
    print("generated")
    
    values = numpy.array(values)/I3Units.nanometer # convert to numpy array and convert units
    print("converted")
    
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

def genMCHistogramsHost(distribution, hist_range, iterations=10000000, numBins=1000):
    print("generating (host)")

    values = []
    for i in range(iterations):
        values.append(distribution.SampleFromDistribution(rng, []))
    samples = len(values)
    print("generated (host)")
    
    values = numpy.array(values)/I3Units.nanometer # convert to numpy array and convert units
    print("converted (host)")
    
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

beta=1.
mediumProps = clsim.MakeIceCubeMediumProperties()
domAcceptance = clsim.GetIceCubeDOMAcceptance()
flatAcceptance = clsim.I3CLSimFunctionConstant(1.)
phaseRefIndex = mediumProps.GetPhaseRefractiveIndex(0)

wlen_range = (mediumProps.GetMinWavelength()/I3Units.nanometer, mediumProps.GetMaxWavelength()/I3Units.nanometer)

genWavelength = clsim.makeCherenkovWavelengthGenerator(domAcceptance, False, mediumProps)
histGenWavelength = genMCHistogramsOpenCL(genWavelength, hist_range=wlen_range)
histGenWavelengthHost = genMCHistogramsHost(genWavelength, hist_range=wlen_range)
numberOfPhotonsPerMeter = clsim.NumberOfPhotonsPerMeter(phaseRefIndex, domAcceptance, wlen_range[0]*I3Units.nanometer, wlen_range[1]*I3Units.nanometer)

genWavelengthFlat = clsim.makeCherenkovWavelengthGenerator(flatAcceptance, False, mediumProps)
histGenWavelengthFlat = genMCHistogramsOpenCL(genWavelengthFlat, hist_range=wlen_range)
histGenWavelengthFlatHost = genMCHistogramsHost(genWavelengthFlat, hist_range=wlen_range)
numberOfPhotonsPerMeterFlat = clsim.NumberOfPhotonsPerMeter(phaseRefIndex, flatAcceptance, wlen_range[0]*I3Units.nanometer, wlen_range[1]*I3Units.nanometer)

#genWavelengthFlatNoDispersion = clsim.makeCherenkovWavelengthGenerator(flatAcceptance, True, mediumProps)
#histGenWavelengthFlatNoDispersion = genMCHistogramsOpenCL(genWavelengthFlatNoDispersion, range=wlen_range)

####


fig = pylab.figure(3)
fig.subplots_adjust(left=0.09, bottom=0.05, top=0.95, right=0.98)

ax = fig.add_subplot(3, 1, 1)
bx = fig.add_subplot(3, 1, 2)
cx = fig.add_subplot(3, 1, 3)


wlens=numpy.linspace(wlen_range[0],wlen_range[1],num=10000)
#energies=numpy.linspace(h_times_c/680.,h_times_c/260.,num=10000)


acceptance_reference = applyOpenCLWlenDependentFunction(wlens, domAcceptance, useReferenceFunction=True)
acceptance_OpenCL = applyOpenCLWlenDependentFunction(wlens, domAcceptance, useReferenceFunction=False)

ax.plot(wlens, qe_dom2007a(wlens), linewidth=6., color='g', linestyle='solid', label=r"pure acceptance (python)")
ax.plot(wlens, acceptance_reference, linewidth=3., color='k', linestyle='solid', label=r"pure acceptance (C++)")
ax.plot(wlens, acceptance_OpenCL, linewidth=1., color='r', linestyle='solid', label=r"pure acceptance (OpenCL)")



addAnnotationToPlot(bx, loc=2, text=r"$\int_{-\infty}^{\infty} \left( \frac{\mathrm{d}N_\mathrm{phot}}{\mathrm{d}l \mathrm{d}\lambda} \times \epsilon_\mathrm{DOM} \right) \mathrm{d}\lambda = %.1f \, \mathrm{m}^{-1}$" % (numberOfPhotonsPerMeter))
bx.plot(histGenWavelengthHost["bins"], histGenWavelengthHost["num"]*numberOfPhotonsPerMeter, linewidth=2, color='g', label="MC generated (C++/CPU)")
bx.plot(histGenWavelength["bins"], histGenWavelength["num"]*numberOfPhotonsPerMeter, linewidth=2, color='r', label="MC generated (OpenCL)")
bx.plot(wlens, qe_dom2007a(wlens)*Cherenkov_dN_dXdwlen(wlens, beta)*1.*1e9, linewidth=2., color='k', linestyle='solid', label=r"acceptance $\times$ cherenkov spectrum (python)")


addAnnotationToPlot(cx, loc=2, text=r"$\int_{-\infty}^{\infty} \frac{\mathrm{d}N_\mathrm{phot}}{\mathrm{d}l \mathrm{d}\lambda} \mathrm{d}\lambda = %.1f \, \mathrm{m}^{-1}$" % (numberOfPhotonsPerMeterFlat))

#cx.plot(histGenWavelengthFlatNoDispersion["bins"], histGenWavelengthFlatNoDispersion["num"]*numberOfPhotonsPerMeterFlat, linewidth=2, color='b', label="MC generated (no dispersion)")
cx.plot(histGenWavelengthFlatHost["bins"], histGenWavelengthFlatHost["num"]*numberOfPhotonsPerMeterFlat, linewidth=2, color='g', label="MC generated (C++/CPU)")
cx.plot(histGenWavelengthFlat["bins"], histGenWavelengthFlat["num"]*numberOfPhotonsPerMeterFlat, linewidth=2, color='r', label="MC generated (OpenCL)")
cx.plot(wlens, Cherenkov_dN_dXdwlen(wlens, beta)*1.*1e9, linewidth=2., color='k', linestyle='solid', label=r"cherenkov spectrum (python)")



#cx.plot(wlens, Cherenkov_dN_dXdwlen(wlens, beta)*1.*1e9, linewidth=3., color='k', linestyle='solid', label=r"cherenkov spectrum beta=1")
#cx.plot(wlens, Cherenkov_dN_dXdwlen(wlens, 0.76)*1.*1e9, linewidth=1., color='r', linestyle='solid', label=r"cherenkov spectrum beta=0.76")





ax.set_xlim(260.,690.)
ax.legend()
ax.grid(True)
ax.set_xlabel("wavelength $\\lambda [\\mathrm{nm}]$")
ax.set_ylabel("DOM acceptance")

bx.set_xlim(260.,690.)
bx.legend()
bx.grid(True)
bx.set_xlabel("wavelength $\\lambda^\prime [\\mathrm{nm}]$")
bx.set_ylabel(r"$\int_{\lambda^\prime}^{\lambda^\prime + 1\,\mathrm{nm}} \left( \frac{\mathrm{d}N_\mathrm{phot}}{\mathrm{d}l \mathrm{d}\lambda} \times \epsilon_\mathrm{DOM} \right) \mathrm{d} \lambda$ [$\mathrm{m}^{-1}$]")

cx.set_xlim(260.,690.)
cx.legend()
cx.grid(True)
cx.set_xlabel("wavelength $\\lambda^\prime [\\mathrm{nm}]$")
cx.set_ylabel(r"$\int_{\lambda^\prime}^{\lambda^\prime + 1\,\mathrm{nm}} \left( \frac{\mathrm{d}N_\mathrm{phot}}{\mathrm{d}l \mathrm{d}\lambda} \right) \mathrm{d}\lambda$ [$\mathrm{m}^{-1}$]")




pylab.savefig("dom_acceptance.pdf", transparent=True)



