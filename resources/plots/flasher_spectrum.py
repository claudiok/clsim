#!/usr/bin/env python

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

def wav2RGB(wavelength):
    w = wavelength/I3Units.nanometer

    # colour
    if w >= 380. and w < 440.:
        R = -(w - 440.) / (440. - 350.)
        G = 0.0
        B = 1.0
    elif w >= 440. and w < 490.:
        R = 0.0
        G = (w - 440.) / (490. - 440.)
        B = 1.0
    elif w >= 490. and w < 510.:
        R = 0.0
        G = 1.0
        B = -(w - 510.) / (510. - 490.)
    elif w >= 510. and w < 580.:
        R = (w - 510.) / (580. - 510.)
        G = 1.0
        B = 0.0
    elif w >= 580. and w < 645.:
        R = 1.0
        G = -(w - 645.) / (645. - 580.)
        B = 0.0
    elif w >= 645. and w <= 780.:
        R = 1.0
        G = 0.0
        B = 0.0
    else:
        R = 0.0
        G = 0.0
        B = 0.0

    # intensity correction
    if w >= 380. and w < 420.:
        SSS = 0.3 + 0.7*(w - 350) / (420 - 350)
    elif w >= 420. and w <= 700.:
        SSS = 1.0
    elif w > 700. and w <= 780.:
        SSS = 0.3 + 0.7*(780 - w) / (780 - 700)
    else:
        SSS = 0.0

    return (SSS*R, SSS*G, SSS*B)


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
print "           using platform:", openCLDevice.platform
print "             using device:", openCLDevice.device
print "            workgroupSize:", workgroupSize
print "    workItemsPerIteration:", workItemsPerIteration



def applyOpenCLWlenDependentFunction(xValues, functionOpenCL, getDerivative=False, useReferenceFunction=False):
    #print "         number of values:", len(xValues)
    
    tester = clsim.I3CLSimWlenDependentValueTester(device=openCLDevice,
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

def genMCHistogramsOpenCL(distribution, range, iterations=1000, numBins=1000):
    tester = clsim.I3CLSimRandomDistributionTester(device=openCLDevice,
                                                   workgroupSize=workgroupSize,
                                                   workItemsPerIteration=workItemsPerIteration,
                                                   randomService=rng,
                                                   randomDistribution=distribution)
    
    values = tester.GenerateRandomNumbers(iterations)
    samples = len(values)
    print "generated"
    
    values = numpy.array(values)/I3Units.nanometer # convert to numpy array and convert units
    print "converted"
    
    range_width=range[1]-range[0]
    
    num_orig, bins = scipy.histogram(values, range=range, bins=numBins)
    print "hist1 complete"
    
    del values # not needed anymore
    print "deleted"
    
    num=[]
    for number in num_orig:
        num.append(float(number)/float(samples)/float(range_width/float(numBins)))
    num=numpy.array(num)
    
    bins = numpy.array(bins[:-1])+(bins[1]-bins[0])/2.
    
    return dict(num=num, bins=bins)


spectrumData = []
spectrumData.append(clsim.GetIceCubeFlasherSpectrumData(colorIndex=0))

spectrumGenerator = []
spectrumGenerator.append(clsim.GetIceCubeFlasherSpectrumGenerator(colorIndex=0))

# normalize
for i in range(len(spectrumData)):
    integral = scipy.integrate.trapz(y=spectrumData[i][1], x=spectrumData[i][0]/I3Units.nanometer)
    spectrumData[i][1] /= integral


mediumProps = clsim.MakeIceCubeMediumProperties()
domAcceptance = clsim.GetIceCubeDOMAcceptance()
flatAcceptance = clsim.I3CLSimWlenDependentValueConstant(1.)

wlen_range = (mediumProps.GetMinWavelength()/I3Units.nanometer, mediumProps.GetMaxWavelength()/I3Units.nanometer)

#genWavelength = clsim.makeWavelengthGenerator(domAcceptance, False, mediumProps)
#histGenWavelength = genMCHistogramsOpenCL(genWavelength, range=wlen_range)

#genWavelengthFlat = clsim.makeWavelengthGenerator(flatAcceptance, False, mediumProps)
#histGenWavelengthFlat = genMCHistogramsOpenCL(genWavelengthFlat, range=wlen_range)

flasherSpectrumGen = []
for i in range(len(spectrumGenerator)):
    print "preparing spectrum", i
    flasherSpectrumGen.append(genMCHistogramsOpenCL(spectrumGenerator[0], range=wlen_range))


####


fig = pylab.figure(3)
fig.subplots_adjust(left=0.09, bottom=0.05, top=0.95, right=0.98)

ax = fig.add_subplot(3, 1, 1)
bx = fig.add_subplot(3, 1, 2)
cx = fig.add_subplot(3, 1, 3)


wlens=numpy.linspace(wlen_range[0],wlen_range[1],num=10000)
#energies=numpy.linspace(h_times_c/680.,h_times_c/260.,num=10000)

peakWlen = spectrumData[0][0][numpy.argmax(spectrumData[0][1])]
color = wav2RGB(peakWlen)
colorDark = (color[0]/2., color[1]/2., color[2]/2.)
ax.plot(flasherSpectrumGen[0]["bins"], flasherSpectrumGen[0]["num"], linewidth=2., color=color, label=r"``405nm'' flasher spectrum (OpenCL MC generated) peak @ %6.2fnm" % (peakWlen/I3Units.nanometer))
ax.plot(spectrumData[0][0]/I3Units.nanometer, spectrumData[0][1], linewidth=1., color=colorDark, linestyle='solid', label=r"``405nm'' flasher spectrum (python) peak @ %6.2fnm" % (peakWlen/I3Units.nanometer))



ax.set_xlim(260.,690.)
ax.set_ylim(0.,0.08)
ax.legend()
ax.grid(True)
ax.set_xlabel("wavelength $\\lambda [\\mathrm{nm}]$")
ax.set_ylabel("DOM acceptance")



pylab.savefig("flasher_spectrum.pdf", transparent=True)



