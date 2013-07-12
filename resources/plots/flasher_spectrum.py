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

def wav2RGB(wavelength):
    w = wavelength/I3Units.nanometer

    # use the first visible color for wavelengths
    # outside the visible spectrum
    
    if w < 380.: w=380.
    if w > 780.: w=780.

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

def genMCHistogramsOpenCL(distribution, range, iterations=1000, numBins=1000):
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
    
    range_width=range[1]-range[0]
    
    num_orig, bins = scipy.histogram(values, range=range, bins=numBins)
    print("hist1 complete")
    
    del values # not needed anymore
    print("deleted")
    
    num=[]
    for number in num_orig:
        num.append(float(number)/float(samples)/float(range_width/float(numBins)))
    num=numpy.array(num)
    
    bins = numpy.array(bins[:-1])+(bins[1]-bins[0])/2.
    
    return dict(num=num, bins=bins)

flasherName = ["340nm", "370nm", "405nm", "450nm", "505nm"]
flasherTypes = [clsim.I3CLSimFlasherPulse.FlasherPulseType.LED340nm,
                clsim.I3CLSimFlasherPulse.FlasherPulseType.LED370nm,
                clsim.I3CLSimFlasherPulse.FlasherPulseType.LED405nm,
                clsim.I3CLSimFlasherPulse.FlasherPulseType.LED450nm,
                clsim.I3CLSimFlasherPulse.FlasherPulseType.LED505nm]

spectrum = []
for i in range(len(flasherName)):
    spectrum.append(clsim.GetIceCubeFlasherSpectrum(spectrumType=flasherTypes[i]))


mediumProps = clsim.MakeIceCubeMediumProperties()
domAcceptance = clsim.GetIceCubeDOMAcceptance()
flatAcceptance = clsim.I3CLSimFunctionConstant(1.)
vectorizedDomAcceptance = numpy.vectorize(lambda x: domAcceptance.GetValue(x))

wlen_range = (mediumProps.GetMinWavelength()/I3Units.nanometer, mediumProps.GetMaxWavelength()/I3Units.nanometer)

spectrumGeneratorNoBias = []
spectrumGeneratorWithBias = []
spectrumData = []
spectrumDataWithBias = []

multiplicatorForBiassedMCSpectrum = []

for i in range(len(flasherName)):
    gen = clsim.makeWavelengthGenerator(spectrum[i], flatAcceptance, mediumProps)
    spectrumGeneratorNoBias.append(gen)

    gen = clsim.makeWavelengthGenerator(spectrum[i], domAcceptance, mediumProps)
    spectrumGeneratorWithBias.append(gen)
    
    vectorizedSpectrum = numpy.vectorize(lambda x: spectrum[i].GetValue(x))
    vectorizedSpectrumWithBias = numpy.vectorize(lambda x: spectrum[i].GetValue(x)*domAcceptance.GetValue(x))
    
    integral_spectrum = scipy.integrate.quad(vectorizedSpectrum, spectrum[i].GetMinWlen(), spectrum[i].GetMaxWlen())[0]/I3Units.nanometer
    print("spectrum integral (w/o bias) [", spectrum[i].GetMinWlen()/I3Units.nanometer, "nm,", spectrum[i].GetMaxWlen()/I3Units.nanometer, "nm] =", integral_spectrum)
    
    bins = numpy.linspace(wlen_range[0]*I3Units.nanometer, wlen_range[1]*I3Units.nanometer, 1000)
    vals = vectorizedSpectrum(bins) / integral_spectrum
    valsWithBias = vectorizedSpectrumWithBias(bins) / integral_spectrum

    correctionFactor = clsim.PhotonNumberCorrectionFactorAfterBias(spectrum[i], domAcceptance, spectrum[i].GetMinWlen(), spectrum[i].GetMaxWlen())
    print("bias photon number correction factor", correctionFactor)

    spectrumData.append([bins, vals])
    spectrumDataWithBias.append([bins, valsWithBias])
    multiplicatorForBiassedMCSpectrum.append(correctionFactor)

biasDataBins = numpy.linspace(wlen_range[0]*I3Units.nanometer, wlen_range[1]*I3Units.nanometer, 10000)
biasDataVals = vectorizedDomAcceptance(biasDataBins)



flasherSpectrumGenNoBias = []
for i in range(len(spectrumGeneratorNoBias)):
    print("preparing spectrum", i, "(no bias)")
    flasherSpectrumGenNoBias.append(genMCHistogramsOpenCL(spectrumGeneratorNoBias[i], range=wlen_range))

flasherSpectrumGenWithBias = []
for i in range(len(spectrumGeneratorWithBias)):
    print("preparing spectrum", i, "(with bias)")
    flasherSpectrumGenWithBias.append(genMCHistogramsOpenCL(spectrumGeneratorWithBias[i], range=wlen_range))


####


fig = pylab.figure(3)
fig.subplots_adjust(left=0.09, bottom=0.05, top=0.95, right=0.98)

ax = fig.add_subplot(3, 1, 1)
bx = fig.add_subplot(3, 1, 2)
cx = fig.add_subplot(3, 1, 3)


wlens=numpy.linspace(wlen_range[0],wlen_range[1],num=10000)
#energies=numpy.linspace(h_times_c/680.,h_times_c/260.,num=10000)

print("plotting...")
legendLabels = []
legendHandles = []

legendLabelsWithBias = []
legendHandlesWithBias = []


for i in range(len(spectrumGeneratorNoBias)):
    peakWlen = spectrumData[i][0][numpy.argmax(spectrumData[i][1])]
    
    peakWlen = spectrumData[i][0][numpy.argmax(spectrumData[i][1])]
    appended=False
    for j in range(len(flasherSpectrumGenNoBias[i]["bins"])-1):
        fromWlen = flasherSpectrumGenNoBias[i]["bins"][j]
        toWlen = flasherSpectrumGenNoBias[i]["bins"][j+1]
        currentWlen = (fromWlen + toWlen)/2.
        currentColor = wav2RGB(currentWlen*I3Units.nanometer)
        h, = ax.plot([fromWlen,toWlen], [flasherSpectrumGenNoBias[i]["num"][j], flasherSpectrumGenNoBias[i]["num"][j+1]], linewidth=2., color=currentColor, linestyle='solid')
        if peakWlen/I3Units.nanometer >= fromWlen and peakWlen/I3Units.nanometer < toWlen and (not appended):
            legendLabels.append(r"``%s'' (peak @ %6.2fnm)" % (flasherName[i], peakWlen/I3Units.nanometer))
            legendHandles.append(h)
            appended=True

    peakWlen = spectrumDataWithBias[i][0][numpy.argmax(spectrumDataWithBias[i][1])]
    appended=False
    for j in range(len(flasherSpectrumGenWithBias[i]["bins"])-1):
        fromWlen = flasherSpectrumGenWithBias[i]["bins"][j]
        toWlen = flasherSpectrumGenWithBias[i]["bins"][j+1]
        currentWlen = (fromWlen + toWlen)/2.
        currentColor = wav2RGB(currentWlen*I3Units.nanometer)
        fromValue = flasherSpectrumGenWithBias[i]["num"][j] * multiplicatorForBiassedMCSpectrum[i]
        toValue = flasherSpectrumGenWithBias[i]["num"][j+1] * multiplicatorForBiassedMCSpectrum[i]
        
        h, = bx.plot([fromWlen,toWlen], [fromValue, toValue], linewidth=2., color=currentColor, linestyle='solid')
        if peakWlen/I3Units.nanometer >= fromWlen and peakWlen/I3Units.nanometer < toWlen and (not appended):
            legendLabelsWithBias.append(r"``%s'' (peak @ %6.2fnm)" % (flasherName[i], peakWlen/I3Units.nanometer))
            legendHandlesWithBias.append(h)
            appended=True

    for j in range(len(spectrumData[i][0])-1):
        fromWlen = spectrumData[i][0][j]/I3Units.nanometer
        toWlen = spectrumData[i][0][j+1]/I3Units.nanometer
        currentWlen = (fromWlen + toWlen)/2.
        currentColor = wav2RGB(currentWlen*I3Units.nanometer)
        currentColorDark = (currentColor[0]/2., currentColor[1]/2., currentColor[2]/2.)
        
        ax.plot([fromWlen,toWlen], [spectrumData[i][1][j], spectrumData[i][1][j+1]], linewidth=1., color=currentColorDark, linestyle='solid')

    for j in range(len(spectrumDataWithBias[i][0])-1):
        fromWlen = spectrumDataWithBias[i][0][j]/I3Units.nanometer
        toWlen = spectrumDataWithBias[i][0][j+1]/I3Units.nanometer
        currentWlen = (fromWlen + toWlen)/2.
        currentColor = wav2RGB(currentWlen*I3Units.nanometer)
        currentColorDark = (currentColor[0]/2., currentColor[1]/2., currentColor[2]/2.)
        
        bx.plot([fromWlen,toWlen], [spectrumDataWithBias[i][1][j], spectrumDataWithBias[i][1][j+1]], linewidth=1., color=currentColorDark, linestyle='solid')


biasHandle, = bx.plot(biasDataBins/I3Units.nanometer, biasDataVals*0.06, linewidth=2., color='k', linestyle='solid')
biasLabel = r"0.06x DOM wavelength acceptance"

ax.set_xlim(260.,690.)
ax.set_ylim(0.,0.08)
ax.legend(legendHandles, legendLabels)
ax.grid(True)
ax.set_xlabel("wavelength $\\lambda [\\mathrm{nm}]$")
ax.set_ylabel("intensity [a.u.]")


bx.set_xlim(260.,690.)
#bx.set_ylim(0.,0.08)
bx.legend(legendHandlesWithBias + [biasHandle], legendLabelsWithBias + [biasLabel])
bx.grid(True)
bx.set_xlabel("wavelength $\\lambda [\\mathrm{nm}]$")
bx.set_ylabel("intensity [a.u.]")


print("saving..")
pylab.savefig("flasher_spectrum.pdf", transparent=True)

print("done.")


