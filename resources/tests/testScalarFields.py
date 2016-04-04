#!/usr/bin/env python

from __future__ import print_function
import numpy
import math

from icecube import icetray, dataclasses, clsim, phys_services
from I3Tray import I3Units

# test parameters
numberOfTrials = 100000
anisotropyDirAzimuthInDeg = 216.
magnitudeAlongDir = 0.04
magnitudePerpToDir = -0.08

maximumRelativeDeviation = 1e-5

# get OpenCL devices
openCLDevices = [device for device in clsim.I3CLSimOpenCLDevice.GetAllDevices()]
if len(openCLDevices)==0:
    raise RuntimeError("No OpenCL devices available!")
openCLDevice = openCLDevices[0]

openCLDevice.useNativeMath=False
workgroupSize = 1
workItemsPerIteration = 10240
print("           using platform:", openCLDevice.platform)
print("             using device:", openCLDevice.device)
print("            workgroupSize:", workgroupSize)
print("    workItemsPerIteration:", workItemsPerIteration)


def evaluateScalarFieldOpenCL(xValues, yValues, zValues, scalarField, useReferenceFunction=False):
    tester = clsim.I3CLSimScalarFieldTester(
        device=openCLDevice,
        workgroupSize=workgroupSize,
        workItemsPerIteration=workItemsPerIteration,
        scalarField=scalarField)
    
    # the function currently only accepts I3VectorFloat as its input type
    vectorX = dataclasses.I3VectorFloat(numpy.array(xValues))
    vectorY = dataclasses.I3VectorFloat(numpy.array(yValues))
    vectorZ = dataclasses.I3VectorFloat(numpy.array(zValues))

    if useReferenceFunction:
        yValues = tester.EvaluateReferenceFunction(vectorX, vectorY, vectorZ)
    else:
        yValues = tester.EvaluateFunction(vectorX, vectorY, vectorZ)

    return numpy.array(yValues)

def DimasAbsLenScalingFactor(x,y,z, thx,logk1,logk2):
    def _DimasAbsLenScalingFactor(
        x,y,z,
        thx=216.,    # direction of ice tilt (perp. to flow)  (thx)
        logk1=0.04,  # magnitude of ice anisotropy along tilt (logk1)
        logk2=-0.08, # magnitude of ice anisotropy along flow (logk2)
        ):

        azx = math.cos(thx*math.pi/180.) # x-component of the "tilt" direction
        azy = math.sin(thx*math.pi/180.) # y-component of the "tilt" direction
        k1  = math.exp(logk1)            # coefficient of anisotropy parallel to "tilt" direction
        k2  = math.exp(logk2)            # coefficient of anisotropy perpendicular to "tilt" direction
        kz  = 1./(k1*k2)                 # a normalizing factor for the z direction

        # stolen from ppc and converted to python
        nr=1.

        n1= azx*x+azy*y
        n2=-azy*x+azx*y
        n3= z

        s1=n1*n1
        l1=k1*k1
        s2=n2*n2
        l2=k2*k2
        s3=n3*n3
        l3=kz*kz

        B2=nr/l1+nr/l2+nr/l3
        nB=s1/l1+s2/l2+s3/l3
        An=s1*l1+s2*l2+s3*l3

        nr=(B2-nB)*An/2.

        return 1./nr
    
    res = numpy.zeros(len(x))
    for i in range(len(x)):
        res[i] = _DimasAbsLenScalingFactor(x[i], y[i], z[i], thx,logk1,logk2)
    return res

DimasScalarFieldOpenCL = clsim.I3CLSimScalarFieldAnisotropyAbsLenScaling(anisotropyDirAzimuth=anisotropyDirAzimuthInDeg*I3Units.deg, magnitudeAlongDir=magnitudeAlongDir, magnitudePerpToDir=magnitudePerpToDir)

zeniths = numpy.arccos(numpy.random.uniform(0.,1.,numberOfTrials)*2. - 1.)
azimuths = numpy.random.uniform(0.,2.*math.pi,numberOfTrials)
xVals = numpy.sin(zeniths)*numpy.cos(azimuths)
yVals = numpy.sin(zeniths)*numpy.sin(azimuths)
zVals = numpy.cos(zeniths)


results_OpenCL = evaluateScalarFieldOpenCL(xVals, yVals, zVals, DimasScalarFieldOpenCL, useReferenceFunction=False)
results_Ref    = evaluateScalarFieldOpenCL(xVals, yVals, zVals, DimasScalarFieldOpenCL, useReferenceFunction=True)
results_Python = DimasAbsLenScalingFactor(xVals, yVals, zVals, thx=anisotropyDirAzimuthInDeg, logk1=magnitudeAlongDir, logk2=magnitudePerpToDir)


deviation_RefFromPython = (results_Python-results_Ref)/results_Python
deviation_OclFromPython = (results_Python-results_OpenCL)/results_Python

maxIndexInRef = numpy.argmax(numpy.abs(deviation_RefFromPython))
maxIndexInOcl = numpy.argmax(numpy.abs(deviation_OclFromPython))
maxDeviationInRef = deviation_RefFromPython[maxIndexInRef]
maxDeviationInOcl = deviation_OclFromPython[maxIndexInOcl]

#print deviation_OclFromPython[10230:10300]

print("maximum relative deviation in reference implementation:", maxDeviationInRef, "@", maxIndexInRef, "xyz:", xVals[maxIndexInRef], yVals[maxIndexInRef], zVals[maxIndexInRef], "python:", results_Python[maxIndexInRef], "ref:", results_Ref[maxIndexInRef], "ocl:", results_OpenCL[maxIndexInRef])
print("maximum relative deviation in OpenCL implementation:   ", maxDeviationInOcl, "@", maxIndexInOcl, "xyz:", xVals[maxIndexInOcl], yVals[maxIndexInOcl], zVals[maxIndexInOcl], "python:", results_Python[maxIndexInOcl], "ref:", results_Ref[maxIndexInOcl], "ocl:", results_OpenCL[maxIndexInOcl])

#print "minimum relative deviation in reference implementation:", minDeviationInRef
#print "minimum relative deviation in OpenCL implementation:   ", minDeviationInOcl

if numpy.abs(maxDeviationInRef) > maximumRelativeDeviation:
    raise RuntimeError("python implementation results differ from C++ reference implementation results!")

if numpy.abs(maxDeviationInOcl) > maximumRelativeDeviation:
    raise RuntimeError("python implementation results differ from OpenCL implementation results!")

print("test successful!")


