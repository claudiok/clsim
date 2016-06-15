#!/usr/bin/env python

from __future__ import print_function
import numpy
import math

from icecube import icetray, dataclasses, clsim, phys_services
from I3Tray import I3Units

# test parameters
numberOfTrials = 100000
renormalizeUnitVector = True

maximumRelativeDeviation = 1e-4

theMatrix = numpy.random.uniform(-10.,10.,(3,3))

# get OpenCL CPU devices
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


def evaluateVectorTransformOpenCL(xValues, yValues, zValues, VectorTransform, useReferenceFunction=False):
    tester = clsim.I3CLSimVectorTransformTester(
        device=openCLDevice,
        workgroupSize=workgroupSize,
        workItemsPerIteration=workItemsPerIteration,
        vectorTransform=VectorTransform)
    
    # the function currently only accepts I3VectorFloat as its input type
    vectorX = dataclasses.I3VectorFloat(numpy.array(xValues))
    vectorY = dataclasses.I3VectorFloat(numpy.array(yValues))
    vectorZ = dataclasses.I3VectorFloat(numpy.array(zValues))

    if useReferenceFunction:
        retX, retY, retZ = tester.EvaluateReferenceFunction(vectorX, vectorY, vectorZ)
    else:
        retX, retY, retZ = tester.EvaluateFunction(vectorX, vectorY, vectorZ)

    return numpy.array([retX, retY, retZ]).T

def calculateDotProducts(matrix, xVals, yVals, zVals, renormalize=False):
    retval = []
    for i in range(len(xVals)):
        thisVec = numpy.dot(matrix, numpy.array([xVals[i], yVals[i], zVals[i]]))
        if renormalize:
            thisVec = thisVec/numpy.sqrt(numpy.sum(thisVec**2))
        retval.append(thisVec)

    return numpy.array(retval)

MatrixVectorTransformOpenCL = clsim.I3CLSimVectorTransformMatrix(dataclasses.I3Matrix(theMatrix), renormalize=renormalizeUnitVector)

zeniths = numpy.arccos(numpy.random.uniform(0.,1.,numberOfTrials)*2. - 1.)
azimuths = numpy.random.uniform(0.,2.*math.pi,numberOfTrials)
xVals = numpy.sin(zeniths)*numpy.cos(azimuths)
yVals = numpy.sin(zeniths)*numpy.sin(azimuths)
zVals = numpy.cos(zeniths)


results_OpenCL = numpy.ravel(evaluateVectorTransformOpenCL(xVals, yVals, zVals, MatrixVectorTransformOpenCL, useReferenceFunction=False))
results_Ref    = numpy.ravel(evaluateVectorTransformOpenCL(xVals, yVals, zVals, MatrixVectorTransformOpenCL, useReferenceFunction=True))
results_Python = numpy.ravel(calculateDotProducts(theMatrix, xVals, yVals, zVals, renormalize=renormalizeUnitVector))



deviation_RefFromPython = results_Python-results_Ref
deviation_OclFromPython = results_Python-results_OpenCL

maxIndexInRef = numpy.argmax(numpy.abs(deviation_RefFromPython))
maxIndexInOcl = numpy.argmax(numpy.abs(deviation_OclFromPython))
maxDeviationInRef = deviation_RefFromPython[maxIndexInRef]
maxDeviationInOcl = deviation_OclFromPython[maxIndexInOcl]

minIndexInRef = numpy.argmin(numpy.abs(deviation_RefFromPython))
minIndexInOcl = numpy.argmin(numpy.abs(deviation_OclFromPython))
minDeviationInRef = deviation_RefFromPython[minIndexInRef]
minDeviationInOcl = deviation_OclFromPython[minIndexInOcl]

print("maximum absolute deviation in reference implementation:", maxDeviationInRef, "@", maxIndexInRef, "python:", results_Python[maxIndexInRef], "ref:", results_Ref[maxIndexInRef], "ocl:", results_OpenCL[maxIndexInRef])
print("maximum absolute deviation in OpenCL implementation:   ", maxDeviationInOcl, "@", maxIndexInOcl, "python:", results_Python[maxIndexInOcl], "ref:", results_Ref[maxIndexInOcl], "ocl:", results_OpenCL[maxIndexInOcl])

print("minimum relative deviation in reference implementation:", minDeviationInRef)
print("minimum relative deviation in OpenCL implementation:   ", minDeviationInOcl)

if numpy.abs(maxDeviationInRef) > maximumRelativeDeviation:
    raise RuntimeError("python implementation results differ from C++ reference implementation results!")

if numpy.abs(maxDeviationInOcl) > maximumRelativeDeviation:
    raise RuntimeError("python implementation results differ from OpenCL implementation results!")

print("test successful!")


