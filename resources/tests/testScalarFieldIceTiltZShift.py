#!/usr/bin/env python

from __future__ import print_function
import numpy
import math

from icecube import icetray, dataclasses, clsim, phys_services
from I3Tray import I3Units

# test parameters
numberOfTrials = 100000

maximumDeviation = 10.*I3Units.cm

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

iceTiltCLSim = clsim.util.GetIceTiltZShift()

xVals = numpy.random.uniform(-1200.,1200.,numberOfTrials)*I3Units.m
yVals = numpy.random.uniform(-1200.,1200.,numberOfTrials)*I3Units.m
zVals = numpy.random.uniform(-1200.,1200.,numberOfTrials)*I3Units.m


results_OpenCL = evaluateScalarFieldOpenCL(xVals, yVals, zVals, iceTiltCLSim, useReferenceFunction=False)
results_Ref    = evaluateScalarFieldOpenCL(xVals, yVals, zVals, iceTiltCLSim, useReferenceFunction=True)

deviation_OclFromRef = results_OpenCL-results_Ref

maxIndexInOcl = numpy.argmax(numpy.abs(deviation_OclFromRef))
maxDeviationInOcl = deviation_OclFromRef[maxIndexInOcl]

print("maximum deviation in OpenCL implementation:   ", maxDeviationInOcl, "@", maxIndexInOcl, "xyz:", xVals[maxIndexInOcl], yVals[maxIndexInOcl], zVals[maxIndexInOcl], "ref:", results_Ref[maxIndexInOcl], "ocl:", results_OpenCL[maxIndexInOcl])

if numpy.abs(maxDeviationInOcl) > maximumDeviation:
    raise RuntimeError("reference implementation results differ from OpenCL implementation results!")

print("test successful!")


