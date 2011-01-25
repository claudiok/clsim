#!/usr/bin/env python

import time
import random
import numpy
import math

from icecube import icetray, dataclasses, dataio, phys_services, clsim
from I3Tray import I3Units

# we need random numbers
rng = phys_services.I3SPRNGRandomService(seed=536781, nstreams=1000, streamnum=0)


conv = clsim.I3CLSimStepToPhotonConverterOpenCL(RandomService=rng,
                                                UseNativeMath=True)
print conv.GetDeviceList()

conv.SetDeviceName("NVIDIA CUDA", "GeForce GTX 580")
#conv.SetDeviceName("ATI Stream", "Intel(R) Core(TM) i5 CPU         760  @ 2.80GHz")
#conv.SetDeviceName("Apple", "GeForce 9600M GT")
#conv.SetDeviceName("Apple", "GeForce 9400M")
#conv.SetDeviceName("Apple", "Intel(R) Core(TM)2 Duo CPU     T9600  @ 2.80GHz")

conv.SetMediumProperties(clsim.MakeIceCubeMediumProperties())
#conv.SetMediumProperties(clsim.MakeAntaresMediumProperties())
conv.SetGeometry(clsim.I3CLSimSimpleGeometryTextFile(43.18*I3Units.cm/2., "resources/Icecube_geometry.20100310.complete.txt"))

conv.Compile()
#print conv.GetFullSource()

conv.workgroupSize = conv.maxWorkgroupSize
conv.maxNumWorkitems = conv.maxWorkgroupSize * 3000

photonsPerStep = 200
repetitions = 10

print "maximum workgroup size is", conv.maxWorkgroupSize
print "configured workgroup size is", conv.workgroupSize
print "maximum number of work items is", conv.maxNumWorkitems

conv.Initialize()

numSteps = conv.maxNumWorkitems
print "making", numSteps, "fake steps"

myStep = clsim.I3CLSimStep()
myStep.x = 20.*I3Units.m
myStep.y = -50.*I3Units.m
myStep.z = 70.*I3Units.m
myStep.SetDirXYZ(1.,1.,0.) # will automatically be normalized
myStep.beta = 1.
myStep.length = 1.*I3Units.mm
myStep.id = 0
myStep.num = photonsPerStep
myStep.weight = 1.
myStep.time = 0.*I3Units.ns

steps = clsim.I3CLSimStepSeries()
for i in range(numSteps):
    steps.append(myStep)
    steps[-1].theta = numpy.arccos(rng.Uniform(-1.,1.))
    steps[-1].phi = rng.Uniform(0.,2.*math.pi)
    steps[-1].x = rng.Uniform(-200.*I3Units.m,200.*I3Units.m)
    steps[-1].y = rng.Uniform(-200.*I3Units.m,200.*I3Units.m)
    steps[-1].z = rng.Uniform(-200.*I3Units.m,200.*I3Units.m)
    steps[-1].id = i

print "sending steps to OpenCL"

starttime = time.time()

for i in range(repetitions):
    conv.EnqueueSteps(steps, i)
for i in range(repetitions):
    res = conv.GetConversionResult()
    photons = res[1]
    print "result #", res[0], ": num photons =", len(res[1])
    if (len(photons)>0):
        print photons[0]
        print "..."
        print photons[-1]
    else:
        print "(empty)"

endtime = time.time()

duration = endtime-starttime
totNumSteps = numSteps*repetitions
totNumPhotons = totNumSteps*photonsPerStep

print "took", duration, "seconds"
print " =>", float(duration*1e9)/float(totNumSteps), "ns per step"
print " =>", float(duration*1e9)/float(totNumPhotons), "ns per photon"



