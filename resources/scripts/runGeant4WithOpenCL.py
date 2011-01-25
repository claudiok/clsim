#!/usr/bin/env python

import time
import random
import numpy
import math
import string

from icecube import icetray, dataclasses, dataio, phys_services, clsim
from I3Tray import I3Units

def initializeOpenCL(rng, geometry, medium):
    conv = clsim.I3CLSimStepToPhotonConverterOpenCL(RandomService=rng,
                                                    UseNativeMath=True,
                                                    CPUOnly=False)
    
    print "available OpenCL devices:"
    deviceList = conv.GetDeviceList()
    for device in deviceList:
        print "platform: \"%s\", device: \"%s\"" % (device[0], device[1])
        
    print ""
    
    # do a semi-smart device selection
    deviceToUse=None
    
    # look for a "GeForce" first
    geForceDevices = [device for device in deviceList if string.lower(device[1]).find("geforce")!=-1]
    if len(geForceDevices)>0:
        if len(geForceDevices)>1:
            deviceToUse=geForceDevices[0]
            print "You seem to have more than one GeForce GPU. Using the first one. (\"%s\")" % (deviceToUse[1])
        else:
            deviceToUse=geForceDevices[0]
            print "You seem to have a GeForce GPU. (\"%s\")" % (deviceToUse[1])
    else:
        print "No GeForce device found. Just using the first available one."
        deviceToUse=deviceList[0]
    
    print " -> using OpenCL device \"%s\" on platform \"%s\"." % (deviceToUse[1], deviceToUse[0])
    
    conv.SetDeviceName(deviceToUse[0], deviceToUse[1])

    conv.SetMediumProperties(medium)
    conv.SetGeometry(geometry)

    conv.Compile()
    #print conv.GetFullSource()

    conv.workgroupSize = conv.maxWorkgroupSize
    
    # use approximately 512000 work items, convert to a multiple of the workgroup size
    conv.maxNumWorkitems = (1024000/conv.workgroupSize)*conv.workgroupSize

    print "maximum workgroup size is", conv.maxWorkgroupSize
    print "configured workgroup size is", conv.workgroupSize
    print "maximum number of work items is", conv.maxNumWorkitems

    conv.Initialize()

    return conv
    
def initializeGeant4(rng, medium, openCLconverter):
    conv = clsim.I3CLSimParticleToStepConverterGeant4()

    conv.SetMediumProperties(medium)
    conv.SetMaxBunchSize(openCLconverter.maxNumWorkitems)
    conv.SetBunchSizeGranularity(openCLconverter.workgroupSize)
    
    conv.Initialize()
    
    return conv
    
# geneartes an example I3Particle
def generateRandomParticle(rng):
    part = dataclasses.I3Particle()
    part.SetPos(rng.Uniform(-200.*I3Units.m,200.*I3Units.m), rng.Uniform(-200.*I3Units.m,200.*I3Units.m), rng.Uniform(-200.*I3Units.m,200.*I3Units.m))
    part.SetDir(numpy.arccos(rng.Uniform(-1.,1.)), rng.Uniform(0.,2.*math.pi))
    part.SetTime(0.)
    part.SetEnergy(1.*I3Units.GeV)
    part.SetType(dataclasses.I3Particle.PiPlus)
    return part
    
# we need random numbers
rng = phys_services.I3SPRNGRandomService(seed=536781, nstreams=1000, streamnum=0)

geometry = clsim.I3CLSimSimpleGeometryTextFile(43.18*I3Units.cm/2., "resources/Icecube_geometry.20100310.complete.txt")
medium = clsim.MakeIceCubeMediumProperties()
#medium = clsim.MakeAntaresMediumProperties()

# initialize OpenCL
openCLStepsToPhotonsConverter = initializeOpenCL(rng, geometry, medium)

# initialize Geant4 (will set bunch sizes according to the OpenCL settings)
geant4ParticleToStepsConverter = initializeGeant4(rng, medium, openCLStepsToPhotonsConverter)


numParticles = 100000

print "sending %u particles to Geant4" % (numParticles)

# send numParticles particles to Geant4
for i in range(numParticles):
    geant4ParticleToStepsConverter.EnqueueParticle(generateRandomParticle(rng), i)
geant4ParticleToStepsConverter.EnqueueBarrier()

print "particles sent to Geant4."

print "receiving steps from Geant4 and sending to OpenCL.."

numBunchesSentToOpenCL=0

while geant4ParticleToStepsConverter.BarrierActive() or geant4ParticleToStepsConverter.MoreStepsAvailable():
    steps = geant4ParticleToStepsConverter.GetConversionResult()
    
    if isinstance(steps, tuple) and isinstance(s[1], dataclasses.I3Particle):
        print "Got a secondary particle from Geant4. This was not configured, something is wrong. ignoring."
        continue
    elif not isinstance(steps, clsim.I3CLSimStepSeries):
        print "Got something unknown from Geant4. ignoring."
        continue
    
    # has to be a bunch of steps at this point
    #print "Got", len(steps), "steps from Geant4, sending them to OpenCL"
    
    openCLStepsToPhotonsConverter.EnqueueSteps(steps, numBunchesSentToOpenCL)
    numBunchesSentToOpenCL += 1
    
    #print "steps sent to OpenCL."
    

# once everything is sent, retrieve the results from OpenCL
print "Geant4 finished. Getting results from OpenCL:"

allPhotons = clsim.I3CLSimPhotonSeries()

for i in range(numBunchesSentToOpenCL):
#while openCLStepsToPhotonsConverter.MorePhotonsAvailable():
    res = openCLStepsToPhotonsConverter.GetConversionResult()
    photons = res[1]
    #print "result #", res[0], ": num photons =", len(res[1])
    if (len(photons)>0):
        allPhotons.extend(photons) # append the photons to the list of all photons
    del photons
    del res

print "results fetched from OpenCL."

print "Got", len(allPhotons), "photons in total for the initial", numParticles, "particles."

print "finished."



