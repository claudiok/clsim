#!/usr/bin/env python

from icecube import icetray, dataclasses, dataio, clsim
from I3Tray import I3Units
conv = clsim.I3CLSimParticleToStepConverterGeant4(randomSeed=748239)

conv.SetMediumProperties(clsim.MakeIceCubeMediumProperties())
conv.SetMaxBunchSize(512000)
conv.SetBunchSizeGranularity(512)

## the module can output em cascades as secondaries
## for a certain secondary energy range
#conv.electronPositronMinEnergyForSecondary = 500.*I3Units.MeV
#conv.electronPositronMaxEnergyForSecondary = 100.*I3Units.GeV

conv.Initialize()

part = dataclasses.I3Particle()
part.SetDir(1.,0.,0.)
part.SetPos(0.,0.,0.)
part.SetTime(0.)
part.SetEnergy(1.*I3Units.GeV)
#part.SetType(dataclasses.I3Particle.EMinus)
part.SetType(dataclasses.I3Particle.PiPlus)
#part.SetType(dataclasses.I3Particle.MuPlus)

for i in range(1000):
    conv.EnqueueParticle(part, i)
conv.EnqueueBarrier()

totalEnergyInEMCascades=0.
totalNumEMCascades=0
totalSteps=0
while conv.BarrierActive() or conv.MoreStepsAvailable():
    s = conv.GetConversionResult()
    if isinstance(s, tuple) and isinstance(s[1], dataclasses.I3Particle):
        #print "Got a particle from Geant4 with id", s[0], "and energy", s[1].GetEnergy()/I3Units.GeV, "GeV"
        totalEnergyInEMCascades += s[1].GetEnergy()
        totalNumEMCascades+=1
    elif isinstance(s, clsim.I3CLSimStepSeries):
        print "Got", len(s), "steps from Geant4"
        totalSteps+=len(s)
    else:
        print "Got something unknown from Geant4."

print "finished!"

print "Got", totalSteps, "Cherenkov steps and a total of", totalEnergyInEMCascades/I3Units.GeV, "GeV in", totalNumEMCascades, "secondary em cascades"
