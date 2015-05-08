#!/usr/bin/env python

from icecube.icetray import I3Units

DOMOversizeFactor = 5
DOMRadius = 0.16510*I3Units.m # 13" diameter

from icecube import icetray, clsim, dataclasses, phys_services
from os.path import expandvars
rng = phys_services.I3GSLRandomService(0)

ppcConverter = clsim.I3CLSimLightSourceToStepConverterPPC(photonsPerStep=200)
particleParameterizations = clsim.GetDefaultParameterizationList(ppcConverter, muonOnly=False)
mediumProperties = clsim.traysegments.common.parseIceModel(expandvars("$I3_SRC/clsim/resources/ice/spice_mie"), disableTilt=False)
domAcceptance = clsim.GetIceCubeDOMAcceptance(domRadius = DOMRadius*DOMOversizeFactor, efficiency=1)

def makeGeometry():
    
    frame = icetray.I3Frame(icetray.I3Frame.Geometry)
    
    omgeo = dataclasses.I3OMGeo()
    omgeo.omtype = dataclasses.I3OMGeo.OMType.IceCube
    omgeo.orientation = dataclasses.I3Orientation(dataclasses.I3Direction(0.,0.,-1.))
    omgeo.position = dataclasses.I3Position(0, 0, 0)
    omgeomap = dataclasses.I3OMGeoMap()
    omgeomap[icetray.OMKey(1,1)] = omgeo
    frame['I3OMGeoMap'] = omgeomap
    
    modgeo = dataclasses.I3ModuleGeo()
    modgeo.module_type = dataclasses.I3ModuleGeo.ModuleType.IceCube
    modgeo.orientation = omgeo.orientation
    modgeo.pos = omgeo.position
    modgeomap = dataclasses.I3ModuleGeoMap()
    modgeomap[dataclasses.ModuleKey(1,1)] = modgeo
    frame['I3ModuleGeoMap'] = modgeomap
    
    subdetectors = dataclasses.I3MapModuleKeyString()
    for k in modgeomap.keys():
        subdetectors[k] = "IceCube"
    frame['Subdetectors'] = subdetectors
    
    return clsim.I3CLSimSimpleGeometryFromI3Geometry(DOMRadius, DOMOversizeFactor, frame)

geant = clsim.I3CLSimLightSourceToStepConverterGeant4()
geant.SetRandomService(rng)
geant.SetWlenBias(domAcceptance)
geant.SetMediumProperties(mediumProperties)
geant.SetGeometry(makeGeometry())
geant.SetLightSourceParameterizationSeries(particleParameterizations)
geant.Initialize()

p = dataclasses.I3Particle()
p.type = p.EMinus
p.time = 0
p.energy = 1
p.pos = dataclasses.I3Position(1,0,0)
p.dir = dataclasses.I3Direction(0,0,1)
p.length = 0

# a really closeby cascade should be treated with no oversizing
geant.EnqueueLightSource(clsim.I3CLSimLightSource(p), 0)
geant.EnqueueBarrier()
steps = geant.GetConversionResult()
for step in steps:
    if step.num > 0:
        assert step.undersizeFactor == DOMOversizeFactor

# a middling cascade should get middling oversizing
p.pos = dataclasses.I3Position(DOMOversizeFactor*3/2.,0,0)
geant.EnqueueLightSource(clsim.I3CLSimLightSource(p), 0)
geant.EnqueueBarrier()
steps = geant.GetConversionResult()
for step in steps:
    if step.num > 0:
        assert step.undersizeFactor > 1. and step.undersizeFactor < DOMOversizeFactor

# a faraway cascade should get the full oversizing
p.pos = dataclasses.I3Position(DOMOversizeFactor*3,0,0)
geant.EnqueueLightSource(clsim.I3CLSimLightSource(p), 0)
geant.EnqueueBarrier()
steps = geant.GetConversionResult()
for step in steps:
    if step.num > 0:
        assert step.undersizeFactor == 1.

