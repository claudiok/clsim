from icecube.load_pybindings import load_pybindings
from icecube import icetray, dataclasses # be nice and pull in our dependencies
load_pybindings(__name__,__path__)

# this is all a hugh hack so I do not have to change the dataclasses
def augmentClassWithPicklingSupport(theClass, setstateFunc, getstateFunc):
    theClass.__setstate__ = setstateFunc
    theClass.__getstate__ = getstateFunc
    theClass.__safe_for_unpickling__ = True

augmentClassWithPicklingSupport(icetray.I3Frame, I3Frame___setstate__, I3Frame___getstate__)
augmentClassWithPicklingSupport(dataclasses.I3Particle, I3Particle___setstate__, I3Particle___getstate__)
augmentClassWithPicklingSupport(dataclasses.I3Direction, I3Direction___setstate__, I3Direction___getstate__)
augmentClassWithPicklingSupport(dataclasses.I3Position, I3Position___setstate__, I3Position___getstate__)


from MakeAntaresMediumProperties import GetPetzoldScatteringCosAngleDistribution, GetAntaresScatteringCosAngleDistribution, MakeAntaresMediumProperties
from MakeIceCubeMediumProperties import MakeIceCubeMediumProperties
from I3CLSimParticleToStepConverterGeant4MP import I3CLSimParticleToStepConverterGeant4MP


# clean up the clsim namespace
del augmentClassWithPicklingSupport
del icetray
del dataclasses

