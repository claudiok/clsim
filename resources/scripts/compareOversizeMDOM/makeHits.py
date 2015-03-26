#!/usr/bin/env python

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-s", "--seed", type=int, default=1)
parser.add_option("--oversize", type=float, default=1)
parser.add_option("--unweighted", default=False, action="store_true")
parser.add_option("--no-gpu", default=True, dest="gpu", action="store_false")

opts, args = parser.parse_args()

infiles, outfile = args[:-1], args[-1]

from icecube import icetray, dataclasses, dataio, phys_services, sim_services
from icecube.icetray import I3Units
from I3Tray import I3Tray

tray = I3Tray()

randomService = phys_services.I3GSLRandomService(opts.seed)
tray.context['I3RandomService'] = randomService

tray.Add("I3Reader", filenamelist=infiles)

icetray.logging.set_level_for_unit('clsim', 'INFO')

from icecube import clsim
import efficiencies

kwargs = dict()

from icecube.clsim.GetIceCubeDOMAcceptance import GetIceCubeDOMAcceptance

if opts.oversize != 1.:
	# kwargs['WavelengthGenerationBias'] = efficiencies.GetAcceptanceEnvelope(oversize=opts.oversize)
	kwargs['WavelengthGenerationBias'] = efficiencies.GetMDOMAcceptance(oversize=opts.oversize, efficiency=1.95*0.99)
	# kwargs['DoNotParallelize'] = True
	# kwargs['WavelengthGenerationBias'] = clsim.I3CLSimFunctionConstant(1e-2)
	# kwargs['WavelengthGenerationBias'] = efficiencies.GetIceCubeDOMAcceptance(oversize=opts.oversize, efficiency=5)
	# kwargs['WavelengthGenerationBias'] = GetIceCubeDOMAcceptance(opts.oversize*0.16510*icetray.I3Units.m, 0.99*1.35*0.75*1.01)

	# import pylab, numpy
	# wvl = numpy.arange(260, 680, 10)*I3Units.nanometer
	# pylab.plot(wvl, numpy.vectorize(kwargs['WavelengthGenerationBias'].GetValue)(wvl))
	#
	# # kwargs['WavelengthGenerationBias'] = efficiencies.GetIceCubeDOMAcceptance(opts.oversize, 0.99*1.35*0.75*1.01)
	#
	#
	# pylab.plot(wvl, numpy.vectorize(kwargs['WavelengthGenerationBias'].GetValue)(wvl))
	# pylab.show()
	# kwargs['UnshadowedFraction'] = 5
	kwargs['DOMOversizeFactor'] = opts.oversize
	kwargs['DOMRadius'] = 356*I3Units.mm/2.
	
	# kwargs['OverrideApproximateNumberOfWorkItems'] = 1
	
	
	# kwargs['ExtraArgumentsToI3CLSimModule'] = dict(EnableDoubleBuffering=False)

if opts.unweighted:
	kwargs['UnweightedPhotons'] = True
	kwargs['DOMOversizeFactor'] = 1.
	kwargs['DOMRadius'] = 356*I3Units.mm/2.
	assert(opts.oversize == 1.)

import copy
class SizeModules(icetray.I3ConditionalModule):
	"""Restore OM sphere sizes to IceCube standard"""
	def __init__(self, ctx):
		super(SizeModules,self).__init__(ctx)
		self.AddParameter("Radius", "", 0.16510*icetray.I3Units.m)
	def Configure(self):
		self.radius = self.GetParameter("Radius")
		pass
	def Geometry(self, frame):
		# modgeomap = copy.copy(frame['I3ModuleGeoMap'])
		modgeomap = frame['I3ModuleGeoMap']
		# del frame['I3ModuleGeoMap']
		for k in modgeomap.keys():
			modgeo = modgeomap[k]
			modgeo.radius = self.radius
			modgeomap[k] = modgeo
		# frame['I3ModuleGeoMap'] = modgeomap
		self.PushFrame(frame)
# if opts.oversize != 1.:
# 	tray.Add(SizeModules, radius=0.16510*icetray.I3Units.m)

icetray.logging.set_level_for_unit('I3CLSimStepToPhotonConverterOpenCL', 'TRACE')
tray.Add(clsim.I3CLSimMakePhotons,
	ParallelEvents=1,
	RandomService=randomService,
	UseGPUs=opts.gpu,
	UseCPUs=not opts.gpu,
	# UseAllCPUCores=True,
	**kwargs
	)

tray.Add("Dump")

# if opts.oversize != 1.:
# 	tray.Add(SizeModules, radius=0.1778*icetray.I3Units.m)

# tray.Add(clsim.I3CLSimMakeHitsFromPhotons, DOMOversizeFactor=opts.oversize, UnshadowedFraction=0.99)

@icetray.traysegment
def I3CLSimMakeMDOMHitsFromPhotons(tray, name,
                               MCTreeName="I3MCTree_sliced",
                               PhotonSeriesName="PhotonSeriesMap",
                               MCPESeriesName="MCPESeriesMap",
                               RandomService=None,
                               DOMOversizeFactor=5.,
                               UnshadowedFraction=0.99,
                               UseHoleIceParameterization=True,
                               If=lambda f: True
                               ):
    """
    Convert I3Photons into I3MCPEs. This applies the DOM
    angular acceptance (and wavenelgth acceptance in case
    you are using the unbiased photon propagation mode.)

    :param MCTreeName:
        The name of the I3MCTree containing the particles to propagate.
    :param PhotonSeriesName:
        Name of the input I3PhotonSeriesMap to be converted.
    :param MCPESeriesName:
        Name of the output I3MCPESeriesMap written by the module.
        Set this to None to prevent generating MCPEs from
        Photons.
    :param RandomService:
        Set this to an instance of a I3RandomService. Alternatively,
        you can specify the name of a configured I3RandomServiceFactory
        added to I3Tray using tray.AddService(). If you don't configure
        this, the default I3RandomServiceFactory will be used.
    :param DOMOversizeFactor:
        Set the DOM oversize factor. To disable oversizing, set this to 1.
    :param UnshadowedFraction:
        Fraction of photocathode available to receive light (e.g. unshadowed by the cable)
    :param UseHoleIceParameterization:
        Use an angular acceptance correction for hole ice scattering.
    :param If:
        Python function to use as conditional execution test for segment modules.        
    """

    from icecube import icetray, dataclasses, clsim

    if False and DOMOversizeFactor != 1.:
        # photons came from standard IceCube simulation (13" DOM spheres)
        DOMRadius = 0.16510*icetray.I3Units.m
        # but we want to make hits on mDOM spheres (14" diameter)
        mDOMRadius = 0.1778*icetray.I3Units.m
    else:
        DOMRadius = 0.16510*icetray.I3Units.m
        mDOMRadius = DOMRadius
    
    # *(DOMRadius/mDOMRadius)
    print "oversize: %.2f pancacke: %.2f" % (DOMOversizeFactor*(DOMRadius/mDOMRadius), DOMOversizeFactor)
    pmtAcceptance = efficiencies.GetMDOMAcceptance(oversize=DOMOversizeFactor, efficiency=UnshadowedFraction)

    tray.AddModule("I3PhotonToMCHitConverterForMDOMs", name + "_clsim_make_hits",
                   RandomService = RandomService,
                   MCTreeName = MCTreeName,
                   InputPhotonSeriesMapName = PhotonSeriesName,
                   OutputMCPESeriesMapName = MCPESeriesName,
                   # Factor by which to squeeze detection spheroid transverse
                   # to photon direction
                   DOMOversizeFactor = DOMOversizeFactor*(DOMRadius/mDOMRadius),
                   # NB: since we reduced the oversize factor by the ratio of
                   # radii, the pancake factor is now larger than the oversize
                   # factor. This will have the effect of stretching the
                   # detection spheroid by the ratio of radii in the parallel
                   # direction
                   DOMPancakeFactor = DOMOversizeFactor,
                   
                   PMTWavelengthAcceptance = pmtAcceptance,
                   PMTAngularAcceptance = efficiencies.angularAcceptance,
                   GlassThickness = efficiencies.glassThickness,
                   GlassAbsorptionLength = efficiencies.glassAbsLen,
                   GelAbsorptionLength = efficiencies.gelAbsLen,
                   # IgnoreDOMsWithoutDetectorStatusEntry = False, # in icesim4 it is the job of the DOM simulation tools to cut out these DOMs
                   If=If)

# FIXME: this is not yet a thing
tray.Add(I3CLSimMakeMDOMHitsFromPhotons, DOMOversizeFactor=opts.oversize, UnshadowedFraction=0.99)

tray.Add("Delete", keys=["PhotonSeriesMap"])

tray.Add("I3NullSplitter", "nullsplit")

tray.Add("I3Writer", Streams=map(icetray.I3Frame.Stream, 'GCDQP'),
    filename=outfile)

tray.Execute()
tray.Finish()
