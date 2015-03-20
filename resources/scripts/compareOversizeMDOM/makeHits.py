#!/usr/bin/env python

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-s", "--seed", type=int, default=1)
parser.add_option("--oversize", type=float, default=1)
parser.add_option("--unweighted", default=False, action="store_true")
opts, args = parser.parse_args()

infiles, outfile = args[:-1], args[-1]

from icecube import icetray, dataclasses, dataio, phys_services, sim_services
from icecube.icetray import I3Units
from I3Tray import I3Tray

tray = I3Tray()

randomService = phys_services.I3GSLRandomService(opts.seed)
tray.context['I3RandomService'] = randomService

tray.Add("I3Reader", filenamelist=infiles)

from icecube import clsim
import efficiencies

kwargs = dict()

if opts.oversize != 1.:
	kwargs['WavelengthGenerationBias'] = efficiencies.GetAcceptanceEnvelope(oversize=opts.oversize)
	kwargs['DOMOversizeFactor'] = opts.oversize

if opts.unweighted:
	kwargs['UnweightedPhotons'] = True
	kwargs['DOMOversizeFactor'] = 1.
	kwargs['DOMRadius'] = 356*I3Units.mm/2.
	assert(opts.oversize == 1.)

tray.Add("Dump")


icetray.logging.set_level_for_unit('I3CLSimStepToPhotonConverterOpenCL', 'TRACE')
tray.Add(clsim.I3CLSimMakePhotons,
	ParallelEvents=1,
	RandomService=randomService,
	UseGPUs=True,
	UseCPUs=False,
	**kwargs
	)


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

    if DOMOversizeFactor != 1.:
        # photons came from standard IceCube simulation (13" DOM spheres)
        DOMRadius = 0.16510*icetray.I3Units.m
        # but we want to make hits on mDOM spheres (14" diameter)
        mDOMRadius = 0.1778*icetray.I3Units.m
    else:
        DOMRadius = 0.16510*icetray.I3Units.m
        mDOMRadius = DOMRadius
    
    pmtAcceptance = efficiencies.GetMDOMAcceptance(oversize=DOMOversizeFactor*(DOMRadius/mDOMRadius), efficiency=UnshadowedFraction)

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

tray.Add("I3Writer", Streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics],
    filename=outfile)

tray.Execute()
tray.Finish()
