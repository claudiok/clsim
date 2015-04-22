from icecube import icetray, dataclasses, dataio
from I3Tray import I3Tray
import numpy, pylab

import sys
infiles = sys.argv[1:]

tray = I3Tray()

tray.AddModule('I3Reader', 'reader', filenamelist=infiles)

def ignore_scan(frame):
	if frame['HealpixPixel'].value == 4141:
		return True
	else:
		return False
tray.AddModule(ignore_scan, "ignore_scan")

from icecube import photonics_service, millipede

def CLShim(tray, name):
	from icecube import icetray, dataclasses, clsim, phys_services
	import os

	from icecube.clsim.traysegments.common import configureOpenCLDevices, parseIceModel

	def getDetectorParameters(IceModel="spice_lea", DisableTilt=False, UnshadowedFraction=1.,
	    UseHoleIceParameterization=True, DOMOversizeFactor=5., UnWeightedPhotons=False):
    
	    DOMRadius = 0.16510*icetray.I3Units.m # 13" diameter
    
	    # ice properties
	    if isinstance(IceModel, str):
	        mediumProperties = parseIceModel(os.path.expandvars('$I3_BUILD/clsim/resources/ice/%s' % IceModel), disableTilt=DisableTilt)
	    else:
	        # get ice model directly if not a string
	        mediumProperties = IceModel

	    # detector properties
	    if UseHoleIceParameterization:
	        # the hole ice acceptance curve peaks at 0.75 instead of 1
	        domEfficiencyCorrection = UnshadowedFraction*0.75*1.35 * 1.01 # DeepCore DOMs have a relative efficiency of 1.35 plus security margin of +1%
	    else:
	        domEfficiencyCorrection = UnshadowedFraction*1.35      * 1.01 # security margin of +1%
	    domAcceptance = clsim.GetIceCubeDOMAcceptance(domRadius = DOMRadius*DOMOversizeFactor, efficiency=domEfficiencyCorrection)
	    domAngularSensitivity = clsim.GetIceCubeDOMAngularSensitivity(holeIce=UseHoleIceParameterization)

	    if not UnWeightedPhotons:
	        wavelengthGenerationBias = domAcceptance
	    else:
	        wavelengthGenerationBias = None

	    wavelengthGenerators = clsim.I3CLSimRandomValuePtrSeries()
	    wavelengthGenerators.append(clsim.makeCherenkovWavelengthGenerator(wavelengthGenerationBias, False, mediumProperties))
	    return mediumProperties, wavelengthGenerationBias, wavelengthGenerators, domAcceptance, domAngularSensitivity, DOMOversizeFactor

	def makeStepConverter(UseGeant4=False):
		#if UseGeant4:
		clsim.AutoSetGeant4Environment()
		converter = clsim.I3CLSimLightSourceToStepConverterGeant4()
		ppcConverter = clsim.I3CLSimLightSourceToStepConverterPPC(photonsPerStep=200)
		parameterizationList = clsim.GetDefaultParameterizationList(ppcConverter, muonOnly=False)
		converter.SetLightSourceParameterizationSeries(parameterizationList)
		
		return converter

	def makeGeometry(positions, DOMOversizeFactor=5., spacing=30):
		import math
		string, dom = 0, 1
		lastx, lasty = float('nan'), float('nan')
	
		DOMRadius = 0.16510*icetray.I3Units.m # 13" diameter
	
		omGeos = dataclasses.I3OMGeoMap()
		moduleGeos = dataclasses.I3ModuleGeoMap()
		subdetectors = dataclasses.I3MapModuleKeyString()
		offsets = dict()
		for i, pos in enumerate(positions):
			if not math.hypot(pos.x-lastx, pos.y-lasty) < spacing:
				lastx, lasty = pos.x, pos.y
				string += 1
				dom = 1
			key = icetray.OMKey(string, dom)
			mkey = dataclasses.ModuleKey(string, dom)
			dom += 1
			omGeo = dataclasses.I3OMGeo()
			moduleGeo = dataclasses.I3ModuleGeo()
			omGeo.position = pos
			omGeo.omtype = omGeo.IceCube
			moduleGeo.pos = pos
			moduleGeo.radius = DOMRadius
			moduleGeo.module_type = moduleGeo.IceCube
			omGeos[key] = omGeo
			moduleGeos[mkey] = moduleGeo
			subdetectors[mkey] = "IceCube"
			offsets[i] = key
		frame = icetray.I3Frame()
		frame['I3OMGeoMap'] = omGeos
		frame['I3ModuleGeoMap'] = moduleGeos
		frame['Subdetectors'] = subdetectors

		simplegeo = clsim.I3CLSimSimpleGeometryFromI3Geometry(DOMRadius, DOMOversizeFactor, frame)
	
		return simplegeo

	def makeKernel(device, rng, detector_params_args=dict()):
		"""
		Configure propagator and step generator, but do not initialize
		"""
		mediumProperties, wavelengthGenerationBias, wavelengthGenerators, domAcceptance, domAngularSensitivity, DOMOversizeFactor = getDetectorParameters(**detector_params_args)
	
		propagator = clsim.I3CLSimStepToPhotonConverterOpenCL(rng)
		propagator.SetDevice(device)
		propagator.SetMediumProperties(mediumProperties)
		propagator.SetWlenBias(wavelengthGenerationBias)
		propagator.SetWlenGenerators(wavelengthGenerators)
		propagator.SetDOMPancakeFactor(DOMOversizeFactor)
		propagator.SetMaxNumWorkitems(768*1400)
		
		stepGenerator = makeStepConverter()
		stepGenerator.SetMediumProperties(mediumProperties)
		stepGenerator.SetRandomService(rng)
		stepGenerator.SetWlenBias(domAcceptance)
		
		return stepGenerator, propagator, domAcceptance, domAngularSensitivity
	
	device = configureOpenCLDevices(UseGPUs=True, UseCPUs=False)[0]
	# device = configureOpenCLDevices(UseGPUs=False, UseCPUs=True, DoNotParallelize=False)[0]
	rng = phys_services.I3GSLRandomService(1337)
	stepGenerator, propagator, domAcceptance, domAngularSensitivity = makeKernel(device, rng)
	
	tray.AddService('I3CLShimFactory', name, AngularAcceptance=domAngularSensitivity,
	    WavelengthAcceptance=domAcceptance, StepConverter=stepGenerator, PhotonConverter=propagator,
	    GeometryFactory=makeGeometry, OversampleFactor=1000, CacheDepth=5)

class MilliPlot(icetray.I3ConditionalModule):
	def __init__(self, context):
		super(MilliPlot, self).__init__(context)
		self.AddOutBox("OutBox")
		self.AddParameter("Seeds", "", [])
		self.AddParameter("NDOMs", "", None)
		self.AddParameter("Differential", "", True)
		
		self.millipede = millipede.PyPyMillipede(context)
		config = self.millipede.GetConfiguration()
		for k in config.keys():
			self.AddParameter(k, config.descriptions[k], config[k])
	
	def Configure(self):
		self.seeds = self.GetParameter("Seeds")
		self.differential = self.GetParameter("Differential")
		self.ndoms = self.GetParameter("NDOMs")
		config = self.millipede.GetConfiguration()
		for k in config.keys():
			self.millipede.SetParameter(k, self.GetParameter(k))
	
	@staticmethod
	def errorbar(axes, domCache, differential=False, **kwargs):
		"""
		Plot the data in this MillipedeDOMCache with counting error bars.
		"""
		from scipy.stats import chi2
		
		left = domCache.time_bin_edges[:-1]
		right = domCache.time_bin_edges[1:]
		widths = right-left
		# construct a 68% credible interval for the mean of a Poisson
		# distribution, given a single sample
		cl = 0.68
		a = (1-cl)/2.
		lo = chi2.ppf(a, 2*domCache.charges)/2.
		lo[numpy.isnan(lo)] = 0
		hi = chi2.ppf(1-a, 2*(domCache.charges+1))/2.
		yerr = domCache.charges-lo, hi-domCache.charges
		y = domCache.charges
		if differential:
			yerr = (yerr[0]/widths, yerr[1]/widths)
			y = y/widths
		#fmt=None, ecolor=color, label='OM(%d,%d): %.1f PE' % (k.string, k.om, dc.charges.sum())
		return axes.errorbar((left+right)/2, y, xerr=(right-left)/2., yerr=yerr, **kwargs)
		
	@staticmethod
	def steps(axes, domCache, expectations, differential=False, **kwargs):
		x = numpy.zeros( 2*len(domCache.time_bin_edges), numpy.float )
		y = numpy.zeros( 2*len(domCache.time_bin_edges), numpy.float )
		
		if differential:
			ex = expectations/numpy.diff(domCache.time_bin_edges)
		else:
			ex = expectations
		
		x[0::2], x[1::2] = domCache.time_bin_edges, domCache.time_bin_edges
		y[1:-1:2], y[2::2] = ex, ex
		return axes.plot(x, y, **kwargs)
	
	def Physics(self, frame):
		from scipy.special import gammaln
		
		self.millipede.DatamapFromFrame(frame)
		
		# Map keys to slices of the response matrix
		slices = dict()
		i = 0
		for k, dc in self.millipede.domCache.iteritems():
			valid = sum(dc.valid)
			slices[k] = slice(i, i+valid)
			i += valid
		
		hypotheses = dict()
		for name in self.seeds:
			source = dataclasses.I3Particle(frame[name])
			if source.type == source.unknown:
				source.type = source.EMinus
			sources = dataclasses.I3VectorI3Particle([source])
			response = self.millipede.GetResponseMatrix(sources)
			response = numpy.asarray(response.to_I3Matrix())
			expectations = numpy.inner(response, [p.energy for p in sources])
			hypotheses[name] = (source, expectations)
		
		# Get 10 highest-charge DOMs
		keys = sorted(self.millipede.domCache.keys(), key=lambda k: self.millipede.domCache[k].charges.sum())[::-1][slice(None,self.ndoms)]
		
		for k in keys:
			dc = self.millipede.domCache[k]
			
			if dc.charges.sum() < 50:
				continue
			
			#, label='OM(%d,%d): %.1f PE' % (k.string, k.om, dc.charges.sum()
			fig = pylab.figure()
			axes = pylab.gca()
			colors = axes._get_lines.color_cycle
			for name in self.seeds:
				source, expectations = hypotheses[name]
				
				ex = numpy.zeros(dc.charges.size)
				ex[numpy.where(dc.valid)[0]] = expectations[slices[k]]
				
				mask = numpy.where(numpy.asarray(dc.valid)&(ex>0))[0]
				q = dc.charges[mask]
				lam = ex[mask]
				logl = q*numpy.log(lam) - lam - gammaln(q+1)
				
				color = colors.next()
				self.steps(axes, dc, ex, self.differential, color=color, linewidth=2, label='%s: $\ln L=%.1f$' % (name, logl.sum()))
			
			if self.differential:
				axes.set_ylabel('Photocathode current [PE/ns]')
			else:
				axes.set_ylabel('Collected charge [PE]')
			axes.set_xlabel('Time [ns]')
			axes.legend(loc='best', prop=dict(size='small'))
			
			self.errorbar(axes, dc, self.differential, fmt='ko', markersize=5, ecolor='k')
			
			start = dc.time_bin_edges[1]-50
			stop = min(start + 1000, dc.time_bin_edges[-3]+50)
			axes.set_xlim(start, stop)
			ylim = axes.get_ylim()
			axes.set_ylim(max(0, ylim[0]), ylim[1])
			
			axes.set_title('OM(%d,%d): %.1f PE' % (k.string, k.om, dc.charges.sum()))
			#fig.savefig('OM%.2d%.2d.pdf' % (k.string, k.om))
		
		pylab.show()
		
		self.PushFrame(frame)

if False:
	base = "/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/ems_mie_z20_a10.%s.fits"
	pxs = photonics_service.I3PhotoSplineService(base % "abs", base % "prob", 0)
	from os.path import basename
	# pxs = photonics_service.I3LoggingPhotonicsService(pxs, basename(base % "abs")+".log", basename(base % "prob")+".log")
else:
	pxs = "CLShim"
	tray.AddSegment(CLShim, pxs)
icetray.logging.I3Logger.global_logger.set_level_for_unit('I3PhotoSplineService', icetray.logging.I3LogLevel.LOG_TRACE);
Pulses='SplitInIcePulses'
tray.AddModule('Delete', keys=['SaturatedDOMs'])

def cleanup(frame):
	del frame["DeepCoreDOMs"]
tray.AddModule(cleanup, "cleanup", streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics, icetray.I3Frame.Geometry])
ExcludedDOMs = tray.AddSegment(millipede.HighEnergyExclusions, 'Excludey', Pulses=Pulses, ExcludeBrightDOMs=False)

Seeds = ['MillipedeStarting2ndPass'] #['MonopodFit']#, 'MonopodWithoutBrights', 'HESE_Cscd_Credo20']
tray.AddModule(MilliPlot, Pulses=Pulses, Seeds=Seeds, Differential=True,
    CascadePhotonicsService=pxs, ExcludedDOMs=ExcludedDOMs, PhotonsPerBin=25, NDOMs=300,
)

tray.AddModule('TrashCan', 'YesWeCan')
tray.Execute()
tray.Finish()
