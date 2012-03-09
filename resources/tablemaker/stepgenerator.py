from icecube import icetray, dataclasses, phys_services, clsim
import os
from I3Tray import I3Units

def setup_converter(useGeant4=False):
	# make a converter
	if useGeant4:
		ppcConverter = clsim.I3CLSimLightSourceToStepConverterGeant4()
	else:
		ppcConverter = clsim.I3CLSimLightSourceToStepConverterPPC(photonsPerStep=200)

	# initialize it
	randomGen = phys_services.I3SPRNGRandomService(
		seed = 123456,
		nstreams = 10000,
		streamnum = 1)
	mediumProperties = clsim.MakeIceCubeMediumProperties()

	#DOMRadius = 0.16510*icetray.I3Units.m # 13" diameter
	#RadiusOverSizeFactor = 5.
	#domAcceptance = clsim.GetIceCubeDOMAcceptance(domRadius = DOMRadius*RadiusOverSizeFactor)
	domAcceptance = clsim.I3CLSimWlenDependentValueConstant(1.)

	# lets set it up
	ppcConverter.SetMediumProperties(mediumProperties)
	ppcConverter.SetRandomService(randomGen)
	ppcConverter.SetWlenBias(domAcceptance)

	ppcConverter.SetMaxBunchSize(10240)
	ppcConverter.SetBunchSizeGranularity(1)

	ppcConverter.Initialize()

	return ppcConverter

def gen_steps(particle, converter, copies=1):
	
	# insert the requested number of copies of the particle
	for i in range(copies):
		# make a light source from the particle
		lightSource = clsim.I3CLSimLightSource(particle)
		
		# put it in the queue
		converter.EnqueueLightSource(lightSource, i)

	# tell the converter that we want all the results now
	converter.EnqueueBarrier()
	
	# retrieve all results
	steps = clsim.I3CLSimStepSeries()
	while True:
		#barrierReset=False
		#stepSeries = converter.GetConversionResultWithBarrierInfo(barrierReset)
		
		stepSeries = converter.GetConversionResult()
		barrierReset = not converter.BarrierActive()
		
		steps.extend(stepSeries)		
		
		# get out of the loop if the barrier has been reset
		if barrierReset:
			break 
	
	return steps
	
def make_source(energy):
	p = dataclasses.I3Particle()
	p.pos = dataclasses.I3Position(0.,0.,0.)
	p.dir = dataclasses.I3Direction(0.,0.,1.)
	p.time = 0.
	p.energy = energy
	p.shape = dataclasses.I3Particle.ParticleShape.Cascade
	p.type = dataclasses.I3Particle.ParticleType.EMinus
	p.length = 0.
	p.location_type = dataclasses.I3Particle.LocationType.InIce
	
	return p
	
def make_steps(energy, copies=1):
	
	ppcConverter = setup_converter()
	
	p = make_source(energy)
	steps = gen_steps(p, ppcConverter, copies=copies)
	return steps
	
def dump_steps(fname, energy, copies):
	import cPickle as pickle
	steps = make_steps(energy, copies)
	f = open('foo.pickle', 'w')
	pickle.dump(steps, f, 2)
	f.close()
	
def read_steps(fname):
	import cPickle as pickle
	import numpy
	f = open(fname, 'r')
	steps = pickle.load(f)
	f.close()
	dtype = zip(('x', 'y', 'z', 'theta', 'phi', 'num'), (float,)*10)
	arr = numpy.empty(len(steps), dtype=dtype)
	for i, step in enumerate(steps):
		arr[i] = (step.pos.x, step.pos.y, step.pos.z, step.theta, step.phi, step.num)
	return arr
	
def hobo_geometry():
	import numpy
	
	from itertools import izip, product, cycle
	DOMRadius = 0.16510*icetray.I3Units.m # 13" diameter
	nx, ny, nz = 5, 5, 5
	geometry = clsim.I3CLSimSimpleGeometryUserConfigurable(DOMRadius, nx*ny*nz)
	gridx = numpy.linspace(-50, 50, nx)
	gridy = numpy.linspace(-50, 50, ny)
	gridz = numpy.linspace(50, -50, nz)
	grid = map(lambda f: enumerate(f), (gridx, gridy, gridz))
	
	for idx, ((i,x),(j,y),(k,z)) in enumerate(product(*grid)):
		geometry.SetStringID(idx, i*nx+j)
		geometry.SetDomID(idx, k)
		geometry.SetPosX(idx, x)
		geometry.SetPosY(idx, y)
		geometry.SetPosZ(idx, z)
		
	return geometry
	
class HoboTracker(clsim.I3CLSimStepToPhotonConverterOpenCL):
	def __init__(self, rng, nativemath):
		clsim.I3CLSimStepToPhotonConverterOpenCL.__init__(self, rng, nativemath)
	def _GetGeometrySource(self):
		return str()
	def GetCollisionDetectionSource(self, get_header):
		import numpy
		def format_array(arr, name):
			content = ',\n'.join(['%.6ff' % f for f in arr])
			return "#define %s_SIZE %d\n__constant floating_t %s[%s_SIZE] = {%s};\n" % (name, arr.size, name, name, content)
		def center(arr):
			return 0.5*(arr[1:]+arr[:-1])
		if get_header:
			return open('spherical_collision_kernel.h.cl').read()
		else:
			source = str()
			r_nbins = 50
			rbins = numpy.linspace(0, 25., 51)
			rbincenters = (0.5*(rbins[1:]+rbins[:-1]))**2
			ctbins = numpy.linspace(-1., 1., 7)
			
			source += format_array(center(numpy.linspace(0, 25., 51))**2, "RBINS");
			source += format_array(center(numpy.linspace(-1., 1., 7)), "COSTHETABINS");
			
			source = open('spherical_collision_kernel.c.cl').read()		    
			
			return source
		    
		   
		
		
	
def configure_devices(geometry, UseGPUs=False, UseCPUs=True, debug=False):
	
	rng = phys_services.I3GSLRandomService(1337)
	
        # it's a PPC ice description directory
	IceModelLocation = os.path.expandvars("$I3_SRC/clsim/resources/ice/spice_mie")
        mediumProperties = clsim.MakeIceCubeMediumProperties(iceDataDirectory=IceModelLocation)
	
	# Set up detector parameters
	DOMRadius = 0.16510*icetray.I3Units.m # 13" diameter
	RadiusOverSizeFactor = 5.
	# detector properties
        domAcceptance = clsim.GetIceCubeDOMAcceptance(domRadius = DOMRadius*RadiusOverSizeFactor)
        domAngularSensitivity = clsim.GetIceCubeDOMAngularSensitivity(holeIce=True)
	
	wavelengthGenerationBias = domAcceptance
	
	wavelengthGenerators = clsim.I3CLSimRandomValuePtrSeries()
	wavelengthGenerators.append(clsim.makeCherenkovWavelengthGenerator(wavelengthGenerationBias, False, mediumProperties))
	
	# get OpenCL devices
	openCLDevices = [device for device in clsim.I3CLSimOpenCLDevice.GetAllDevices() if (device.gpu and UseGPUs) or (device.cpu and UseCPUs)]

	converters = []

        # (auto-)configure OpenCL devices
	for device in openCLDevices:
		if 'Tesla' in device.device or 'GTX' in device.device:
		# assume these are "fast", all others are "slow"
			device.useNativeMath=True
			device.approximateNumberOfWorkItems=1024000
		else:
			device.useNativeMath=False
			device.approximateNumberOfWorkItems=10240
		if debug:
			device.approximateNumberOfWorkItems=1
		# mc = clsim.I3CLSimStepToPhotonConverterOpenCL(rng, device.GetUseNativeMath())
		mc = HoboTracker(rng, device.GetUseNativeMath())
		
		mc.SetDevice(device)
		
	        mc.SetWlenGenerators(wavelengthGenerators)
	        mc.SetWlenBias(wavelengthGenerationBias)

	        mc.SetMediumProperties(mediumProperties);
	        mc.SetGeometry(geometry);

	        mc.SetEnableDoubleBuffering(True);
	        mc.SetDoublePrecision(False);
	        mc.SetStopDetectedPhotons(False);

	        try:
			mc.Compile();
		finally:
			f = open('thebad.c', 'w')
			f.write(mc.GetFullSource())
			f.close()
        
	        maxWorkgroupSize = mc.GetMaxWorkgroupSize();
	        mc.SetWorkgroupSize(maxWorkgroupSize);

	        workgroupSize = mc.GetWorkgroupSize();
		
	        # use approximately the given number of work items, convert to a multiple of the workgroup size
		maxNumWorkitems = (device.approximateNumberOfWorkItems/workgroupSize)*workgroupSize
	        if maxNumWorkitems==0:
			maxNumWorkitems=workgroupSize;
        
	        mc.SetMaxNumWorkitems(maxNumWorkitems);

	        mc.Initialize();
		
		converters.append(mc)
	return converters
	
def photomc(steps, geometry, UseGPUs=True, debug=False):
	import itertools, operator
	mc = configure_devices(geometry, UseCPUs=not UseGPUs, UseGPUs=UseGPUs, debug=debug)[0]
	
	norm = sum(itertools.imap(operator.attrgetter("num"), steps))
	
	DOMRadius = 0.16510*icetray.I3Units.m # 13" diameter
	RadiusOverSizeFactor = 5.
        domAcceptance = clsim.GetIceCubeDOMAcceptance(domRadius = DOMRadius*RadiusOverSizeFactor)
        domAngularSensitivity = clsim.GetIceCubeDOMAngularSensitivity(holeIce=True)
	
	# Normalize the photon weight
	def weight_photon(photon):
		photon.weight *= domAcceptance.GetValue(photon.wavelength)*domAngularSensitivity.GetValue(-photon.dir.z)
		photon.weight /= norm
	
	blocksize = mc.maxNumWorkitems
	nblocks = len(steps)/blocksize
	for i in xrange(nblocks):
		mc.EnqueueSteps(steps[i*blocksize:(i+1)*blocksize], i)
	
	photons = list()
	for i in xrange(nblocks):
		# wait for results to come back
		try:
			result = mc.GetConversionResult()
			for ph in result[1]:
				weight_photon(ph)
				photons.append(ph)
		except IndexError, e:
			print e
			continue			

	return photons

def single_step(num):
	step = clsim.I3CLSimStep()
	step.dir = dataclasses.I3Direction(0,0)
	step.pos = dataclasses.I3Position(0,0,0)
	step.time = 0
	step.weight = 1
	step.beta = 1
	step.sourceType = dataclasses.I3Particle.EMinus
	step.length = 0
	step.num = int(num)
	return [step]
	
if __name__ == "__main__":
	#steps = make_steps(2., 1)
	steps = single_step(num=1e2)
	geometry = hobo_geometry()
	
	photons = photomc(steps, geometry, UseGPUs=False, debug=True)
	print len(photons), 'photons'
		
		
	
	
