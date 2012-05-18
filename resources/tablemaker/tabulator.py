
from icecube.icetray import I3Units
from icecube.dataclasses import I3Position, I3Particle, I3Direction, I3Constants
from icecube.phys_services import I3Calculator
from icecube.clsim import I3Photon
import numpy, numpy_extensions, math

def potemkin_photon():
	ph = I3Photon()
	ph.numScattered = 1
	ph.groupVelocity = I3Constants.c/I3Constants.n_ice_group
	ph.startPos = I3Position(0,0,0)
	ph.startTime = 0.
	ph.pos = I3Position(100., 100., 100.,)
	ph.AppendToIntermediatePositionList(I3Position(100., 0., 0.,))
	
	return ph

def potemkin_source():
	p = I3Particle()
	p.pos = I3Position(0,0,0)
	p.dir = I3Direction(0,0,1)
	p.time = 0
	p.shape = p.Cascade
	p.energy = 1e3
	
	return p
	
class Tabulator(object):
	
	def __init__(self, source, nbins=(10, 10, 10, 10)):
		self.binedges = [
		    numpy.linspace(0, numpy.sqrt(600), nbins[0]+1)**2,
		    numpy.linspace(-1, 1, nbins[1]+1),
		    numpy.linspace(0, 180, nbins[2]+1),
		    numpy.linspace(0, numpy.sqrt(7e3), nbins[3]+1)**2,
		]
		self.values = numpy.zeros(tuple([len(v)-1 for v in self.binedges]), dtype=numpy.double)
		self.weights = numpy.zeros(self.values.shape)
		self.n_photons = 0
		self.source = source
		self._step_length = 1.*I3Units.m
		self._c = I3Constants.c/I3Constants.n_ice_group # FIXME: actually depends on local medium and wavelength!
		
	def _normalize(self):
		
		near = numpy.meshgrid_nd(*[b[:-1] for b in self.binedges[:-1]])
		far = numpy.meshgrid_nd(*[b[1:] for b in self.binedges[:-1]])
		
		# The effective area of the bin is its volume divided by
		# the sampling frequency
		area = ((far[0]**3-near[0]**3)/3.)*(far[1]-near[1])*(far[2]-near[2])/self._step_length
		fluxconst = area*self.n_photons
		
		# convert each weight (in units of photons * m^2) to a probability
		self.values /= fluxconst
		self.weights /= fluxconst*fluxconst
	
	def _getDirection(self, pos1, pos2):
		d = I3Direction(pos2.x-pos1.x, pos2.y-pos1.y, pos2.z-pos1.z)
		return numpy.array([d.x, d.y, d.z])
		
	def _getPerpendicularDirection(self, d):
		perpz = math.hypot(d.x, d.y)
		if perpz > 0:
			return numpy.array([-d.x*d.z/perpz, -d.y*d.z/perpz, perpz])
		else:
			return numpy.array([1., 0., 0.])

	def _getBinIndices(self, pos, dt):
		"""
		Sperical case only!
		
		TODO: implement detector-centered (degenerate) cylindrical coords
		"""
		p0 = numpy.array(self.source.pos)
		p1 = numpy.array(pos)
		
		source_dir = numpy.array([self.source.dir.x, self.source.dir.y, self.source.dir.z])
		perp_dir = self._getPerpendicularDirection(self.source.dir)
		
		norm = lambda arr: numpy.sqrt((arr**2).sum())
		
		displacement = p1-p0
		r = norm(displacement)
		l = numpy.dot(displacement, source_dir)
		rho = displacement - l*source_dir
		if norm(rho) > 0:
			azimuth = numpy.arccos(numpy.dot(rho, perp_dir)/norm(rho))/I3Units.degree
		else:
			azimuth = 0
		cosZen = math.hypot(l, norm(rho))
		idx = [numpy.searchsorted(edges, v) for v, edges in zip((r, cosZen, azimuth, dt), self.binedges)]
		for j, i in enumerate(idx):
			if i > len(self.binedges[j])-2:
				idx[j] = len(self.binedges[j])-2
		return tuple(idx)

	def _getTimeResidual(self, pos, t):
		return I3Calculator.time_residual(self.source, pos, t, I3Constants.n_ice_group, I3Constants.n_ice_phase)

	def recordPhoton(self, photon):
		positions = photon.positionList
		t = photon.startTime
		# TODO: get the actual absorption lengths for each intermediate step
		abs_lengths = numpy.arange(len(positions))
		for start, start_lengths, stop, stop_lengths in zip(positions[:-1], abs_lengths[:-1], positions[1:], abs_lengths[1:]):
			pdir = self._getDirection(start, stop)
			# XXX HACK: the cosine of the impact angle with the
			# DOM is the same as the z-component of the photon
			# direction if the DOM is pointed straight down.
			impact = pdir[2]
			distance = start.calc_distance(stop)
			# segment length in units of the local absorption length
			abs_distance = stop_lengths-start_lengths 
			
			# FIXME: this samples at fixed distances. Randomize instead?
			for d in numpy.arange(0, distance, self._step_length):
				# Which bin did we land in?
				pos = I3Position(*(numpy.array(start) + d*pdir))
				dt = self._getTimeResidual(pos, t + d/photon.groupVelocity)
				indices = self._getBinIndices(pos, dt)
			
				# TODO: add extra weights:
				#       - wavelength
				#       - DOM area
				#       - DOM angular sensitivity
				# weight will be multiplied by the DOM effective area, so units will be photons * m^2
				d_abs = start_lengths + (d/distance)*abs_distance
				weight = numpy.exp(-d_abs)
				self.values[indices] += weight
				self.weights[indices] += weight*weight
			
			t += distance/photon.groupVelocity
		self.n_photons += 1