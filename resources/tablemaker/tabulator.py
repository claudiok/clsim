#!/usr/bin/env python

from icecube.icetray import I3Units, I3Module
from icecube.dataclasses import I3Position, I3Particle, I3Direction, I3Constants
from icecube.phys_services import I3Calculator
from icecube.clsim import I3Photon, GetIceCubeDOMAcceptance, GetIceCubeDOMAngularSensitivity
import numpy, math
from icecube.photospline import numpy_extensions # meshgrid_nd
	
class PhotoTable(object):
	
	def __init__(self, binedges, values, weights, n_photons=0):
		self.binedges = binedges
		self.values = values
		self.weights = weights
		self.bincenters = [(edges[:1]+edges[:-1])/2. for edges in self.binedges]
		self.n_photons = n_photons
		
	def __iadd__(self, other):
		if not isinstance(other, self.__class__):
			raise TypeError, "Can only add another %s" % (self.__class__.__name__)
		self.values += other.values
		self.weights += other.weights
		self.n_photons += other.n_photons
		
		return self
	
	def __idiv__(self, num):
		return self.__imul__(1./num)
		
	def __imul__(self, num):
		self.values *= num
		self.weights *= num*num
		
		return self
		
	def normalize(self):
		self /= self.n_photons
		
	def save(self, fname):
		import pyfits
		
		data = pyfits.PrimaryHDU(self.values)
		data.header.update('TYPE', 'Photon detection probability table')
		
		data.header.update('NPHOTONS', self.n_photons)
		
		hdulist = pyfits.HDUList([data])
		
		if self.weights is not None:
			errors = pyfits.ImageHDU(self.weights, name='ERRORS')
			hdulist.append(errors)
		
		for i in xrange(self.values.ndim):
			edgehdu = pyfits.ImageHDU(self.binedges[i],name='EDGES%d' % i)
			hdulist.append(edgehdu)
			
		hdulist.writeto(fname)
		
	@classmethod
	def load(cls, fname):
		import pyfits
		
		hdulist = pyfits.open(fname)
		data = hdulist[0]
		values = data.data
		binedges = []
		for i in xrange(values.ndim):
			binedges.append(hdulist['EDGES%d' % i].data)
		
		try:
			weights = hdulist['ERRORS'].data
		except KeyError:
			weights = None
			
		n_photons = data.header['NPHOTONS']
			
		return cls(binedges, values, weights, n_photons)
	
class Tabulator(object):
	
	def __init__(self, nbins=(10, 10, 10, 10)):
		self.binedges = [
		    numpy.linspace(0, numpy.sqrt(600), nbins[0]+1)**2,
		    numpy.linspace(-1, 1, nbins[1]+1),
		    numpy.linspace(0, 180, nbins[2]+1),
		    numpy.linspace(0, numpy.sqrt(7e3), nbins[3]+1)**2,
		]
		self.values = numpy.zeros(tuple([len(v)-1 for v in self.binedges]), dtype=numpy.double)
		self.weights = numpy.zeros(self.values.shape)
		self.n_photons = 0
		self._step_length = 1.*I3Units.m
		self._angularEfficiency = GetIceCubeDOMAngularSensitivity()
		self._effectiveArea = GetIceCubeDOMAcceptance()
		
	def save(self, fname):
		# normalize to a photon flux, but not the number of photons (for later stacking)
		area = self._getBinAreas()
		PhotoTable(self.binedges, self.values/area, self.weights/(area*area), self.n_photons).save(fname)
	
	def _getBinAreas(self):
		
		near = numpy.meshgrid_nd(*[b[:-1] for b in self.binedges[:-1]])
		far = numpy.meshgrid_nd(*[b[1:] for b in self.binedges[:-1]])
		
		# The effective area of the bin is its volume divided by
		# the sampling frequency
		area = ((far[0]**3-near[0]**3)/3.)*(far[1]-near[1])*(far[2]-near[2])/self._step_length
		return area
	
	def _normalize(self):
		
		area = self._getBinAreas()
		
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

	def _getBinIndices(self, source, pos, t):
		"""
		This is the only function that needs to be significantly modified
		to support different table geometries.
		
		NB: the currently implements sperical, source-centered coordinates only!
		
		TODO: implement detector-centered (degenerate) cylindrical coords
		"""
		p0 = numpy.array(source.pos)
		p1 = numpy.array(pos)
		
		dt = I3Calculator.time_residual(source, pos, t, I3Constants.n_ice_group, I3Constants.n_ice_phase)
		
		source_dir = numpy.array([source.dir.x, source.dir.y, source.dir.z])
		perp_dir = self._getPerpendicularDirection(source.dir)
		
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

	def recordPhoton(self, source, photon, n_photons=1):
		"""
		:param nphotons: number of photons to add to the total tally
		after recording. Set this to zero if you're going to do your own counting.
		"""
		# position of photon at each scatter
		positions = photon.positionList
		# path length in units of absorption length at each scatter
		abs_lengths = [photon.GetDistanceInAbsorptionLengthsAtPositionListEntry(i) for i in xrange(len(positions))]
		
		# The various weights are constant for different bits of the photon track.
		# Constant for a given photon:
		# - Photon weight: 1 / wavelength bias
		# - DOM effective area (as a function of wavelength)
		# Constant for paths between scatters:
		# - relative angular efficiency (as a function of impact angle)
		# Varies along the track:
		# - survival probability (exp(-path/absorption length))
		# The units of the combined weight are m^2 / photon
		wlen_weight = self._effectiveArea.GetValue(photon.wavelength)*photon.weight
		
		t = photon.startTime
		
		for start, start_lengths, stop, stop_lengths in zip(positions[:-1], abs_lengths[:-1], positions[1:], abs_lengths[1:]):
			pdir = self._getDirection(start, stop)
			# XXX HACK: the cosine of the impact angle with the
			# DOM is the same as the z-component of the photon
			# direction if the DOM is pointed straight down.
			# This can be modified for detectors with other values
			# of pi.
			impact = pdir[2]
			distance = start.calc_distance(stop)
			# segment length in units of the local absorption length
			abs_distance = stop_lengths-start_lengths 
			
			impact_weight = wlen_weight*self._angularEfficiency.GetValue(impact)
			
			# FIXME: this samples at fixed distances. Randomize instead?
			for d in numpy.arange(0, distance, self._step_length):
				# Which bin did we land in?
				pos = I3Position(*(numpy.array(start) + d*pdir))
				indices = self._getBinIndices(source, pos, t + d/photon.groupVelocity)
				
				d_abs = start_lengths + (d/distance)*abs_distance
				weight = impact_weight*numpy.exp(-d_abs)
				
				self.values[indices] += weight
				self.weights[indices] += weight*weight
			
			t += distance/photon.groupVelocity
		self.n_photons += n_photons
		
class I3TabulatorModule(I3Module, Tabulator):
	def __init__(self, ctx):
		I3Module.__init__(self, ctx)
		Tabulator.__init__(self)
		
		self.AddParameter("Filename", "Output filename", "foo.fits")
		self.AddParameter("Photons", "Name of I3PhotonSeriesMap in the frame", "I3PhotonSeriesMap")
		self.AddParameter("Statistics", "Name of the CLSimStatistics object in the frame", "")
		self.AddParameter("Source", "Name of the source I3Particle in the frame", "")
		self.AddOutBox("OutBox")
		
	def Configure(self):
		
		self.fname = self.GetParameter("Filename")
		self.photons = self.GetParameter("Photons")
		self.stats = self.GetParameter("Statistics")
		self.source = self.GetParameter("Source")
		
	def DAQ(self, frame):
		source = frame[self.source]
		photonmap = frame[self.photons]
		
		for photons in photonmap.itervalues():
			for photon in photons:
				self.recordPhoton(source, photon, 0)
		
		# Each I3Photon can only carry a fixed number of intermediate
		# steps, so there may be more than one I3Photon for each generated
		# photon.
		stats = frame[self.stats]
		self.n_photons += stats.GetTotalNumberOfPhotonsGenerated()
		
		self.PushFrame(frame)
		
	def Finish(self):
		
		self.save(self.fname)

def potemkin_photon():
	ph = I3Photon()
	ph.numScattered = 1
	ph.groupVelocity = I3Constants.c/I3Constants.n_ice_group
	ph.startPos = I3Position(0,0,0)
	ph.startTime = 0.
	ph.pos = I3Position(100., 100., 100.,)
	ph.AppendToIntermediatePositionList(I3Position(100., 0., 0.,), 1)
	ph.distanceInAbsorptionLengths = 2
	ph.wavelength = 450.
	ph.weight = 1.
	
	return ph

def potemkin_source():
	p = I3Particle()
	p.pos = I3Position(0,0,0)
	p.dir = I3Direction(0,0,1)
	p.time = 0
	p.shape = p.Cascade
	p.energy = 1e3
	
	return p

def test_tray():
	from I3Tray import I3Tray
	from icecube import dataio, phys_services
	from icecube.clsim import I3PhotonSeries, I3PhotonSeriesMap, I3CLSimEventStatistics
	from icecube.icetray import OMKey, I3Frame
	
	tray = I3Tray()
	
	tray.AddModule('I3InfiniteSource', 'maw')
	
	def inserty(frame):
		source = potemkin_source()
		psm = I3PhotonSeriesMap()
		psm[OMKey(0,0)] = I3PhotonSeries([potemkin_photon()])
		
		stats = I3CLSimEventStatistics()
		stats.AddNumPhotonsGeneratedWithWeights(numPhotons=1, weightsForPhotons=1, majorID=0, minorID=0)
		
		frame['Source'] = source
		frame['Photons'] = psm
		frame['Statistics'] = stats
		
	tray.AddModule(inserty, 'inserty', Streams=[I3Frame.DAQ])
	
	import os
	if os.path.exists('potemkin.fits'):
		os.unlink('potemkin.fits')
	tray.AddModule(I3TabulatorModule, 'tabbycat',
	    Photons='Photons', Source='Source', Statistics='Statistics',
	    Filename='potemkin.fits')
	
	tray.AddModule('TrashCan', 'YesWeCan')
	tray.Execute(1)
	tray.Finish()
	
if __name__ == "__main__":
	test_tray()
	 