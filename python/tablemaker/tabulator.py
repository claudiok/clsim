# Copyright (c) 2012
# Jakob van Santen <jvansanten@icecube.wisc.edu>
# Claudio Kopper <claudio.kopper@icecube.wisc.edu>
# and the IceCube Collaboration <http://www.icecube.wisc.edu>
#
# Permission to use, copy, modify, and/or distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
# SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
# OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
# CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#
#
# $Id$
#
# @file tabulator.py
# @version $Revision$
# @date $Date$
# @author Jakob van Santen

from __future__ import print_function

from icecube.icetray import I3Units, I3Module, traysegment
from icecube.dataclasses import I3Position, I3Particle, I3MCTree, I3Direction, I3Constants
from icecube.phys_services import I3Calculator, I3GSLRandomService
from icecube.clsim import I3Photon, I3CLSimTabulator, GetIceCubeDOMAcceptance, GetIceCubeDOMAngularSensitivity
from icecube.clsim import FlasherInfoVectToFlasherPulseSeriesConverter, I3CLSimFlasherPulse, I3CLSimFlasherPulseSeries
from icecube.clsim.traysegments.common import parseIceModel
from icecube.clsim.util import GetRefractiveIndexRange
import numpy, math
from icecube.photospline import numpy_extensions # meshgrid_nd
from icecube.photospline.photonics import Table, Efficiency, Geometry, Parity

# A default header. This contains the same keys as the one created in photo2numpy from photospline.
empty_header = {
    'n_photons':         0,
    'efficiency':        Efficiency.NONE,
    'geometry':          Geometry.SPHERICAL,
    'parity':            Parity.EVEN,
    'zenith':            0.,
    'z':                 0.,
    'energy':            0.,
    'type':              int(I3Particle.ParticleType.unknown),
    'level':             1,
    'n_group':           I3Constants.n_ice_group,
    'n_phase':           I3Constants.n_ice_phase,
}

class PhotoTable(Table, object):
    """
    A (mostly) drop-in replacement for the Photonics table wrapper from photospline.
    """
    
    def __init__(self, binedges, values, weights, header=empty_header):
        self.bin_edges = binedges
        self.values = values
        self.weights = weights
        self.bin_centers = [(edges[1:]+edges[:-1])/2. for edges in self.bin_edges]
        self.bin_widths = [numpy.diff(edges) for edges in self.bin_edges]
        self.header = header
        
    def __iadd__(self, other):
        self.raise_if_incompatible(other)
        self.values += other.values
        self.weights += other.weights
        self.header['n_photons'] += other.header['n_photons']
        
        return self

    def __idiv__(self, num):
        return self.__imul__(1./num)
        
    def __imul__(self, num):
        self.values *= num
        self.weights *= num*num
        
        return self
        
    def raise_if_incompatible(self, other):
        """
        Check for generalized brokenness.
        """
        if not isinstance(other, self.__class__):
            raise TypeError("Can't combine a %s with this %s" % (other.__class__.__name__, self.__class__.__name__))
        if self.values.shape != other.values.shape:
            raise ValueError("Shape mismatch in data arrays!")
        nans = self.values.size - numpy.isfinite(self.values).sum()
        if nans != 0:
            raise ValueError("This table has %d NaN values. You might want to see to that.")
        nans = other.values.size - numpy.isfinite(other.values).sum()
        if nans != 0:
            raise ValueError("Other table has %d NaN values. You might want to see to that.")
        for k, v in self.header.items():
            if k == 'n_photons':
                continue
            if other.header[k] != v:
                raise ValueError("Can't combine tables with %s=%s and %s" % (k, v, other.header[k]))
        
    def normalize(self):
        if not self.header['efficiency'] & Efficiency.N_PHOTON:
            self /= self.header['n_photons']
            self.header['efficiency'] |= Efficiency.N_PHOTON
        
    def save(self, fname, overwrite=False):
        import pyfits, os
        
        if os.path.exists(fname):
            if overwrite:
                os.unlink(fname)
            else:
                raise IOError("File '%s' exists!" % fname)
        
        data = pyfits.PrimaryHDU(self.values)
        data.header.update('TYPE', 'Photon detection probability table')
        
        for k, v in self.header.items():
            # work around 8-char limit in FITS keywords
            tag = 'hierarch _i3_' + k
            data.header.update(tag, v)
        
        hdulist = pyfits.HDUList([data])
        
        if self.weights is not None:
            errors = pyfits.ImageHDU(self.weights, name='ERRORS')
            hdulist.append(errors)
        
        for i in range(self.values.ndim):
            edgehdu = pyfits.ImageHDU(self.bin_edges[i],name='EDGES%d' % i)
            hdulist.append(edgehdu)
            
        hdulist.writeto(fname)
        
    @classmethod
    def load(cls, fname):
        import pyfits
        
        hdulist = pyfits.open(fname)
        data = hdulist[0]
        values = data.data
        binedges = []
        for i in range(values.ndim):
            binedges.append(hdulist['EDGES%d' % i].data)
        
        try:
            weights = hdulist['ERRORS'].data
        except KeyError:
            weights = None
        
        header = dict()
        for k in map(str.lower, data.header.keys()):
            if k.startswith('_i3_'):
                header[k[4:]] = data.header[k]
            
        return cls(binedges, values, weights, header)
    
class Tabulator(object):
    _dtype = numpy.float32
    
    def SetBins(self, binedges, step_length=1.*I3Units.m):
        self.binedges = binedges
        self.values = numpy.zeros(tuple([len(v)-1 for v in self.binedges]), dtype=self._dtype)
        self.weights = numpy.zeros(self.values.shape, dtype=self._dtype)
        
        self._step_length = step_length
        
    def SetEfficiencies(self, wavelengthAcceptance=GetIceCubeDOMAcceptance(), angularAcceptance=GetIceCubeDOMAngularSensitivity(), domRadius=0.16510*I3Units.m):
        self._angularEfficiency = angularAcceptance
        domArea = numpy.pi*domRadius**2
        
        self._effectiveArea = lambda wvl: domArea*wavelengthAcceptance.GetValue(wvl)
    
    def SetRandomService(self, rng):
        self.rng = rng
    
    def GetValues(self):
        return self.values, self.weights
        
    def save(self, fname):
        # normalize to a photon flux, but not the number of photons (for later stacking)
        area = self._getBinAreas()
        print("Total weights:", self.n_photons)
        PhotoTable(self.binedges, self.values/area, self.weights/(area*area), self.n_photons).save(fname)
    
    def _getBinAreas(self):
        
        near = numpy.meshgrid_nd(*[b[:-1].astype(self._dtype) for b in self.binedges[:-1]], lex_order=True)
        far = numpy.meshgrid_nd(*[b[1:].astype(self._dtype) for b in self.binedges[:-1]], lex_order=True)
        
        # The effective area of the bin is its volume divided by
        # the sampling frequency. 
        # NB: since we're condensing onto a half-sphere, the azimuthal extent of each bin is doubled.
        area = ((far[0]**3-near[0]**3)/3.)*(2*(far[1]-near[1])*I3Units.degree)*(far[2]-near[2])
        print('Total volume:', area.sum())
        area = area.reshape(area.shape + (1,))/self._dtype(self._step_length)
        return area
        
    @staticmethod
    def GetBinVolumes(binedges, dtype=numpy.float32):
        
        near = numpy.meshgrid_nd(*[b[:-1].astype(dtype) for b in binedges[:-1]], lex_order=True)
        far = numpy.meshgrid_nd(*[b[1:].astype(dtype) for b in binedges[:-1]], lex_order=True)
        
        # The effective area of the bin is its volume divided by
        # the sampling frequency. 
        # NB: since we're condensing onto a half-sphere, the azimuthal extent of each bin is doubled.
        volumes = ((far[0]**3-near[0]**3)/3.)*(2*(far[1]-near[1])*I3Units.degree)*(far[2]-near[2])
        volumes = volumes.reshape(volumes.shape + (1,))
        return volumes
    
    def Normalize(self):
        
        fluxconst = self._getBinAreas()/self._step_length
        
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
        if r > 0:
            cosZen = l/r
        else:
            cosZen = 1
        idx = [numpy.searchsorted(edges, v) for v, edges in zip((r, azimuth, cosZen, dt), self.binedges)]
        for j, i in enumerate(idx):
            if i > len(self.binedges[j])-2:
                idx[j] = len(self.binedges[j])-2
        return tuple(idx)

    def RecordPhoton(self, source, photon):
        """
        :param nphotons: number of photons to add to the total tally
        after recording. Set this to zero if you're going to do your own counting.
        """
        # position of photon at each scatter
        positions = photon.positionList
        # path length in units of absorption length at each scatter
        abs_lengths = [photon.GetDistanceInAbsorptionLengthsAtPositionListEntry(i) for i in range(len(positions))]
        
        # The various weights are constant for different bits of the photon track.
        # Constant for a given photon:
        # - Photon weight: 1 / wavelength bias
        # - DOM effective area (as a function of wavelength)
        # Constant for paths between scatters:
        # - relative angular efficiency (as a function of impact angle)
        # Varies along the track:
        # - survival probability (exp(-path/absorption length))
        # The units of the combined weight are m^2 / photon
        wlen_weight = self._effectiveArea(photon.wavelength)*photon.weight
        
        t = photon.startTime
        
        for start, start_lengths, stop, stop_lengths in zip(positions[:-1], abs_lengths[:-1], positions[1:], abs_lengths[1:]):
            if start is None or stop is None:
                continue
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
            nsamples = int(numpy.floor(distance/self._step_length))
            nsamples += int(self.rng.uniform(0,1) < (distance/self._step_length - nsamples))
            for i in range(nsamples):
                d = distance*self.rng.uniform(0,1)
            # for d in numpy.arange(0, distance, self._step_length):
                # Which bin did we land in?
                pos = I3Position(*(numpy.array(start) + d*pdir))
                indices = self._getBinIndices(source, pos, t + d/photon.groupVelocity)
                # We've run off the edge of the recording volume. Bail.
                if indices[0] == len(self.binedges[0]) or indices[3] == len(self.binedges[3]):
                    break
                d_abs = start_lengths + (d/distance)*abs_distance
                weight = impact_weight*numpy.exp(-d_abs)
                
                self.values[indices] += weight
                self.weights[indices] += weight*weight
            
            t += distance/photon.groupVelocity

# tabulator = Tabulator
tabulator = I3CLSimTabulator

class I3TabulatorModule(I3Module, tabulator):
    def __init__(self, ctx):
        I3Module.__init__(self, ctx)
        tabulator.__init__(self)
        
        self.AddParameter("Filename", "Output filename", "foo.fits")
        self.AddParameter("Photons", "Name of I3PhotonSeriesMap in the frame", "I3PhotonSeriesMap")
        self.AddParameter("Statistics", "Name of the CLSimStatistics object in the frame", "")
        self.AddParameter("Source", "Name of the source I3Particle in the frame", "")
        self.AddParameter("RandomService", "A random number service", None)
        self.AddParameter("StepLength", "The mean step size for volume sampling", 1*I3Units.m)
        self.AddParameter("TableHeader", "A dictionary of source depth, orientation, etc", empty_header)
        
        nbins=(200, 36, 100, 105)
        self.binedges = [
            numpy.linspace(0, numpy.sqrt(580), nbins[0]+1)**2,
            numpy.linspace(0, 180, nbins[1]+1),
            numpy.linspace(-1, 1, nbins[2]+1),
            numpy.linspace(0, numpy.sqrt(7e3), nbins[3]+1)**2,
        ]
        self.AddParameter("BinEdges", "A list of the bin edges in each dimension", self.binedges)
        
        self.AddOutBox("OutBox")
        
    def Configure(self):
        
        self.fname = self.GetParameter("Filename")
        self.photons = self.GetParameter("Photons")
        self.stats = self.GetParameter("Statistics")
        self.source = self.GetParameter("Source")
        self.rng = self.GetParameter("RandomService")
        self.stepLength = self.GetParameter("StepLength")
        self.header = self.GetParameter("TableHeader")
        
        self.n_photons = 0
        
        self.binedges = self.GetParameter("BinEdges")
        
        self.domRadius = 0.16510*I3Units.m
        self.stepLength = 1*I3Units.m
        
        self.SetBins(self.binedges, self.stepLength)
        # XXX magic scaling factor used in CLSim weighted-photon generation
        self.SetEfficiencies(GetIceCubeDOMAcceptance(efficiency=0.75*1.35 * 1.01), GetIceCubeDOMAngularSensitivity(), self.domRadius)
        self.SetRandomService(self.rng)
        
    def DAQ(self, frame):
        source = frame[self.source]
        photonmap = frame[self.photons]
        
        for photons in photonmap.values():
            for photon in photons:
                self.RecordPhoton(source, photon)
        
        # Each I3Photon can only carry a fixed number of intermediate
        # steps, so there may be more than one I3Photon for each generated
        # photon. Since the generation spectrum is biased, we want to
        # normalize to the sum of weights, not the number of photons.
        stats = frame[self.stats]
        print('%f photons at DOMs (total weight %f), %f generated (total weight %f)' % (
            stats.GetTotalNumberOfPhotonsAtDOMs(), stats.GetTotalSumOfWeightsPhotonsAtDOMs(),
            stats.GetTotalNumberOfPhotonsGenerated(), stats.GetTotalSumOfWeightsPhotonsGenerated()))
        self.n_photons += stats.GetTotalSumOfWeightsPhotonsAtDOMs()
        
        self.PushFrame(frame)
        
    def Finish(self):
        self.Normalize()
        values, weights = self.GetValues()
        self.header['n_photons'] = self.n_photons
        table = PhotoTable(self.binedges, values, weights, self.header)
        table.save(self.fname)

def generate_seed():
    import struct
    with open('/dev/random') as rand:
        return struct.unpack('i', rand.read(4))[0]

class MakeParticle(I3Module):
    def __init__(self, context):
        I3Module.__init__(self, context)
        self.AddParameter("I3RandomService", "the service", None)
        self.AddParameter("PhotonSource", "", "CASCADE")
        self.AddParameter("ParticleType", "", I3Particle.ParticleType.EMinus)
        self.AddParameter("Energy", "", 1.*I3Units.GeV)
        self.AddParameter("Zenith", "", 90.*I3Units.degree)
        self.AddParameter("Azimuth", "", 0.*I3Units.degree)
        self.AddParameter("ZCoordinate", "", 0.*I3Units.m)
        self.AddParameter("FlasherWidth", "", 127)
        self.AddParameter("FlasherBrightness", "", 127)
        self.AddParameter("NEvents", "", 100)

        self.AddOutBox("OutBox")        

    def Configure(self):
        self.rs = self.GetParameter("I3RandomService")
        self.photonSource = self.GetParameter("PhotonSource")
        self.particleType = self.GetParameter("ParticleType")
        self.energy = self.GetParameter("Energy")
        self.zenith = self.GetParameter("Zenith")
        self.azimuth = self.GetParameter("Azimuth")
        self.zCoordinate = self.GetParameter("ZCoordinate")
        self.flasherWidth = self.GetParameter("FlasherWidth")
        self.flasherBrightness = self.GetParameter("FlasherBrightness")
        self.nEvents = self.GetParameter("NEvents")
        
        self.emittedEvents=0

    def DAQ(self, frame):

        # create source particle
        source = I3Particle()
        source.type = self.particleType
        source.energy = self.energy
        source.pos = I3Position(0., 0., self.zCoordinate)
        source.dir = I3Direction(self.zenith, self.azimuth)
        source.time = 0.
        source.location_type = I3Particle.LocationType.InIce

        # the table-making module prefers plain I3Particles
        frame["Source"] = source

        if self.photonSource.upper() == "CASCADE":

            primary = I3Particle()
            primary.type = I3Particle.ParticleType.NuE
            primary.energy = self.energy
            primary.pos = I3Position(0., 0., self.zCoordinate)
            primary.dir = I3Direction(self.zenith, self.azimuth)
            primary.time = 0.
            primary.location_type = I3Particle.LocationType.Anywhere

            mctree = I3MCTree()
            mctree.add_primary(primary)
            mctree.append_child(primary, source)
    
            # clsim likes I3MCTrees
            frame["I3MCTree"] = mctree
            
        elif self.photonSource.upper() == "FLASHER":

            outputSeries = I3CLSimFlasherPulseSeries()
            flasherPulse = self.getFlasherPulse(0., 0., self.zCoordinate, self.zenith, self.azimuth)
            outputSeries.append(flasherPulse)

            # clsim likes I3CLSimFlasherPulseSeries
            frame["I3FlasherPulseSeriesMap"] =  outputSeries

        self.emittedEvents += 1
        
        self.PushFrame(frame)
        
        print(self.emittedEvents)
        if self.emittedEvents >= self.nEvents:
            self.RequestSuspension()
    
    def getFlasherPulse(self, x, y, z, zenith, azimuth):

        pulse = I3CLSimFlasherPulse()
        pulse.type = I3CLSimFlasherPulse.FlasherPulseType.LED405nm
        pulse.pos = I3Position(x, y, z)
        pulse.dir = I3Direction(zenith, azimuth)
        pulse.time = 0.
        pulse.pulseWidth = (float(self.flasherWidth)/2.)*I3Units.ns
        scale = 32582*5.21 # scale down to match 1 GeV equivalent electromagnetic cascade energy
        pulse.numberOfPhotonsNoBias = 1.17e10/scale*(0.0006753+0.00005593*float(self.flasherBrightness))*(float(self.flasherWidth)+13.9-(57.5/(1.+float(self.flasherBrightness)/34.4)))

        if numpy.abs(zenith - 90.*I3Units.degree) > 22.5*I3Units.degree:
            tiltedFlasher = True # this is only a rough approximation to describe a tilted flasher
        else:
            tiltedFlasher = False

        pulse.angularEmissionSigmaPolar = FlasherInfoVectToFlasherPulseSeriesConverter.LEDangularEmissionProfile[(pulse.type, tiltedFlasher)][0]
        pulse.angularEmissionSigmaAzimuthal = FlasherInfoVectToFlasherPulseSeriesConverter.LEDangularEmissionProfile[(pulse.type, tiltedFlasher)][1]

        return pulse

def get_minimum_refractive_index(IceModelLocation, DisableTilt):
    mediumProperties = parseIceModel(IceModelLocation, DisableTilt)
    n_group_min, n_group_max = GetRefractiveIndexRange.GetGroupRefractiveIndexRange(mediumProperties)
    n_phase_min, n_phase_max = GetRefractiveIndexRange.GetPhaseRefractiveIndexRange(mediumProperties)
    return (n_group_min, n_phase_min)

@traysegment
def PhotonGenerator(tray, name, PhotonSource="CASCADE", Zenith=90.*I3Units.degree, ZCoordinate=0.*I3Units.m,
    Energy=1.*I3Units.GeV, FlasherWidth=127, FlasherBrightness=127, Seed=12345, NEvents=100,
    IceModel='spice_mie', DisableTilt=False):
    
    """

    :param PhotonSource: the type of photon source (CASCADE or FLASHER)
    :param Zenith: the orientation of the source
    :param ZCoordinate: the depth of the source
    :param Energy: the energy of the source (only for cascade tables)
    :param FlasherWidth: the width of the flasher pulse (only for flasher tables)
    :param FlasherBrightness: the brightness of the flasher pulse (only for flasher tables)
    :param Seed: the seed for the random number service
    :param NEvents: the number of events to simulate
    :param IceModel: the path to an ice model in $I3_SRC/clsim/resources/ice. Likely values include:
        'spice_mie' ppc-style SPICE-Mie parametrization
        'photonics_spice_1/Ice_table.spice.i3coords.cos080.10feb2010.txt' Photonics-style SPICE1 table
        'photonics_wham/Ice_table.wham.i3coords.cos094.11jul2011.txt' Photonics-style WHAM! table
    :param DisableTilt: if true, disable tilt in ice model
    
    """

    # check sanity of args
    if PhotonSource.upper() != "CASCADE" and PhotonSource.upper() != "FLASHER":
        print("photon source %s is unknown. Please specify either CASCADE or FLASHER!" % PhotonSource)
        sys.exit(1)
    
    from icecube import icetray, dataclasses, dataio, phys_services, sim_services, clsim
    from os.path import expandvars
    
    # a random number generator
    randomService = phys_services.I3SPRNGRandomService(
        seed = Seed,
        nstreams = 10000,
        streamnum = 0)
        
    tray.AddModule("I3InfiniteSource",name+"streams",
                   Stream=icetray.I3Frame.DAQ)

    tray.AddModule("I3MCEventHeaderGenerator",name+"gen_header",
                   Year=2009,
                   DAQTime=158100000000000000,
                   RunNumber=1,
                   EventID=1,
                   IncrementEventID=True)

    # empty GCD
    tray.AddService("I3EmptyStreamsFactory", name+"empty_streams",
        InstallCalibration=True,
        InstallGeometry=True,
        InstallStatus=True)
    tray.AddModule("I3MetaSynth", name+"synth")

    ptype = I3Particle.ParticleType.EMinus

    tray.AddModule(MakeParticle, name+"MakeParticle", PhotonSource=PhotonSource, Zenith=Zenith, ZCoordinate=ZCoordinate, 
        ParticleType=ptype, Energy=Energy, FlasherWidth=FlasherWidth, FlasherBrightness=FlasherBrightness, NEvents=NEvents)

    if PhotonSource.upper() == "FLASHER":
        flasherpulse = "I3FlasherPulseSeriesMap"
        mctree = None
    elif PhotonSource.upper() == "CASCADE":
        flasherpulse = None
        mctree = "I3MCTree"
                
    tray.AddSegment(clsim.I3CLSimMakePhotons, name+"makeCLSimPhotons",
        PhotonSeriesName = "PropagatedPhotons",
        MCTreeName = mctree,                        # if source is a cascade this will point to the I3MCTree
        FlasherPulseSeriesName = flasherpulse,      # if source is a flasher this will point to the I3CLSimFlasherPulseSeries
        MMCTrackListName = None,                    # do NOT use MMC
        ParallelEvents = 1,                         # only work at one event at a time (it'll take long enough)
        RandomService = randomService,
        # UnWeightedPhotons=True,
        UseGPUs=False,                              # table-making is not a workload particularly suited to GPUs
        UseCPUs=True,                               # it should work fine on CPUs, though
        UnshadowedFraction=1.0,                     # no cable shadow
        DOMOversizeFactor=1.0,                      # no oversizing (there are no DOMs, so this is pointless anyway)
        StopDetectedPhotons=False,                  # do not stop photons on detection (also somewhat pointless without DOMs)
        PhotonHistoryEntries=10000,                 # record all photon paths
        DoNotParallelize=True,                      # no multithreading
        OverrideApproximateNumberOfWorkItems=1,     # if you *would* use multi-threading, this would be the maximum number of jobs to run in parallel (OpenCL is free to split them)
        ExtraArgumentsToI3CLSimModule=dict(SaveAllPhotons=True,                 # save all photons, regardless of them hitting anything
                                           SaveAllPhotonsPrescale=1.,           # do not prescale the generated photons
                                           StatisticsName="I3CLSimStatistics",  # save a statistics object (contains the initial number of photons)
                                           FixedNumberOfAbsorptionLengths=46.,  # this is approx. the number used by photonics (it uses -ln(1e-20))
                                           LimitWorkgroupSize=1,                # this effectively disables all parallelism (there should be only one OpenCL worker thread)
                                                                                #  it also should save LOTS of memory
                                          ),
        IceModelLocation=expandvars("$I3_SRC/clsim/resources/ice/" + IceModel),
        DisableTilt=DisableTilt,
    )
    
    n_group, n_phase = get_minimum_refractive_index(expandvars("$I3_SRC/clsim/resources/ice/" + IceModel), DisableTilt)
    
    header = dict(empty_header)
    header['zenith'] = Zenith/I3Units.degree
    header['z'] = ZCoordinate
    header['energy'] = Energy
    header['type'] = int(ptype)
    header['n_group'] = n_group
    header['n_phase'] = n_phase
    header['efficiency'] = Efficiency.RECEIVER | Efficiency.WAVELENGTH
    if PhotonSource.upper() == "FLASHER":
        header['flasherwidth'] = FlasherWidth
        header['flasherbrightness'] = FlasherBrightness
        
    return randomService, header

