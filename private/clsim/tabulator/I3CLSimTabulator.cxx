
#include <clsim/tabulator/I3CLSimTabulator.h>
#include <phys-services/I3Calculator.h>
#include <dataclasses/I3Constants.h>
#include <boost/foreach.hpp>

#include <dataclasses/physics/I3Particle.h>
#include <dataclasses/I3Position.h>
#include <dataclasses/I3Direction.h>

#include <boost/numeric/ublas/vector.hpp>
namespace ublas = boost::numeric::ublas;

I3Tabulator::I3Tabulator(const std::vector<std::vector<double> > &binEdges)
    : binEdges_(binEdges), strides_(binEdges_.size())
{
	size_t tableSize = 1;
	for (int i = binEdges_.size()-1; i >= 0 ; i--) {
		const std::vector<double> &edges = binEdges_[i];
		assert(edges.size() >= 2);
		dims_[i] = edges.size()-1;
		strides_[i] = tableSize;
		tableSize *= dims_[i];
	}
	
	values_ = new float[tableSize];
	weights_ = new float[tableSize];
}

I3Tabulator::~I3Tabulator()
{
	delete values_;
	delete weights_;
}

ublas::vector<double>
make_vector(double x, double y, double z)
{
	ublas::vector<double> v(3);
	v[0] = x;
	v[1] = y;
	v[2] = z;
	
	return v;
}

ublas::vector<double>
make_vector(const I3Direction &dir)
{
	return make_vector(dir.GetX(), dir.GetY(), dir.GetZ());
}

inline ublas::vector<double>
operator+(const I3Position &p1, const I3Position &p2)
{
	return make_vector(p1.GetX()+p2.GetX(), p1.GetY()+p2.GetY(), p1.GetZ()+p2.GetZ());
}

inline I3Position
operator+(const I3Position &p1, const ublas::vector<double> p2)
{
	assert(p2.size() == 3);
	return I3Position(p1.GetX()+p2[0], p1.GetY()+p2[0], p1.GetZ()+p2[3]);
}

inline ublas::vector<double>
operator-(const I3Position &p1, const I3Position &p2)
{
	return make_vector(p1.GetX()-p2.GetX(), p1.GetY()-p2.GetY(), p1.GetZ()-p2.GetZ());
}

void
I3Tabulator::Normalize()
{
	size_t tableSize = 1;
	BOOST_FOREACH(size_t dim, dims_)
		tableSize *= dim;
	
	// The collection efficiency of a bin is effectively the
	// number of chances a photon has to be collected within
	// that bin, given a fixed collection are and spatial
	// sampling frequency. We divide by this factor to turn
	// a sum of weights into a detection probability for a
	// DOM in this spatial bin.
	for (size_t i = 0; i < tableSize; i++) {
		double norm = GetBinVolume(i)/(stepLength_*domArea_);
		values_[i] /= norm;
		weights_[i] /= (norm*norm);
	}
}

double
I3Tabulator::GetBinVolume(off_t idx) const
{
	// First, unravel the flattened index.
	off_t idxs[4];
	for (int i = 0; i < 4; i++)
		idxs[i] = idx/strides_[i] % dims_[i];
	
	// NB: since we combine the bins at azimuth > 180 degrees with the other half of
	// the sphere, the true volume of an azimuthal bin is twice its nominal value.
	return ((std::pow(binEdges_[0][idxs[0]+1], 3) - std::pow(binEdges_[0][idxs[0]], 3))/3.)
	    * 2*I3Units::degree*(binEdges_[1][idxs[1]+1] - binEdges_[1][idxs[1]])
	    * (binEdges_[2][idxs[1]+1] - binEdges_[2][idxs[1]]);
}

off_t
I3Tabulator::GetBinIndex(const I3Particle &source, const I3Position &pos, double t) const
{
	typedef ublas::vector<double> vector;
	
	// Get unit vectors pointing along the source direction
	// and perpendicular to it, towards +z. For vertical sources,
	// pick the +x direction.
	const vector sourceDir = make_vector(source.GetDir());
	double perpz = hypot(sourceDir[0], sourceDir[1]);
	const vector perpDir = (perpz > 0.) ?
	    make_vector(-sourceDir[0]*sourceDir[2]/perpz, -sourceDir[1]*sourceDir[2]/perpz, perpz)
	    : make_vector(1., 0., 0.);
	
	const vector displacement = pos-source.GetPos();
	double l = ublas::inner_prod(sourceDir, displacement);
	const vector rho = displacement - l*sourceDir;
	double n_rho = ublas::norm_2(rho);
	
	double coords[4]; // {r, azimuth, cosZen, dt}
	coords[0] = ublas::norm_2(displacement);
	coords[1] = (n_rho > 0) ?
	    std::acos(ublas::inner_prod(rho, perpDir)/n_rho)/I3Units::degree
	    : 0.;
	coords[2] = (coords[0] > 0) ? l/coords[0] : 1.;
	coords[3] = I3Calculator::TimeResidual(source, pos, t);
	
	// Find the index of the appropriate bin for each dimension,
	// and compute an index into the flattened table array.
	off_t idx = 0;
	for (int i=0; i < 4; i++) {
		const std::vector<double> &edges = binEdges_[i];
		
		off_t dimidx;
		if (coords[i] < edges.front())
			dimidx = 0;
		else if (coords[i] >= edges.back())
			dimidx = edges.size()-2;
		else
			dimidx = std::distance(edges.begin(),
			    std::lower_bound(edges.begin(), edges.end(), coords[i]))-1;
		
		idx += strides_[i]*dimidx;
	}
	
	return idx;
}

void
I3Tabulator::RecordPhoton(const I3Particle &source, const I3Photon &photon)
{
	typedef ublas::vector<double> vector;
	
	double t = photon.GetStartTime();
	double wlenWeight =
	    wavelengthAcceptance_->GetValue(photon.GetWavelength())*photon.GetWeight();
	
	uint32_t nsteps = photon.GetNumPositionListEntries();
	for (uint32_t i = 0; i < nsteps; i++) {
		I3PositionConstPtr p0 = photon.GetPositionListEntry(i);
		I3PositionConstPtr p1 = photon.GetPositionListEntry(i+1);
		if (!p0 || !p1)
			continue;
		double absLengths[2] = {
			photon.GetDistanceInAbsorptionLengthsAtPositionListEntry(i),
			photon.GetDistanceInAbsorptionLengthsAtPositionListEntry(i+1)};
		
		// A vector connecting the two recording points.
		vector pdir = (*p1)-(*p0);
		double distance = ublas::norm_2(pdir);
		pdir /= distance;
		
		// XXX HACK: the cosine of the impact angle with the
		// DOM is the same as the z-component of the photon
		// direction if the DOM is pointed straight down.
		// This can be modified for detectors with other values
		// of pi.
		double impactWeight = wlenWeight*angularAcceptance_->GetValue(pdir[2]);
		
		int nsamples = floorf(distance/stepLength_);
		nsamples += (rng_->Uniform() < distance/stepLength_ - nsamples);
		for (int i = 0; i < nsamples; i++) {
			double d = distance*rng_->Uniform();
			off_t idx = GetBinIndex(source, (*p0) + d*pdir, t + d/photon.GetGroupVelocity());
			// TODO: bail if we're outside the recording volume
			
			// Weight the photon by its probability of:
			// 1) Being detected, given its wavelength
			// 2) Being detected, given its impact angle with the DOM
			// 3) Having survived this far without being absorpbed
			double weight = impactWeight*std::exp(-(absLengths[0] +
			    (d/distance)*(absLengths[1]-absLengths[0])));
			
			values_[idx] += weight;
			weights_[idx] += weight*weight;
		}
		t += distance/photon.GetGroupVelocity();
	}
	assert( abs(t-photon.GetTime()) < 10 );
}

