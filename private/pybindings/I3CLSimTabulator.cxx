/**
 * Copyright (c) 2012
 * Jakob van Santen <jvansanten@icecube.wisc.edu>
 * and the IceCube Collaboration <http://www.icecube.wisc.edu>
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
 * SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
 * OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
 * CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 *
 * $Id$
 *
 * @file I3CLSimTabulator.cxx
 * @version $Revision$
 * @date $Date$
 * @author Jakob van Santen
 */

#include <phys-services/I3Calculator.h>
#include <dataclasses/I3Constants.h>
#include <boost/foreach.hpp>

#include <dataclasses/physics/I3Particle.h>
#include <dataclasses/I3Position.h>
#include <dataclasses/I3Direction.h>

#include <boost/python/tuple.hpp>
namespace bp = boost::python;

#ifdef USE_NUMPY
#define PY_ARRAY_UNIQUE_SYMBOL clsim_PyArray_API
#include <numpy/ndarrayobject.h>
#endif

#include "I3CLSimTabulator.h"

void
I3CLSimTabulator::SetBins(boost::python::object binEdges, double stepLength)
{
	stepLength_ = stepLength;
	
	Py_XDECREF(values_);
	Py_XDECREF(weights_);
	binEdges_.clear();
	
	for (int i = 0; i < bp::len(binEdges); i++)
		binEdges_.push_back(bp::extract<std::vector<double> >(binEdges[i]));	
	
	npy_intp dims[binEdges_.size()];
	for (unsigned i = 0; i < binEdges_.size(); i++) {
		const std::vector<double> &edges = binEdges_[i];
		if (!(edges.size() >= 2))
			log_fatal("Edge array in dimension %d has only %zu entries!", i, edges.size());
		
		dims[i] = edges.size()-1;
	}
	
	#ifdef USE_NUMPY
	values_ = PyArray_ZEROS(binEdges_.size(), dims, NPY_FLOAT, 0);
	weights_ = PyArray_ZEROS(binEdges_.size(), dims, NPY_FLOAT, 0);
	Py_INCREF(values_);
	Py_INCREF(weights_);
	#endif
}

I3CLSimTabulator::~I3CLSimTabulator()
{
	Py_XDECREF(values_);
	Py_XDECREF(weights_);
}

bp::object
I3CLSimTabulator::GetValues() const
{
#ifdef USE_NUMPY
	return bp::make_tuple(bp::object(bp::handle<>(values_)),
	    bp::object(bp::handle<>(weights_)));
#else
	return bp::object()
#endif	
}

void
I3CLSimTabulator::SetEfficiencies(I3CLSimFunctionConstPtr wavelengthAcceptance,
    I3CLSimFunctionConstPtr angularAcceptance, double domRadius)
{
	wavelengthAcceptance_ = wavelengthAcceptance;
	angularAcceptance_ = angularAcceptance;
	domArea_ = M_PI*domRadius*domRadius;
}

void
I3CLSimTabulator::SetRandomService(I3RandomServicePtr rng)
{
	rng_ = rng;
}

void
I3CLSimTabulator::Normalize()
{
	size_t tableSize = PyArray_SIZE(values_);
	
	// The collection efficiency of a bin is effectively the
	// number of chances a photon has to be collected within
	// that bin, given a fixed collection are and spatial
	// sampling frequency. We divide by this factor to turn
	// a sum of weights into a detection probability for a
	// DOM in this spatial bin.

	for (size_t i = 0; i < tableSize; i++) {
		double norm = GetBinVolume(i)/(stepLength_*domArea_);
		((float*)PyArray_DATA(values_))[i] /= norm;
		((float*)PyArray_DATA(weights_))[i] /= (norm*norm);
	}
}

double
I3CLSimTabulator::GetBinVolume(off_t idx) const
{
	// First, unravel the flattened index.
	off_t idxs[4];
	for (int i = 0; i < 4; i++) {
		idxs[i] = idx/(PyArray_STRIDE(values_, i)/sizeof(float)) % PyArray_DIM(values_, i);
		// printf("%d stride: %zd dim: %zd\n", i, PyArray_STRIDE(values_, i)/sizeof(float), PyArray_DIM(values_, i));
	}
	
	// NB: since we combine the bins at azimuth > 180 degrees with the other half of
	// the sphere, the true volume of an azimuthal bin is twice its nominal value.
	return ((std::pow(binEdges_[0][idxs[0]+1], 3) - std::pow(binEdges_[0][idxs[0]], 3))/3.)
	    * 2*I3Units::degree*(binEdges_[1][idxs[1]+1] - binEdges_[1][idxs[1]])
	    * (binEdges_[2][idxs[2]+1] - binEdges_[2][idxs[2]]);
}

off_t
I3CLSimTabulator::GetBinIndex(const I3CLSimTabulator::Source &source, const I3Position &pos, double t) const
{
	const I3Position displacement = pos-source.pos;
	double l = source.dir*displacement;
	const I3Position rho = displacement - l*source.dir;
	double n_rho = rho.Magnitude();
	
	double coords[4]; // {r, azimuth, cosZen, dt}
	coords[0] = displacement.Magnitude();
	coords[1] = (n_rho > 0) ?
	    std::acos((rho*source.perpdir)/n_rho)/I3Units::degree
	    : 0.;
	coords[2] = (coords[0] > 0) ? l/coords[0] : 1.;
	coords[3] = t - coords[0]*I3Constants::n_ice_group/I3Constants::c;
	
	// Bail if the photon is too delayed at this point to be recorded
	if (coords[3] > binEdges_[3].back())
		return -1;

        // Bail if any of the dimensions is nan
        for (int i=0; i < 4; i++) {
		if (!std::isfinite(coords[i])) {
			log_error("Coordinate %d is %f!", i, coords[i]);
			return -1;
		}
	}
	
	// Find the index of the appropriate bin for each dimension,
	// and compute an index into the flattened table array.
	off_t idx = 0;
	for (int i=0; i < 4; i++) {
		const std::vector<double> &edges = binEdges_[i];
		
		off_t dimidx;
		if (coords[i] <= edges.front())
			dimidx = 0;
		else if (coords[i] >= edges.back())
			dimidx = edges.size()-2;
		else
			dimidx = std::distance(edges.begin(),
			    std::lower_bound(edges.begin(), edges.end(), coords[i]))-1;
		
		assert(dimidx >= 0);
		assert(dimidx < PyArray_DIM(values_, i));
		
		idx += (PyArray_STRIDE(values_, i)/sizeof(float))*dimidx;
	}
	return idx;
}

I3CLSimTabulator::Source::Source(const I3Particle &p)
    : pos(p.GetPos()), time(p.GetTime()), dir(p.GetDir())
{
	// Get unit vectors pointing along the source direction
	// and perpendicular to it, towards +z. For vertical sources,
	// pick the +x direction.
	double perpz = hypot(dir.GetX(), dir.GetY());
	perpdir = (perpz > 0.) ?
	    I3Direction(-dir.GetX()*dir.GetZ()/perpz, -dir.GetY()*dir.GetZ()/perpz, perpz)
	    : I3Direction(1., 0., 0.);
}

void
I3CLSimTabulator::RecordPhotons(const I3Particle &source_p, const I3Map<ModuleKey, I3Vector<I3Photon> > &photon_map)
{
	Source source(source_p);
	typedef I3Map<ModuleKey, I3Vector<I3Photon> >::value_type photon_pair;
	BOOST_FOREACH(const photon_pair &photons, photon_map)
		BOOST_FOREACH(const I3Photon &photon, photons.second)
			RecordPhoton(source, photon);
}

void
I3CLSimTabulator::RecordPhoton(const I3CLSimTabulator::Source &source, const I3Photon &photon)
{
	double t = photon.GetStartTime();
	double wlenWeight =
	    wavelengthAcceptance_->GetValue(photon.GetWavelength())*photon.GetWeight();
	
	boost::optional<I3Position> p0 = photon.GetPositionListEntry(0);
	boost::optional<I3Position> p1;
	double absLengths[2] = {
		photon.GetDistanceInAbsorptionLengthsAtPositionListEntry(0),
		0.
	};
	uint32_t nsteps = photon.GetNumPositionListEntries();
	for (uint32_t i = 1; i < nsteps; i++, p0 = p1, absLengths[0] = absLengths[1]) {
		p1 = photon.GetPositionListEntry(i);
		absLengths[1] =
			photon.GetDistanceInAbsorptionLengthsAtPositionListEntry(i);
		if (!p0 || !p1)
			continue;
		
		// A vector connecting the two recording points.
<<<<<<< .working
		I3Position pdir = (*p1)-(*p0);
		double distance = pdir.Magnitude();
		pdir /= distance;
=======
		I3Position displacement(*p1-*p0);
		double distance = displacement.Magnitude();
		I3Direction pdir(displacement);
>>>>>>> .merge-right.r131588
		
		// XXX HACK: the cosine of the impact angle with the
		// DOM is the same as the z-component of the photon
		// direction if the DOM is pointed straight down.
		// This can be modified for detectors with other values
		// of pi.
		double impactWeight = wlenWeight*angularAcceptance_->GetValue(pdir.GetZ());
		
		int nsamples = floorf(distance/stepLength_);
		nsamples += (rng_->Uniform() < distance/stepLength_ - nsamples);
		for (int i = 0; i < nsamples; i++) {
			double d = distance*rng_->Uniform();
			off_t idx = GetBinIndex(source, *p0 + d*pdir, t + d/photon.GetGroupVelocity());
			// Once the photon has accumulated enough delay time
			// to run off the end of the table, there's no going back. Bail.
			if (idx < 0) {
				i = nsteps;
				break;
			}
			
			// Weight the photon by its probability of:
			// 1) Being detected, given its wavelength
			// 2) Being detected, given its impact angle with the DOM
			// 3) Having survived this far without being absorpbed
			double weight = impactWeight*std::exp(-(absLengths[0] +
			    (d/distance)*(absLengths[1]-absLengths[0])));
			
			((float*)PyArray_DATA(values_))[idx] += weight;
			((float*)PyArray_DATA(weights_))[idx] += weight*weight;
		}
		t += distance/photon.GetGroupVelocity();
	}
	assert( abs(t-photon.GetTime()) < 10 );
}

namespace bp = boost::python;

#ifdef USE_NUMPY
#if PY_MAJOR_VERSION >= 3
static PyObject *hack_import_array() {import_array(); return NULL;}
#else
static void hack_import_array() {import_array();}
#endif
#endif

void
register_I3CLSimTabulator()
{
#ifdef USE_NUMPY
	hack_import_array();
#endif
	bp::class_<I3CLSimTabulator, boost::shared_ptr<I3CLSimTabulator> >("I3CLSimTabulator")
	    .def("SetBins", &I3CLSimTabulator::SetBins)
	    .def("SetEfficiencies", &I3CLSimTabulator::SetEfficiencies)
	    .def("SetRandomService", &I3CLSimTabulator::SetRandomService)
	    .def("RecordPhotons", &I3CLSimTabulator::RecordPhotons)
	    .def("Normalize", &I3CLSimTabulator::Normalize)
	    .def("GetValues", &I3CLSimTabulator::GetValues)
	;	
}
