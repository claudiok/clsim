/**
 * Copyright (c) 2015
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
 * @file Axes.cxx
 * @version $LastChangedRevision$
 * @date $Date$
 * @author Jakob van Santen
 */

#include "icetray/I3Units.h"
#include "clsim/tabulator/Axes.h"
#include "clsim/I3CLSimHelperToFloatString.h"
#include "opencl/I3CLSimHelperLoadProgramSource.h"

#include <boost/foreach.hpp>

namespace {

std::string 
loadKernel(const std::string& name, bool header=false)
{
    const std::string I3_BUILD(getenv("I3_BUILD"));
    const std::string kernelBaseDir = I3_BUILD+"/clsim/resources/kernels/";
    const std::string ext = header ? ".h.cl" : ".c.cl";
    return I3CLSimHelper::LoadProgramSource(kernelBaseDir+name+ext);
}

}

namespace clsim {

namespace tabulator {

Axes::Axes(const std::vector<value_type> &axes) : axes_(axes), n_dim_(axes_.size()),
    shape_(n_dim_), strides_(n_dim_)
{
	int i = n_dim_-1;
	// NB: every axis has an over- and an under-flow bin.
	shape_[i] = axes_[i]->GetNBins()+2;
	strides_[i] = 1;
	for (i--; i >= 0; i--) {
		shape_[i] = axes_[i]->GetNBins()+2;
		strides_[i] = strides_[i+1]*shape_[i+1];
	}
	n_bins_ = strides_[0]*shape_[0];
}

Axes::~Axes()
{}

std::string
Axes::GetBinIndexFunction() const
{
	std::ostringstream ss;
	ss << "inline uint getBinIndex(coordinate_t coords)";
	ss << "\n{\n";
	ss << "    return ";
	
	for (uint i = 0; i < axes_.size(); i++) {
		std::ostringstream var;
		var << "coords.s" << i;
		ss << strides_[i] << "*" << axes_[i]->GetIndexCode(var.str());
		if (i+1 < axes_.size())
			ss << "\n         + ";
	}
	
	ss << ";\n}\n";
	
	return ss.str();
}

std::string
Axes::GenerateBinningCode() const
{
	return GetCoordinateFunction() + "\n"
	    + GetBoundsCheckFunction() + "\n"
	    + GetBinIndexFunction() + "\n";
}

std::string
SphericalAxes::GetCoordinateFunction() const
{
	std::ostringstream ss;
	if (at(1)->GetMax() > 180.)
		ss << "#define HAS_FULL_AZIMUTH_EXTENSION\n";
	ss << loadKernel("spherical_coordinates");
	return ss.str();
}

std::string
SphericalAxes::GetBoundsCheckFunction() const
{
	std::ostringstream ss;
	ss << "inline bool isOutOfBounds(const coordinate_t coords)";
	ss << "\n{\n";
	ss << "    return (coords.s3 > "<<I3CLSimHelper::ToFloatString(this->at(3)->GetMax())<<")"
	    << "|| (coords.s0 > "<<I3CLSimHelper::ToFloatString(this->at(0)->GetMax())<<")"
	    << ";";
	ss << "\n}\n";

	return ss.str();
}

double
SphericalAxes::GetBinVolume(const std::vector<size_t> &idxs) const
{
	// NB: since we combine the bins at azimuth > 180 degrees with the
	// other half of the sphere, the true volume of an azimuthal bin is
	// twice its nominal value.
	// In case of a full table the azimuthal bin volume is simply the nominal volume.
	assert(idxs.size() >= 3);
	double scalefactor = (at(1)->GetMax() > 180.) ? 1 : 2;
	return ((std::pow(at(0)->GetBinEdge(idxs[0]+1), 3) - std::pow(at(0)->GetBinEdge(idxs[0]), 3))/3.)
	    * scalefactor*I3Units::degree*(at(1)->GetBinEdge(idxs[1]+1) - at(1)->GetBinEdge(idxs[1]))
	    * (at(2)->GetBinEdge(idxs[2]+1) - at(2)->GetBinEdge(idxs[2]));

}

std::string
CylindricalAxes::GetCoordinateFunction() const
{
	return loadKernel("cylindrical_coordinates");
}

std::string
CylindricalAxes::GetBoundsCheckFunction() const
{
	std::ostringstream ss;
	ss << "inline bool isOutOfBounds(const coordinate_t coords)";
	ss << "\n{\n";
	ss << "    return (coords.s3 > "<<I3CLSimHelper::ToFloatString(this->at(3)->GetMax())<<");";
	ss << "\n}\n";

	return ss.str();
}

double
CylindricalAxes::GetBinVolume(const std::vector<size_t> &idxs) const
{
	// NB: since we combine the bins at azimuth > pi with the
	// other half of the cylinder, the true volume of an azimuthal bin is
	// twice its nominal value.
	assert(idxs.size() >= 3);
	return ((std::pow(at(0)->GetBinEdge(idxs[0]+1), 2) - std::pow(at(0)->GetBinEdge(idxs[0]), 2))/2.)
	    * 2*(at(1)->GetBinEdge(idxs[1]+1) - at(1)->GetBinEdge(idxs[1]))
	    * (at(2)->GetBinEdge(idxs[2]+1) - at(2)->GetBinEdge(idxs[2]));
}

}

}
