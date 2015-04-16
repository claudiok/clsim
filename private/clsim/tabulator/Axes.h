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
 * @file Axes.h
 * @version $LastChangedRevision$
 * @date $Date$
 * @author Jakob van Santen
 */

#ifndef CLSIM_TABULATOR_AXES_H_INCLUDED
#define CLSIM_TABULATOR_AXES_H_INCLUDED

#include "icetray/I3PointerTypedefs.h"
#include "clsim/tabulator/Axis.h"

#include <vector>
#include <string>

#include <boost/noncopyable.hpp>

namespace clsim {

namespace tabulator {

/// The Axes encapsulates the coordinate system used for binning photon paths.
class Axes : boost::noncopyable {
public:
	typedef boost::shared_ptr<Axis> value_type;
	
	Axes(const std::vector<value_type> &);
	virtual ~Axes();
	
	const value_type& at(unsigned i) const { return axes_.at(i); };
	
	size_t GetNDim() const { return n_dim_; };
	/// Total number of bins in the bin content array
	size_t GetNBins() const { return n_bins_; };
	/// Number of bins in each dimension
	std::vector<size_t> GetShape() const { return shape_; };
	/// Offset between adjacent bins in each dimension
	std::vector<size_t> GetStrides() const { return strides_; };

	/// Generate the OpenCL code required to transform photon positions
	/// into source-relative positions and calculate the corresponding index
	/// in the flattened bin-content array.
	std::string GenerateBinningCode() const;

	/// Calculate the volume (in m^3) of the cell at index multiIndex in
	/// the bin-content array
	virtual double GetBinVolume(const std::vector<size_t> &multiIndex) const = 0;

protected:
	/// Generate an OpenCL function that calculates the source-relative
	/// coordinates from the photon position and time
	virtual std::string GetCoordinateFunction() const = 0;
	/// Generate an OpenCL function that returns true if the photon has
	/// exited the recording volume and should be stopped
	virtual std::string GetBoundsCheckFunction() const = 0;
	

private:
	/// Generate an OpenCL function that transforms a set of coordinates
	/// into an index in the flattened bin content array. This is
	/// synthesized from the properties of the individual Axis instances.
	std::string GetBinIndexFunction() const;
	
	std::vector<value_type> axes_;
	size_t n_dim_;
	size_t n_bins_;
	std::vector<size_t> shape_;
	std::vector<size_t> strides_;
};

I3_POINTER_TYPEDEFS(Axes);

/// [Half]-Spherical coordinate system centered on the source position,
/// appropriate for approximately point-like sources
class SphericalAxes : public Axes {
public:
	SphericalAxes(const std::vector<value_type> &axes) : Axes(axes) {}
	
	virtual double GetBinVolume(const std::vector<size_t> &multiIndex) const;
	
protected:
	virtual std::string GetCoordinateFunction() const;
	virtual std::string GetBoundsCheckFunction() const;
	
};

/// [Half]-Cylindrical coordinate system centered on the source axis,
/// appropriate for a source moving at the speed of light with infinite range.
/// NB: since the position of the source is degenerate with time, the reciever
///     depth is used as a coordinate instead of the source depth.
class CylindricalAxes : public Axes {
public:
	CylindricalAxes(const std::vector<value_type> &axes) : Axes(axes) {}
	
	virtual double GetBinVolume(const std::vector<size_t> &multiIndex) const;
protected:
	virtual std::string GetCoordinateFunction() const;
	virtual std::string GetBoundsCheckFunction() const;
	
};

}

}

#endif // CLSIM_TABULATOR_AXES_H_INCLUDED
