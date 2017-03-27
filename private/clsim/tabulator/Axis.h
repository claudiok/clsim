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
 * @file Axis.h
 * @version $LastChangedRevision$
 * @date $Date$
 * @author Jakob van Santen
 */

#ifndef CLSIM_TABULATOR_AXIS_H_INCLUDED
#define CLSIM_TABULATOR_AXIS_H_INCLUDED

#include <vector>
#include <string>

namespace clsim {

namespace tabulator {

class Axis {
public:
	Axis(double min, double max, unsigned n_bins);
	virtual ~Axis();
	
	/// lowest bin edge
	double GetMin() const { return min_; };
	/// highest bin edge
	double GetMax() const { return max_; };
	/// the number of bins *without under/overflow*
	unsigned GetNBins() const { return n_bins_; };
	
	std::string GetIndexCode(const std::string &varName) const;
	std::vector<double> GetBinEdges() const;
	double GetBinEdge(unsigned i) const;
	
	/// Transform value from linear into nonlinear space
	virtual double Transform(double value) const = 0;
	/// Transform value from nonlinear into linear space
	virtual double InverseTransform(double value) const = 0;
	
	/// Generate OpenCL code to transform variable from linear into nonlinear space
	virtual std::string GetTransformCode(const std::string &varName) const = 0;
	/// Generate OpenCL code to transform variable from nonlinear into linear space
	virtual std::string GetInverseTransformCode(const std::string &varName) const = 0;
private:
	double min_, max_;
	unsigned n_bins_;
};

class LinearAxis : public Axis {
public:
	LinearAxis(double min, double max, unsigned n_bins);
	virtual ~LinearAxis();
	
	virtual double Transform(double value) const;
	virtual double InverseTransform(double value) const;
	
	virtual std::string GetTransformCode(const std::string &varName) const;
	virtual std::string GetInverseTransformCode(const std::string &varName) const;
};

class PowerAxis : public Axis {
public:
	PowerAxis(double min, double max, unsigned n_bins, unsigned power=1);
	virtual ~PowerAxis();
	
	unsigned GetPower() const { return power_; };
	
	double Transform(double value) const;
	double InverseTransform(double value) const;
	
	virtual std::string GetTransformCode(const std::string &varName) const;
	virtual std::string GetInverseTransformCode(const std::string &varName) const;
private:
	unsigned power_;
};

}

}

#endif // CLSIM_TABULATOR_AXIS_H_INCLUDED
