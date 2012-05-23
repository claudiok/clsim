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
 * @file I3CLSimTabulator.h
 * @version $Revision$
 * @date $Date$
 * @author Jakob van Santen
 */

#ifndef CLSIM_I3CLSIMTABULATOR_H_INCLUDED
#define CLSIM_I3CLSIMTABULATOR_H_INCLUDED

#include <phys-services/I3RandomService.h>

#include <clsim/I3Photon.h>
#include <clsim/function/I3CLSimFunction.h>

#include <boost/python/object.hpp>

class I3Position;
class I3Particle;

class I3CLSimTabulator {
public:
	I3CLSimTabulator() : values_(NULL), weights_(NULL), stepLength_(0), domArea_(0)
	{
		#ifndef USE_NUMPY
		log_fatal("I can't work unless built with Numpy support!");
		#endif
	};
	
	virtual ~I3CLSimTabulator();
	
	// Each of these must be called (in no particular order) before using
	// any more of the public interface. 
	void SetBins(boost::python::object binEdges, double stepLength);
	void SetEfficiencies(I3CLSimFunctionConstPtr wavelengthAcceptance,
	    I3CLSimFunctionConstPtr angularAcceptance, double domRadius);
	void SetRandomService(I3RandomServicePtr rng);
	
	void RecordPhoton(const I3Particle &source, const I3Photon &photon);
	
	void Normalize();
	boost::python::object GetValues() const;
	
private:
	off_t GetBinIndex(const I3Particle &source, const I3Position &pos, double time) const;
	double GetBinVolume(off_t idx) const;
	
	std::vector<std::vector<double> > binEdges_;
	
	PyObject *values_;
	PyObject *weights_;
	
	double stepLength_;
	double domArea_;
	I3CLSimFunctionConstPtr angularAcceptance_;
	I3CLSimFunctionConstPtr wavelengthAcceptance_;
	I3RandomServicePtr rng_;
	
	SET_LOGGER("I3CLSimTabulator");
};

#endif /* CLSIM_I3CLSIMTABULATOR_H_INCLUDED */
