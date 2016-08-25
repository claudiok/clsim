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
 * @file I3CLSimStepToTableConverter.h
 * @version $LastChangedRevision$
 * @date $Date$
 * @author Jakob van Santen
 */

#ifndef CLSIM_TABULATOR_STEPTOTABLECONVERTER_H_INCLUDED
#define CLSIM_TABULATOR_STEPTOTABLECONVERTER_H_INCLUDED

#include "clsim/I3CLSimStep.h"
#include "clsim/I3CLSimQueue.h"
#include "clsim/I3CLSimMediumProperties.h"
#include "clsim/I3CLSimSpectrumTable.h"
#include "clsim/random_value/I3CLSimRandomValue.h"
#include "clsim/I3CLSimOpenCLDevice.h"
#include "dataclasses/physics/I3Particle.h"

#include "clsim/tabulator/Axes.h"

#define __CL_ENABLE_EXCEPTIONS
#include "clsim/cl.hpp"

#include <boost/noncopyable.hpp>
#include <boost/thread.hpp>

class I3CLSimStepToTableConverter : boost::noncopyable {
public:
	I3CLSimStepToTableConverter(I3CLSimOpenCLDevice device,
	    clsim::tabulator::AxesConstPtr axes, size_t entriesPerStream,
	    bool storeSquaredWeights,
	    I3CLSimMediumPropertiesConstPtr medium,
	    I3CLSimSpectrumTableConstPtr spectrumTable,
	    double referenceArea,
	    I3CLSimFunctionConstPtr wavelengthAcceptance,
	    I3CLSimFunctionConstPtr angularAcceptance,
	    I3RandomServicePtr rng);
	virtual ~I3CLSimStepToTableConverter();
	void EnqueueSteps(I3CLSimStepSeriesConstPtr, I3ParticleConstPtr);
	void Finish();
	
	size_t GetBunchSize() const { return maxNumWorkitems_; }
	
	void WriteFITSFile(const std::string &fname,
	    boost::python::dict tableHeader);
private:
	
	void FetchSteps(cl::Kernel, I3RandomServicePtr);
	
	float GetBinVolume(size_t i);
	void Normalize();
	
	cl::Context context_;
	cl::CommandQueue commandQueue_;
	size_t maxWorkgroupSize_, maxNumWorkitems_, entriesPerStream_;
	
	typedef std::pair<I3CLSimStepSeriesConstPtr, I3ParticleConstPtr> bunch_t;
	I3CLSimQueue<bunch_t> stepQueue_;
	boost::thread harvesterThread_;
	bool run_;
	
	double domArea_;
	double stepLength_;
	/// group, phase
	std::pair<double, double> minimumRefractiveIndex_;
	
	clsim::tabulator::AxesConstPtr axes_;
	std::vector<float> binContent_;
	std::vector<float> squaredWeights_;
	// double rather than an integer because steps have weights
	uint64_t numPhotons_;
	double sumOfPhotonWeights_;
	/// number of Photonics photons represented by each clsim photon
	double spectralBiasFactor_;
	
	SET_LOGGER("I3CLSimStepToTableConverter");
};

I3_POINTER_TYPEDEFS(I3CLSimStepToTableConverter);

#endif // CLSIM_TABULATOR_STEPTOTABLECONVERTER_H_INCLUDED

