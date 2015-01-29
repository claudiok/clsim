
#ifndef CLSIM_TABULATOR_STEPTOTABLECONVERTER_H_INCLUDED
#define CLSIM_TABULATOR_STEPTOTABLECONVERTER_H_INCLUDED

#include "clsim/I3CLSimStep.h"
#include "clsim/I3CLSimQueue.h"
#include "clsim/I3CLSimMediumProperties.h"
#include "clsim/random_value/I3CLSimRandomValue.h"
#include "clsim/I3CLSimOpenCLDevice.h"

#define __CL_ENABLE_EXCEPTIONS
#include "clsim/cl.hpp"

#include <boost/noncopyable.hpp>
#include <boost/thread.hpp>

class I3CLSimStepToTableConverter : boost::noncopyable {
public:
	I3CLSimStepToTableConverter(I3CLSimOpenCLDevice device,
	    I3CLSimMediumPropertiesConstPtr medium,
	    I3CLSimFunctionConstPtr wavelengthBias,
	    I3CLSimFunctionConstPtr angularAcceptance,
	    I3RandomServicePtr rng);
	virtual ~I3CLSimStepToTableConverter();
	void EnqueueSteps(I3CLSimStepSeriesConstPtr);
	void Finish();
private:
	
	void FetchSteps();
	void FetchEntries(size_t nsteps);
	
	struct DeviceBuffers {
		DeviceBuffers() {};
		DeviceBuffers(cl::Context, I3RandomServicePtr, size_t streams,
		    size_t entriesPerStream);
		struct {
			cl::Buffer x, a;
		} mwc; 
		cl::Buffer inputSteps;
		cl::Buffer outputEntries;
		cl::Buffer numEntries;
	};
	DeviceBuffers buffers_;
	
	cl::Context context_;
	cl::CommandQueue commandQueue_;
	cl::Kernel propagationKernel_;
	size_t maxWorkgroupSize_, maxNumWorkitems_, entriesPerStream_;
	
	I3CLSimQueue<I3CLSimStepSeriesConstPtr> stepQueue_;
	boost::thread harvesterThread_;
	bool run_;
	
	std::vector<float> binContent_;
	// double rather than an integer because steps have weights
	uint64_t numPhotons_;
	double sumOfPhotonWeights_;
	
	SET_LOGGER("I3CLSimStepToTableConverter");
};

I3_POINTER_TYPEDEFS(I3CLSimStepToTableConverter);

#endif // CLSIM_TABULATOR_STEPTOTABLECONVERTER_H_INCLUDED

