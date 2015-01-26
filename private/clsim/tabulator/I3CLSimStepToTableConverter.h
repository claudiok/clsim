
#ifndef CLSIM_TABULATOR_STEPTOTABLECONVERTER_H_INCLUDED
#define CLSIM_TABULATOR_STEPTOTABLECONVERTER_H_INCLUDED

#include "clsim/I3CLSimStep.h"
#include "clsim/I3CLSimQueue.h"
#include "clsim/I3CLSimMediumProperties.h"
#include "clsim/random_value/I3CLSimRandomValue.h"
#include "clsim/I3CLSimOpenCLDevice.h"
#include "clsim/cl.hpp"

#include <boost/noncopyable.hpp>

class I3CLSimStepToTableConverter : boost::noncopyable {
public:
	I3CLSimStepToTableConverter(I3CLSimOpenCLDevice device,
	    I3CLSimMediumPropertiesConstPtr medium,
	    I3CLSimFunctionConstPtr wavelengthBias,
	    I3CLSimFunctionConstPtr angularAcceptance,
	    I3RandomServicePtr rng);
	void EnqueueSteps(I3CLSimStepSeriesConstPtr);
private:
	void BuildKernel(I3CLSimOpenCLDevice device,
	    I3CLSimMediumPropertiesConstPtr mediumProperties,
	    I3CLSimFunctionConstPtr wavelengthBias);
	
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
	size_t maxWorkgroupSize_, entriesPerStream_;
	
	I3CLSimQueue<I3CLSimStepSeriesConstPtr> stepQueue_;
};

I3_POINTER_TYPEDEFS(I3CLSimStepToTableConverter);

#endif // CLSIM_TABULATOR_STEPTOTABLECONVERTER_H_INCLUDED

