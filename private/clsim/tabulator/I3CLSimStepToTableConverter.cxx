
#include "clsim/tabulator/I3CLSimStepToTableConverter.h"
#include "opencl/I3CLSimHelperMath.h"
#include "opencl/I3CLSimHelperGenerateMediumPropertiesSource.h"
#include "opencl/I3CLSimHelperLoadProgramSource.h"
#include "opencl/mwcrng_init.h"
#include "clsim/I3CLSimModuleHelper.h"
#include "clsim/I3CLSimHelperToFloatString.h"
#include "clsim/cl.hpp"

#include <boost/foreach.hpp>

static std::string 
loadKernel(const std::string& name, bool header)
{
    const std::string I3_SRC(getenv("I3_SRC"));
    const std::string kernelBaseDir = I3_SRC+"/clsim/resources/kernels/";
    const std::string ext = header ? ".h.cl" : ".c.cl";
    return I3CLSimHelper::LoadProgramSource(kernelBaseDir+name+ext);
}

struct I3CLSimTableEntry {
	uint32_t index;
	float weight;
};

I3CLSimStepToTableConverter::DeviceBuffers::DeviceBuffers(cl::Context context,
    I3RandomServicePtr rng, size_t streams, size_t entriesPerStream)
{
	std::vector<uint64_t> xv(streams);
	std::vector<uint32_t> av(streams);
	
	init_MWC_RNG(&xv[0], &av[0], streams, rng);
	
	mwc.x = cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
	    streams*sizeof(uint64_t), &xv[0]);
	mwc.a = cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
	    streams*sizeof(uint32_t), &av[0]);
	
	inputSteps = cl::Buffer(context, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR,
	    streams*sizeof(I3CLSimStep));
	outputEntries = cl::Buffer(context, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR,
	    streams*entriesPerStream*sizeof(I3CLSimTableEntry));
	numEntries = cl::Buffer(context, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR,
	    streams*sizeof(uint32_t));
}

struct Axis {
	Axis(double min, double max, size_t n_bins, uint power = 1)
		: min_(min), max_(max), n_bins_(n_bins), power_(power)
	{};
	std::string transform(std::string var) const
	{
		std::ostringstream ss;
		if (power_ < 4) {
			ss << var;
			for (uint i = 1; i < power_; i++)
				ss << "*" << var;
		} else {
			ss << "pow(" << var << "," << I3CLSimHelper::ToFloatString(power_) << ")";
		}
		return ss.str();
	}
	std::string itransform(std::string var) const
	{
		std::ostringstream ss;
		if (power_ == 1) {
			ss << var;
		} else if (power_ == 2) {
			ss << "sqrt(" << var << ")";
		} else if (power_ == 3) {
			ss << "cbrt(" << var << ")";
		} else {
			ss << "pow(" << var << "," << I3CLSimHelper::ToFloatString(1./power_) << ")";
		}
		return ss.str();
	}
	std::string GetIndexFunction(std::string var) const
	{
		std::ostringstream ss;
		ss << (n_bins_-1)
		    << "*(" << itransform(var) << " - " << itransform(I3CLSimHelper::ToFloatString(min_)) << ")"
		    << "/(" << itransform(I3CLSimHelper::ToFloatString(max_)) << " - " << itransform(I3CLSimHelper::ToFloatString(min_)) << ")"
		;
		return ss.str();
	}
	double min_, max_;
	size_t n_bins_;
	uint power_;
};


std::string
GenerateBinningCode(const std::string &functionName, const std::vector<Axis> &axes)
{
	std::ostringstream ss;
	ss << "inline uint " << functionName << "(";
	for (uint i = 0; i < axes.size(); i++) {
		ss << "floating_t c" << i;
		if (i+1 < axes.size())
			ss << ", ";
	}
	ss << ")\n{\n";
	ss << "return ";
	
	std::vector<size_t> strides(axes.size(), 1);
	for (int i = axes.size()-2; i >= 0; i--)
		strides[i] *= strides[i+1]*axes[i+1].n_bins_;
	for (uint i = 0; i < axes.size(); i++) {
		std::ostringstream var;
		var << "c" << i;
		ss << strides[i] << "*" << axes[i].GetIndexFunction(var.str());
		if (i+1 < axes.size())
			ss << "\n + ";
	}
	
	ss << ";\n}\n";
	
	return ss.str();
}

I3CLSimStepToTableConverter::I3CLSimStepToTableConverter(I3CLSimOpenCLDevice device,
    I3CLSimMediumPropertiesConstPtr mediumProperties,
    I3CLSimFunctionConstPtr wavelengthBias, I3CLSimFunctionConstPtr angularAcceptance,
    I3RandomServicePtr rng) : entriesPerStream_(4096)
{
	std::vector<I3CLSimRandomValueConstPtr> wavelengthGenerators;
	wavelengthGenerators.push_back(I3CLSimModuleHelper::makeCherenkovWavelengthGenerator
	                                (wavelengthBias,
	                                 false /*generateCherenkovPhotonsWithoutDispersion_*/,
	                                 mediumProperties
	                                )
	                               );
	
	std::vector<std::string> sources;
	
	const std::string I3_SRC(getenv("I3_SRC"));
	const std::string kernelBaseDir = I3_SRC+"/clsim/resources/kernels/";
	
	std::ostringstream preamble;
	preamble <<                                                               \
	    "#define SAVE_ALL_PHOTONS\n"                                          \
	    "#define SAVE_ALL_PHOTONS_PRESCALE 1\n"                               \
	    "#define PROPAGATE_FOR_FIXED_NUMBER_OF_ABSORPTION_LENGTHS 42\n"       \
	    "#define TABULATE\n"                                                  \
	;
	preamble << "#define TABLE_ENTRIES_PER_STREAM " << entriesPerStream_ << "\n";
	preamble << "#define VOLUME_MODE_STEP 1.0f\n";
	
	sources.push_back(I3CLSimHelper::GetMathPreamble(device, false /* single precision for now */));
	sources.push_back(preamble.str());
	sources.push_back(I3CLSimHelper::LoadProgramSource(kernelBaseDir+"/mwcrng_kernel.cl"));
	sources.push_back(I3CLSimHelper::GenerateWavelengthGeneratorSource(wavelengthGenerators));
	sources.push_back(wavelengthBias->GetOpenCLFunction("getWavelengthBias"));
	sources.push_back(I3CLSimHelper::GenerateMediumPropertiesSource(*mediumProperties));
	sources.push_back(angularAcceptance->GetOpenCLFunction("getAngularAcceptance"));
	
	std::vector<Axis> axes;
	axes.push_back(Axis(0, 580, 200, 2));
	axes.push_back(Axis(0, 180, 36, 1));
	axes.push_back(Axis(-1, 1, 100, 1));
	axes.push_back(Axis(0, 7e3, 105, 2));
	sources.push_back(GenerateBinningCode("getBinIndex", axes));
	
	sources.push_back(loadKernel("propagation_kernel", true));
	// if (!saveAllPhotons_) {
	// 	propagationKernelSource += this->GetCollisionDetectionSource(true);
	// 	propagationKernelSource += this->GetCollisionDetectionSource(false);
	// }
	sources.push_back(loadKernel("propagation_kernel", false));
	
	{
		std::stringstream source;
		BOOST_FOREACH(std::string &part, sources)
			source << part;
		std::string line;
		unsigned lineno = 1;
		while (std::getline(source, line)) {
			std::cout << std::setw(4) << lineno << " " << line << std::endl;
			lineno++;
		}
	}

	{
		VECTOR_CLASS<cl::Device> devices(1, *(device.GetDeviceHandle()));
		cl_context_properties properties[] = 
		{ CL_CONTEXT_PLATFORM, (cl_context_properties)(*(device.GetPlatformHandle()))(), 0};
		
		context_ = cl::Context(devices, properties);
		
		std::string BuildOptions;
		BuildOptions += "-cl-mad-enable ";

		cl::Program::Sources source;
		BOOST_FOREACH(const std::string &part, sources)
			source.push_back(std::make_pair(part.c_str(), part.size()));
		
		cl::Program program(context_, source);
		program.build(devices, BuildOptions.c_str());
		
		cl::Device device = devices[0];
        std::string deviceName = device.getInfo<CL_DEVICE_NAME>();
        log_error("  * build status on %s\"", deviceName.c_str());
        log_error("==============================");
        log_error("Build Status: %s", boost::lexical_cast<std::string>(program.getBuildInfo<CL_PROGRAM_BUILD_STATUS>(device)).c_str());
        log_error("Build Options: %s", boost::lexical_cast<std::string>(program.getBuildInfo<CL_PROGRAM_BUILD_OPTIONS>(device)).c_str());
        log_error("Build Log: %s", boost::lexical_cast<std::string>(program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device)).c_str());
        log_error("==============================");
		
		commandQueue_ = cl::CommandQueue(context_, device, 0);
		propagationKernel_ = cl::Kernel(program, "propKernel");
		
		maxWorkgroupSize_ = propagationKernel_.getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(device);
		log_debug_stream("max work group size " << maxWorkgroupSize_);
		
		
		buffers_ = DeviceBuffers(context_, rng, maxWorkgroupSize_, entriesPerStream_);
	}
	
	
	exit(1);
}

void
I3CLSimStepToTableConverter::EnqueueSteps(I3CLSimStepSeriesConstPtr steps)
{
	stepQueue_.Put(steps);
}