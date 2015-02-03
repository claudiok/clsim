
#include "clsim/tabulator/I3CLSimStepToTableConverter.h"
#include "clsim/tabulator/Axes.h"

#include "opencl/I3CLSimHelperMath.h"
#include "opencl/I3CLSimHelperGenerateMediumPropertiesSource.h"
#include "opencl/I3CLSimHelperLoadProgramSource.h"
#include "opencl/mwcrng_init.h"
#include "clsim/I3CLSimModuleHelper.h"
#include "clsim/I3CLSimHelperToFloatString.h"
#include "clsim/cl.hpp"

#include <fitsio.h>
#include <fitsio2.h>

#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>

namespace  {

std::string 
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
} __attribute__ ((packed));

struct I3CLSimReferenceParticle {
	I3CLSimReferenceParticle(const I3Particle &source) {
		((cl_float *)(&posAndTime))[0] = source.GetPos().GetX();
		((cl_float *)(&posAndTime))[1] = source.GetPos().GetY();
		((cl_float *)(&posAndTime))[2] = source.GetPos().GetZ();
		((cl_float *)(&posAndTime))[3] = source.GetTime();
		
		((cl_float *)(&dir))[0] = source.GetDir().GetX();
		((cl_float *)(&dir))[1] = source.GetDir().GetY();
		((cl_float *)(&dir))[2] = source.GetDir().GetZ();
		((cl_float *)(&dir))[3] = 0.;
		
		const I3Direction &dir = source.GetDir();
		double perpz = hypot(dir.GetX(), dir.GetY());
		I3Direction perpdir = (perpz > 0.) ?
		    I3Direction(-dir.GetX()*dir.GetZ()/perpz, -dir.GetY()*dir.GetZ()/perpz, perpz)
		    : I3Direction(1., 0., 0.);
		((cl_float *)(&perpDir))[0] = perpdir.GetX();
		((cl_float *)(&perpDir))[1] = perpdir.GetY();
		((cl_float *)(&perpDir))[2] = perpdir.GetZ();
		((cl_float *)(&perpDir))[3] = 0.;
		
	}
	cl_float4 posAndTime;   // x,y,z,time
	cl_float4 dir;          // dx,dy,dz,0
	cl_float4 perpDir;
} __attribute__ ((packed));

// Brute-force search for the minimum refractive index
std::pair<double, double>
GetMinimumRefractiveIndex(const I3CLSimMediumProperties &med)
{
	std::pair<double, double> n_min(
	    std::numeric_limits<double>::infinity(),
	    std::numeric_limits<double>::infinity());
	
	for (unsigned j = 0; j < med.GetLayersNum(); j++) {
		if (!med.GetGroupRefractiveIndexOverride(j))
			log_fatal("Medium properties don't know how to calculate group refractive indices");
		const I3CLSimFunction &groupIndex = *(med.GetGroupRefractiveIndexOverride(j));
		const I3CLSimFunction &phaseIndex = *(med.GetPhaseRefractiveIndex(j));
		
		double wmin = std::max(med.GetMinWavelength(), groupIndex.GetMinWlen());
		double wmax = std::min(med.GetMaxWavelength(), groupIndex.GetMaxWlen());
		for (unsigned i = 0; i < 1000; i++) {
			double n = groupIndex.GetValue(wmin + i*(wmax-wmin));
			if (n > 1 && n < n_min.first) {
				n_min = std::make_pair(n, phaseIndex.GetValue(wmin + i*(wmax-wmin)));
			}
		}
	}
	
	return n_min;
}

}

I3CLSimStepToTableConverter::I3CLSimStepToTableConverter(I3CLSimOpenCLDevice device,
    clsim::tabulator::AxesConstPtr axes,
    I3CLSimMediumPropertiesConstPtr mediumProperties,
    I3CLSimFunctionConstPtr wavelengthAcceptance, I3CLSimFunctionConstPtr angularAcceptance,
    I3RandomServicePtr rng) : entriesPerStream_(1048576), run_(true),
    domArea_(M_PI*std::pow(0.16510*I3Units::m, 2)), stepLength_(1.), axes_(axes),
    numPhotons_(0), sumOfPhotonWeights_(0.)
{
	std::vector<I3CLSimRandomValueConstPtr> wavelengthGenerators;
	wavelengthGenerators.push_back(I3CLSimModuleHelper::makeCherenkovWavelengthGenerator
	                                (wavelengthAcceptance,
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
	    "//#define PRINTF_ENABLED\n"                                            \
	;
	preamble << "#define TABLE_ENTRIES_PER_STREAM " << entriesPerStream_ << "\n";
	preamble << "#define VOLUME_MODE_STEP "<<I3CLSimHelper::ToFloatString(stepLength_)<<"\n";
	minimumRefractiveIndex_ = GetMinimumRefractiveIndex(*mediumProperties);
	
	preamble << "__constant floating_t min_invPhaseVel = " << I3CLSimHelper::ToFloatString(
	    minimumRefractiveIndex_.first/I3Constants::c) << ";\n";
	preamble << "__constant floating_t tan_thetaC = " << I3CLSimHelper::ToFloatString(
	    std::sqrt(minimumRefractiveIndex_.second*minimumRefractiveIndex_.second-1.)) << ";\n";
	
	sources.push_back(I3CLSimHelper::GetMathPreamble(device, false /* single precision for now */));
	sources.push_back(preamble.str());
	sources.push_back(I3CLSimHelper::LoadProgramSource(kernelBaseDir+"/mwcrng_kernel.cl"));
	sources.push_back(I3CLSimHelper::GenerateWavelengthGeneratorSource(wavelengthGenerators));
	sources.push_back(wavelengthAcceptance->GetOpenCLFunction("getWavelengthBias"));
	sources.push_back(I3CLSimHelper::GenerateMediumPropertiesSource(*mediumProperties));
	sources.push_back(angularAcceptance->GetOpenCLFunction("getAngularAcceptance"));
	
	sources.push_back(loadKernel("propagation_kernel", true));
	sources.push_back(axes_->GenerateBinningCode());
	sources.push_back(loadKernel("propagation_kernel", false));
	
	binContent_.resize(axes_->GetNBins());
	
	if (1) {
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
		BuildOptions += "-cl-mad-enable";

		cl::Program::Sources source;
		BOOST_FOREACH(const std::string &part, sources)
			source.push_back(std::make_pair(part.c_str(), part.size()));
		
		cl::Device device = devices[0];
		cl::Program program(context_, source);
		try {
			program.build(devices, BuildOptions.c_str());
		} catch (cl::Error &err) {
	        std::string deviceName = device.getInfo<CL_DEVICE_NAME>();
	        log_error("  * build status on %s\"", deviceName.c_str());
	        log_error("==============================");
	        log_error("Build Status: %s", boost::lexical_cast<std::string>(program.getBuildInfo<CL_PROGRAM_BUILD_STATUS>(device)).c_str());
	        log_error("Build Options: %s", boost::lexical_cast<std::string>(program.getBuildInfo<CL_PROGRAM_BUILD_OPTIONS>(device)).c_str());
	        log_error("Build Log: %s", boost::lexical_cast<std::string>(program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device)).c_str());
	        log_error("==============================");
			log_fatal_stream(err.what() << ": " << err.errstr());
		}
		
		commandQueue_ = cl::CommandQueue(context_, device, 0);
		cl::Kernel kernel(program, "propKernel");
		
		maxNumWorkitems_ = kernel.getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(device);
		log_debug_stream("max work group size " << maxNumWorkitems_);
		
		harvesterThread_ = boost::thread(boost::bind(&I3CLSimStepToTableConverter::FetchSteps, this, kernel, rng));
	}
	
}

I3CLSimStepToTableConverter::~I3CLSimStepToTableConverter()
{
	Finish();
}

void
I3CLSimStepToTableConverter::EnqueueSteps(I3CLSimStepSeriesConstPtr steps, I3ParticleConstPtr reference)
{
	if (!steps)
		return;
	BOOST_FOREACH(const I3CLSimStep &step, *steps) {
		numPhotons_ += step.GetNumPhotons();
		sumOfPhotonWeights_ += step.GetNumPhotons()*step.GetWeight();
	}
	
	stepQueue_.Put(bunch_t(steps, reference));
}

void
I3CLSimStepToTableConverter::Finish()
{
	run_ = false;
	log_debug("Finish");
	harvesterThread_.join();
	
	log_info_stream("Propagated " << numPhotons_ << " photons");
}

namespace {

struct DeviceBuffers {
	DeviceBuffers() {};
	DeviceBuffers(cl::Context, I3RandomServicePtr, size_t streams,
	    size_t entriesPerStream);
	struct {
		cl::Buffer x, a;
	} mwc; 
	cl::Buffer inputSteps;
	cl::Buffer referenceSource;
	cl::Buffer outputEntries;
	cl::Buffer numEntries;
};

DeviceBuffers::DeviceBuffers(cl::Context context,
    I3RandomServicePtr rng, size_t streams, size_t entriesPerStream)
{
	std::vector<uint64_t> xv(streams);
	std::vector<uint32_t> av(streams);
	
	init_MWC_RNG(&xv[0], &av[0], streams, rng);
	
	mwc.x = cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
	    streams*sizeof(uint64_t), &xv[0]);
	mwc.a = cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
	    streams*sizeof(uint32_t), &av[0]);
	
	inputSteps = cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR,
	    streams*sizeof(I3CLSimStep));
	referenceSource = cl::Buffer(context, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR,
	    sizeof(I3CLSimReferenceParticle));
	outputEntries = cl::Buffer(context, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR,
	    streams*entriesPerStream*sizeof(I3CLSimTableEntry));
	numEntries = cl::Buffer(context, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR,
	    streams*sizeof(uint32_t));
}

}

void
I3CLSimStepToTableConverter::FetchSteps(cl::Kernel kernel, I3RandomServicePtr rng)
{
	
	DeviceBuffers buffers(context_, rng, maxNumWorkitems_, entriesPerStream_);

	// Set kernel arguments
	uint args = 0;
	kernel.setArg(args++, buffers.inputSteps);
	kernel.setArg(args++, buffers.referenceSource);
	kernel.setArg(args++, buffers.outputEntries);
	kernel.setArg(args++, buffers.numEntries);
	kernel.setArg(args++, buffers.mwc.x);
	kernel.setArg(args++, buffers.mwc.a);
	
	I3CLSimStepSeries osteps(maxNumWorkitems_);
	std::vector<uint32_t> numEntries(maxNumWorkitems_);
	std::vector<I3CLSimTableEntry> tableEntries(maxNumWorkitems_*entriesPerStream_);
	
	while (1) {
		bunch_t bunch;
	
		if (!stepQueue_.GetNonBlocking(bunch)) {
			if (run_) {
				continue;
			} else {
				return;
			}
		}
		
		VECTOR_CLASS<cl::Event> buffersFilled(3);
		VECTOR_CLASS<cl::Event> kernelFinished(1);
		VECTOR_CLASS<cl::Event> buffersRead(3);
		
#if 0
		log_trace_stream("got " << steps->size() << " steps");
		BOOST_FOREACH(const I3CLSimStep &step, *steps)
			log_trace_stream(step.GetNumPhotons() << " photons");
#endif
		
		const size_t items = bunch.first->size();
		assert(items <= maxNumWorkItems_);
		commandQueue_.enqueueWriteBuffer(buffers.inputSteps, CL_FALSE, 0,
		    items*sizeof(I3CLSimStep), &(*bunch.first)[0], NULL, &buffersFilled[0]);
		
		I3CLSimReferenceParticle ref(*bunch.second);
		commandQueue_.enqueueWriteBuffer(buffers.referenceSource, CL_FALSE, 0,
		    sizeof(I3CLSimReferenceParticle), &ref, NULL, &buffersFilled[1]);
		
		commandQueue_.enqueueFillBuffer<uint32_t>(buffers.numEntries, 0u /*pattern*/,
		    0 /*offset*/, items*sizeof(uint32_t) /*size*/, NULL, &buffersFilled[2]);
		commandQueue_.flush();
		
		try {
		commandQueue_.enqueueNDRangeKernel(kernel, cl::NullRange,
		    cl::NDRange(items), cl::NDRange(maxNumWorkitems_),
		    &buffersFilled, &kernelFinished[0]);
		} catch (cl::Error &err) {
			log_error_stream(err.what() << " " << err.errstr());
			throw;
		}
	
		// TODO: enqueue task to consume steps when done
		commandQueue_.enqueueReadBuffer(buffers.inputSteps, CL_FALSE, 0,
		    items*sizeof(I3CLSimStep), &osteps[0], &kernelFinished, &buffersRead[0]);
		commandQueue_.enqueueReadBuffer(buffers.numEntries, CL_FALSE, 0,
		    items*sizeof(uint32_t), &numEntries[0], &kernelFinished, &buffersRead[1]);
		commandQueue_.enqueueReadBuffer(buffers.outputEntries, CL_FALSE, 0,
		    items*entriesPerStream_*sizeof(I3CLSimTableEntry), &tableEntries[0], &kernelFinished, &buffersRead[2]);

		commandQueue_.flush();
	
		cl::Event::waitForEvents(buffersRead);
		I3CLSimStepSeriesPtr rsteps;
		for (size_t i = 0; i < items; i++) {
			// If any steps ran out of space, return them to the queue to finish
			if (osteps[i].GetNumPhotons() > 0) {
				log_trace_stream(osteps[i].GetNumPhotons() << " left");
				if (!rsteps)
					rsteps = I3CLSimStepSeriesPtr(new I3CLSimStepSeries);
				rsteps->push_back(osteps[i]);
			}
		}

		if (rsteps && rsteps->size() > 0) {
			assert(rsteps->size() <= maxNumWorkItems_);
			bunch.first = rsteps;
			stepQueue_.Put(bunch);
		}

		for (size_t i = 0; i < items; i++) {
			size_t size = numEntries[i];
			size_t offset = i*entriesPerStream_;
			for (size_t j = 0; j < size; j++) {
				binContent_[tableEntries[offset+j].index] += tableEntries[offset+j].weight;
			}
		}
	} // while (1)
}

void
I3CLSimStepToTableConverter::Normalize()
{
	const unsigned ndim = axes_->GetNDim();
	const std::vector<size_t> shape = axes_->GetShape();
	const std::vector<size_t> strides = axes_->GetStrides();
	std::vector<size_t> idxs(ndim);

	// NB: assume that the first 3 dimensions are spatial
	const size_t spatial_stride = strides[2];
	for (size_t offset = 0; offset < binContent_.size(); offset += spatial_stride) {
		// unravel index
		for (unsigned j=0; j < ndim; j++)
			idxs[j] = offset/strides[j] % shape[j];
		assert(idxs[ndim-1] == 0);
		
		// apply volume normalization to each spatial cell
		double norm = axes_->GetBinVolume(idxs)/(stepLength_*domArea_);
		for (size_t i=0; i < spatial_stride; i++)
			binContent_[i+offset] /= norm;
	}
}


void I3CLSimStepToTableConverter::WriteFITSFile(const std::string &path, boost::python::dict tableHeader)
{
	fitsfile *fits;
	int error = 0;

	fits_create_diskfile(&fits, path.c_str(), &error);
	if (error != 0) {
		char err_text[30];
		fits_get_errstatus(error, err_text);
		log_fatal_stream("Could not create " << path << ": " << err_text);
	}
	
	/*
	 * Create the bin content array with transposed axis
	 * counts, like PyFITS does.
	 */
	{
		std::vector<size_t> shape(axes_->GetShape());
		std::vector<long> naxes(shape.size());
		std::reverse_copy(shape.begin(), shape.end(), naxes.begin());
		fits_create_img(fits, FLOAT_IMG, axes_->GetNDim(), &naxes[0], &error);
		if (error != 0) {
			char err_text[30];
			fits_get_errstatus(error, err_text);
			log_fatal_stream("Could not create image: " << err_text);
		}
	}
	
	/*
	 * Write bin content
	 */
	Normalize();
	{
		std::vector<long> fpixel(axes_->GetNDim(), 1);
		fits_write_pix(fits, TFLOAT, &fpixel[0], binContent_.size(),
		    &binContent_[0], &error);
		if (error != 0) {
			char err_text[30];
			fits_get_errstatus(error, err_text);
			log_fatal_stream("Could not fill image: " << err_text);
		}
	}
	
	// Fill in things that only we know
	tableHeader["n_photons"] = sumOfPhotonWeights_;
	tableHeader["n_group"] = minimumRefractiveIndex_.first;
	tableHeader["n_phase"] = minimumRefractiveIndex_.second;
	// Write header keywords
	{
		namespace bp = boost::python;

		bp::list keys = tableHeader.keys();
		for (int i = 0; i < bp::len(keys); i++) {
			bp::object key = keys[i];
			bp::object value = tableHeader[key];
			
			std::ostringstream name;
			name << "hierarch _i3_" << bp::extract<std::string>(key)();
			
			bp::extract<int> inty(value);
			bp::extract<double> doubly(value);
			if (inty.check()) {
				int v = inty();
				fits_write_key(fits, TINT, name.str().c_str(), &v,
				    NULL, &error);
			} else if (doubly.check()) {
				double v = doubly();
				fits_write_key(fits, TDOUBLE, name.str().c_str(), &v,
				    NULL, &error);
			}
			if (error != 0) {
				char err_text[30];
				fits_get_errstatus(error, err_text);
				log_fatal_stream("Could not write header keyword "<<name.str()<<": " << err_text);
			}
		}
	}
	
	/*
	 * Write each of the bin edge vectors in an extension HDU
	 */
	for (unsigned i = 0; i < axes_->GetNDim(); i++) {
		std::ostringstream name;
		name << "EDGES" << i;
		long fpixel = 1;
		std::vector<double> edges = axes_->at(i)->GetBinEdges();
		long size = edges.size();
		
		fits_create_img(fits, DOUBLE_IMG, 1, &size, &error);
		if (error != 0) {
			char err_text[30];
			fits_get_errstatus(error, err_text);
			log_fatal_stream("Could not create edge array "<<i<<": " << err_text);
		}
		fits_write_key(fits, TSTRING, "EXTNAME", (void*)(name.str().c_str()),
		    NULL, &error);
		if (error != 0) {
			char err_text[30];
			fits_get_errstatus(error, err_text);
			log_fatal_stream("Could not name HDU "<<name.str()<<": " << err_text);
		}
		
		fits_write_pix(fits, TDOUBLE, &fpixel, size,
		    &edges[0], &error);
		if (error != 0) {
			char err_text[30];
			fits_get_errstatus(error, err_text);
			log_fatal_stream("Could not write edge array "<<i<<": " << err_text);
		}
	}
	
	fits_close_file(fits, &error);
	if (error != 0) {
		char err_text[30];
		fits_get_errstatus(error, err_text);
		log_fatal_stream("Could not close " << path << ": " << err_text);
	}
}
