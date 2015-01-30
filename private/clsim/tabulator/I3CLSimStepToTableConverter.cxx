
#include "clsim/tabulator/I3CLSimStepToTableConverter.h"
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
	
	inputSteps = cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR,
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
	std::string GetIndexFunction(std::string var) const {
		std::ostringstream ss;
		std::string scale = GetIndex(var);
		ss << "("<<var<<" >= "<<I3CLSimHelper::ToFloatString(max_)<<" ? "<<(n_bins_-1)<<"u"
		    <<" : "<<"("<<var<<" < "<<I3CLSimHelper::ToFloatString(min_)<<" ? 0u : (uint)("<<scale<<"))" <<")"
		;
		return ss.str();
	}
	std::string GetIndex(std::string var) const
	{
		std::ostringstream ss;
		ss << (n_bins_) 
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
GenerateIndex(const std::string &functionName, const std::vector<Axis> &axes)
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

std::string
PrintMultiIndex(const std::string &functionName, const std::vector<Axis> &axes)
{
	std::ostringstream ss;
	ss << "inline uint " << functionName << "(";
	for (uint i = 0; i < axes.size(); i++) {
		ss << "floating_t c" << i;
		if (i+1 < axes.size())
			ss << ", ";
	}
	ss << ")\n{\n";
	ss << "dbg_printf(\"[ ";
	for (uint i = 0; i < axes.size(); i++)
		ss << "%d ";
	ss << "]\\n\", ";
	for (uint i = 0; i < axes.size(); i++) {
		std::ostringstream var;
		var << "c" << i;
		ss << axes[i].GetIndexFunction(var.str());
		if (i+1 < axes.size())
			ss << ", ";
	}
	ss << ");\n}\n";
	
	return ss.str();
}

std::string
GenerateBoundsCheck(const std::string &functionName, const std::vector<Axis> &axes)
{
	std::ostringstream ss;
	ss << "inline bool " << functionName << "(";
	for (uint i = 0; i < axes.size(); i++) {
		ss << "floating_t c" << i;
		if (i+1 < axes.size())
			ss << ", ";
	}
	ss << ")\n{\n";
	ss << "return (c3 > "<<I3CLSimHelper::ToFloatString(axes.at(3).max_)<<");";
	ss << "\n}\n";
	
	return ss.str();
}

std::string
GenerateSourcePosition(const I3Particle &source)
{
	std::ostringstream ss;
	ss << "static __constant floating4_t sourcePos = {"
	    <<I3CLSimHelper::ToFloatString(source.GetPos().GetX())<<","
	    <<I3CLSimHelper::ToFloatString(source.GetPos().GetY())<<","
	    <<I3CLSimHelper::ToFloatString(source.GetPos().GetZ())<<","
	    <<I3CLSimHelper::ToFloatString(source.GetTime())<<"};\n";
	const I3Direction &dir = source.GetDir();
	ss << "static __constant floating4_t sourceDir = {"
	    <<I3CLSimHelper::ToFloatString(dir.GetX())<<","
	    <<I3CLSimHelper::ToFloatString(dir.GetY())<<","
	    <<I3CLSimHelper::ToFloatString(dir.GetZ())<<","
	    <<I3CLSimHelper::ToFloatString(0.)<<"};\n";
	double perpz = hypot(dir.GetX(), dir.GetY());
	I3Direction perpdir = (perpz > 0.) ?
	    I3Direction(-dir.GetX()*dir.GetZ()/perpz, -dir.GetY()*dir.GetZ()/perpz, perpz)
	    : I3Direction(1., 0., 0.);
	ss << "static __constant floating4_t perpDir = {"
	    <<I3CLSimHelper::ToFloatString(perpdir.GetX())<<","
	    <<I3CLSimHelper::ToFloatString(perpdir.GetY())<<","
	    <<I3CLSimHelper::ToFloatString(perpdir.GetZ())<<","
	    <<I3CLSimHelper::ToFloatString(0.)<<"};\n";
	return ss.str();
}

// Brute-force search for the minimum refractive index
double
GetMinimumRefractiveIndex(const I3CLSimMediumProperties &med)
{
	double n_min = std::numeric_limits<double>::infinity();
	

	
	for (unsigned j = 0; j < med.GetLayersNum(); j++) {
		if (!med.GetGroupRefractiveIndexOverride(j))
			log_fatal("Medium properties don't know how to calculate group refractive indices");
		const I3CLSimFunction &groupIndex = *(med.GetGroupRefractiveIndexOverride(j));
		double wmin = std::max(med.GetMinWavelength(), groupIndex.GetMinWlen());
		double wmax = std::min(med.GetMaxWavelength(), groupIndex.GetMaxWlen());
		for (unsigned i = 0; i < 1000; i++) {
			double n = groupIndex.GetValue(wmin + i*(wmax-wmin));
			if (n > 1 && n < n_min) {
				n_min = n;
			}
		}
	}
	
	return n_min;
}

I3CLSimStepToTableConverter::I3CLSimStepToTableConverter(I3CLSimOpenCLDevice device,
    const I3Particle &referenceSource,
    I3CLSimMediumPropertiesConstPtr mediumProperties,
    I3CLSimFunctionConstPtr wavelengthAcceptance, I3CLSimFunctionConstPtr angularAcceptance,
    I3RandomServicePtr rng) : entriesPerStream_(1048576), run_(true),
    domArea_(M_PI*std::pow(0.16510*I3Units::m, 2)), stepLength_(1.),
    referenceSource_(referenceSource),
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
	
	preamble << "#define MIN_INV_GROUPVEL " << I3CLSimHelper::ToFloatString(
	    minimumRefractiveIndex_/I3Constants::c) << "\n";
	
	sources.push_back(I3CLSimHelper::GetMathPreamble(device, false /* single precision for now */));
	sources.push_back(preamble.str());
	sources.push_back(I3CLSimHelper::LoadProgramSource(kernelBaseDir+"/mwcrng_kernel.cl"));
	sources.push_back(I3CLSimHelper::GenerateWavelengthGeneratorSource(wavelengthGenerators));
	sources.push_back(wavelengthAcceptance->GetOpenCLFunction("getWavelengthBias"));
	sources.push_back(I3CLSimHelper::GenerateMediumPropertiesSource(*mediumProperties));
	sources.push_back(angularAcceptance->GetOpenCLFunction("getAngularAcceptance"));
	sources.push_back(GenerateSourcePosition(referenceSource_));
	std::cout << GenerateSourcePosition(referenceSource_) << std::endl;;
	
	std::vector<Axis> axes;
	axes.push_back(Axis(0, 580, 200, 2));
	axes.push_back(Axis(0, 180, 36, 1));
	axes.push_back(Axis(-1, 1, 100, 1));
	axes.push_back(Axis(0, 7e3, 105, 2));
	
	size_t total_size = 1;
	BOOST_FOREACH(const Axis &axis, axes) {
		binEdges_.push_back(std::vector<double>());
		double imin = std::pow(axis.min_, 1./axis.power_);
		double imax = std::pow(axis.max_, 1./axis.power_);
		double istep = (imax-imin)/(axis.n_bins_);
		for (unsigned i = 0; i <= axis.n_bins_; i++) {
			binEdges_.back().push_back(std::pow(imin + i*istep, axis.power_));
		}
		
		total_size *= axis.n_bins_;
	}
	binContent_.resize(total_size);
	
	sources.push_back(GenerateBoundsCheck("isOutOfBounds", axes));
	sources.push_back(GenerateIndex("getBinIndex", axes));
	
	sources.push_back(loadKernel("propagation_kernel", true));
	sources.push_back(loadKernel("propagation_kernel", false));
	
	if (0) {
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
		propagationKernel_ = cl::Kernel(program, "propKernel");
		
		maxWorkgroupSize_ = propagationKernel_.getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(device);
		log_debug_stream("max work group size " << maxWorkgroupSize_);
		maxWorkgroupSize_ = 1;
		
		
		buffers_ = DeviceBuffers(context_, rng, maxWorkgroupSize_, entriesPerStream_);
		
		// Set arguments
		uint args = 0;
		propagationKernel_.setArg(args++, buffers_.inputSteps);
		propagationKernel_.setArg(args++, buffers_.outputEntries);
		propagationKernel_.setArg(args++, buffers_.numEntries);
		propagationKernel_.setArg(args++, buffers_.mwc.x);
		propagationKernel_.setArg(args++, buffers_.mwc.a);
	}
	
	harvesterThread_ = boost::thread(boost::bind(&I3CLSimStepToTableConverter::FetchSteps, this));
	// exit(1);
}

I3CLSimStepToTableConverter::~I3CLSimStepToTableConverter()
{
	Finish();
}

void
I3CLSimStepToTableConverter::EnqueueSteps(I3CLSimStepSeriesConstPtr steps)
{
	if (!steps)
		return;
	BOOST_FOREACH(const I3CLSimStep &step, *steps) {
		numPhotons_ += step.GetNumPhotons();
		sumOfPhotonWeights_ += step.GetNumPhotons()*step.GetWeight();
	}
	
	stepQueue_.Put(steps);
}

void
I3CLSimStepToTableConverter::Finish()
{
	run_ = false;
	log_debug("Finish");
	harvesterThread_.join();
	
	log_info_stream("Propagated " << numPhotons_ << " photons");
}

void
I3CLSimStepToTableConverter::FetchSteps()
{
	I3CLSimStepSeries osteps(maxWorkgroupSize_);
	std::vector<uint32_t> numEntries(maxWorkgroupSize_);
	std::vector<I3CLSimTableEntry> tableEntries(maxWorkgroupSize_*entriesPerStream_);
	
	while (1) {
		I3CLSimStepSeriesConstPtr steps;
	
		if (!stepQueue_.GetNonBlocking(steps)) {
			if (run_) {
				continue;
			} else {
				return;
			}
		}
		
		VECTOR_CLASS<cl::Event> buffersFilled(2);
		VECTOR_CLASS<cl::Event> kernelFinished(1);
		VECTOR_CLASS<cl::Event> buffersRead(3);
		
#if 0
		log_trace_stream("got " << steps->size() << " steps");
		BOOST_FOREACH(const I3CLSimStep &step, *steps)
			log_trace_stream(step.GetNumPhotons() << " photons");
#endif
		
		const size_t items = steps->size();
		assert(items <= maxWorkgroupSize_);
		commandQueue_.enqueueWriteBuffer(buffers_.inputSteps, CL_FALSE, 0,
		    items*sizeof(I3CLSimStep), &(*steps)[0], NULL, &buffersFilled[0]);
		commandQueue_.enqueueFillBuffer<uint32_t>(buffers_.numEntries, 0u /*pattern*/,
		    0 /*offset*/, items*sizeof(uint32_t) /*size*/, NULL, &buffersFilled[1]);
		commandQueue_.flush();

		commandQueue_.enqueueNDRangeKernel(propagationKernel_, cl::NullRange,
		    cl::NDRange(items), cl::NDRange(maxWorkgroupSize_),
		    &buffersFilled, &kernelFinished[0]);
	
		// TODO: enqueue task to consume steps when done
		commandQueue_.enqueueReadBuffer(buffers_.inputSteps, CL_FALSE, 0,
		    items*sizeof(I3CLSimStep), &osteps[0], &kernelFinished, &buffersRead[0]);
		commandQueue_.enqueueReadBuffer(buffers_.numEntries, CL_FALSE, 0,
		    items*sizeof(uint32_t), &numEntries[0], &kernelFinished, &buffersRead[1]);
		commandQueue_.enqueueReadBuffer(buffers_.outputEntries, CL_FALSE, 0,
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
#if 1
		if (rsteps && rsteps->size() > 0) {
			assert(rsteps->size() <= maxWorkgroupSize_);
			stepQueue_.Put(rsteps);
		}
		
		for (size_t i = 0; i < items; i++) {
			size_t size = numEntries[i];
			size_t offset = i*entriesPerStream_;
			for (size_t j = 0; j < size; j++) {
				binContent_[tableEntries[offset+j].index] += tableEntries[offset+j].weight;
			}
		}
		
		
#endif
	}
}

// FIXME: factor binning-related calculations out into a worker class
void
I3CLSimStepToTableConverter::Normalize()
{
	std::vector<size_t> strides(binEdges_.size(), 1);
	std::vector<size_t> dims(binEdges_.size(), 1);
	for (unsigned i = 0; i < binEdges_.size(); i++)
		dims[i] = binEdges_[i].size()-1;
	for (int i = binEdges_.size()-2; i >= 0; i--)
		strides[i] *= strides[i+1]*dims[i+1];
	
	for (size_t idx = 0; idx < binContent_.size(); idx++) {
		// unravel index
		size_t idxs[4];
		for (int j = 0; j < 4; j++) {
			idxs[j] = idx/strides[j] % dims[j];
		}
		
		// NB: since we combine the bins at azimuth > 180 degrees with the other half of
		// the sphere, the true volume of an azimuthal bin is twice its nominal value.
		double volume = ((std::pow(binEdges_[0][idxs[0]+1], 3) - std::pow(binEdges_[0][idxs[0]], 3))/3.)
		    * 2*I3Units::degree*(binEdges_[1][idxs[1]+1] - binEdges_[1][idxs[1]])
		    * (binEdges_[2][idxs[2]+1] - binEdges_[2][idxs[2]]);
		double norm = volume/(stepLength_*domArea_);		
		binContent_[idx] /= norm;
	}
}


void I3CLSimStepToTableConverter::WriteFITSFile(const std::string &path, boost::python::dict tableHeader)
{
	fitsfile *fits;
	int error = 0;
	char *err_text;
	
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
		std::vector<long> naxes(binEdges_.size());
		for (unsigned i = 0; i < binEdges_.size(); i++)
			naxes[i] = binEdges_[binEdges_.size() - i - 1].size()-1;
	
		fits_create_img(fits, FLOAT_IMG, binEdges_.size(), &naxes[0], &error);
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
		std::vector<long> fpixel(binEdges_.size(), 1);
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
	tableHeader["n_group"] = minimumRefractiveIndex_;
	tableHeader["z"] = referenceSource_.GetPos().GetZ()/I3Units::m;
	tableHeader["zenith"] = referenceSource_.GetDir().GetZenith()/I3Units::degree;
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
	for (unsigned i = 0; i < binEdges_.size(); i++) {
		std::ostringstream name;
		name << "EDGES" << i;
		long fpixel = 1;
		long size = binEdges_[i].size();
		
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
		    &binEdges_[i][0], &error);
		if (error != 0) {
			char err_text[30];
			fits_get_errstatus(error, err_text);
			log_fatal_stream("Could not write edge array "<<i<<": " << err_text);
		}
	}
	
	// TODO write header keywords
	
	fits_close_file(fits, &error);
	if (error != 0) {
		char err_text[30];
		fits_get_errstatus(error, err_text);
		log_fatal_stream("Could not close " << path << ": " << err_text);
	}
}


// void
// I3CLSimStepToTableConverter::FetchEntries(size_t nsteps)
// {
// 	std::vector<I3CLSimStep> steps(nsteps);
// 	commandQueue_.enqueueReadBuffer(buffers_.inputSteps, CL_TRUE, 0, nsteps*sizeof(I3CLSimStep),
// 	    )
// }