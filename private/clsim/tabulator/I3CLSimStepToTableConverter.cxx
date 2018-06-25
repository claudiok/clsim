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
 * @file I3CLSimStepToTableConverter.cxx
 * @version $LastChangedRevision$
 * @date $Date$
 * @author Jakob van Santen
 */

#include "clsim/tabulator/I3CLSimStepToTableConverter.h"
#include "clsim/tabulator/Axes.h"

#include "clsim/function/I3CLSimFunctionConstant.h"

#include "opencl/I3CLSimHelperMath.h"
#include "opencl/I3CLSimHelperGenerateMediumPropertiesSource.h"
#include "opencl/I3CLSimHelperLoadProgramSource.h"
#include "opencl/mwcrng_init.h"
#include "clsim/I3CLSimModuleHelper.h"
#include "clsim/I3CLSimLightSourceToStepConverterUtils.h"
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
    const std::string I3_BUILD(getenv("I3_BUILD"));
    const std::string kernelBaseDir = I3_BUILD+"/clsim/resources/kernels/";
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
    clsim::tabulator::AxesConstPtr axes, size_t entriesPerStream, bool storeSquaredWeights,
    I3CLSimMediumPropertiesConstPtr mediumProperties, I3CLSimSpectrumTableConstPtr spectrumTable,
    double referenceArea,
    I3CLSimFunctionConstPtr wavelengthAcceptance, I3CLSimFunctionConstPtr angularAcceptance,
    I3RandomServicePtr rng) : entriesPerStream_(entriesPerStream), stepQueue_(1), run_(true),
    domArea_(referenceArea), stepLength_(1.), axes_(axes),
    numPhotons_(0), sumOfPhotonWeights_(0.)
{
	std::vector<I3CLSimRandomValueConstPtr> wavelengthGenerators;
	
	wavelengthGenerators.push_back(I3CLSimModuleHelper::makeCherenkovWavelengthGenerator
	                                (wavelengthAcceptance,
	                                 false /*generateCherenkovPhotonsWithoutDispersion_*/,
	                                 mediumProperties
	                                )
	                               );
	// Photonics tables are normalized "per photon", but those photons were
	// drawn from a Cherenkov spectrum between 300 and 600 nm and then
	// weighted by the DOM quantum efficiency (QE) at recording time. Since 
	// we draw photon wavelengths directly from a QE-weighted Cherenkov 
	// spectrum, each of our photons represents on average more than one
	// Photonics photon.
	{
		using I3CLSimLightSourceToStepConverterUtils::NumberOfPhotonsPerMeter;
		
		I3CLSimFunctionConstant uno(1.);
		I3CLSimFunctionConstPtr nPhase = mediumProperties->GetPhaseRefractiveIndex(0);
		spectralBiasFactor_ = NumberOfPhotonsPerMeter(*nPhase, uno, 300*I3Units::nanometer, 600*I3Units::nanometer)
		    / NumberOfPhotonsPerMeter(*nPhase, *wavelengthAcceptance, mediumProperties->GetMinWavelength(), mediumProperties->GetMaxWavelength());
		log_info_stream("Each photon is worth "<<spectralBiasFactor_<<" Photonics photons");
	}
	if ((spectrumTable) && (spectrumTable->size() > 1)) {
	    // a spectrum table has been configured and it contains more than the
	    // default Cherenkov spectrum at index #0.

	    for (std::size_t i=1;i<spectrumTable->size();++i)
	    {
	        wavelengthGenerators.push_back(I3CLSimModuleHelper::makeWavelengthGenerator
	                                        ((*spectrumTable)[i],
	                                         wavelengthAcceptance,
	                                         mediumProperties
	                                         )
	                                        );
	    }

	    log_debug("%zu additional (non-Cherenkov) wavelength generators (spectra) have been configured.",
	             spectrumTable->size()-1);
	}
	
	std::vector<std::string> sources;
	
	const std::string I3_BUILD(getenv("I3_BUILD"));
	const std::string kernelBaseDir = I3_BUILD+"/clsim/resources/kernels/";
	
	std::ostringstream preamble;
	preamble <<                                                               \
	    "#define SAVE_ALL_PHOTONS\n"                                          \
	    "#define SAVE_ALL_PHOTONS_PRESCALE 1\n"                               \
	    "#define PROPAGATE_FOR_FIXED_NUMBER_OF_ABSORPTION_LENGTHS 42\n"       \
	    "#define TABULATE\n"                                                  \
	    "//#define DOM_RADIUS "<<I3CLSimHelper::ToFloatString(0.16510*I3Units::m)<<"\n"\
	    "//#define PRINTF_ENABLED\n"                                            \
	;
	if (axes_->GetNDim() > 4)
		preamble << "#define TABULATE_IMPACT_ANGLE\n";
	preamble << "#define TABLE_ENTRIES_PER_STREAM " << entriesPerStream_ << "\n";
	preamble << "#define VOLUME_MODE_STEP "<<I3CLSimHelper::ToFloatString(stepLength_)<<"\n";
	minimumRefractiveIndex_ = GetMinimumRefractiveIndex(*mediumProperties);
	
	preamble << "__constant floating_t min_invGroupVel = " << I3CLSimHelper::ToFloatString(
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
	
	binContent_.resize(axes_->GetNBins(), 0);
	if (storeSquaredWeights)
		squaredWeights_.resize(axes->GetNBins(), 0);
	
#ifndef NDEBUG
	std::stringstream source;
	BOOST_FOREACH(std::string &part, sources)
		source << part;
	std::string line;
	unsigned lineno = 1;
	while (std::getline(source, line)) {
		std::cout << std::setw(4) << lineno << " " << line << std::endl;
		lineno++;
	}
#endif

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
		
		commandQueue_ = cl::CommandQueue(context_, device, CL_QUEUE_PROFILING_ENABLE);
		
		cl::Kernel kernel(program, "propKernel");
		
		maxNumWorkitems_ = std::min(kernel.getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(device), size_t(1));
		log_debug_stream("max work group size " << maxNumWorkitems_);
		log_debug_stream(device.getInfo<CL_DEVICE_NAME>() << " max memory "<<device.getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>());
		
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
	if (harvesterThread_.joinable()) {
		run_ = false;
		log_debug("Finish");
		harvesterThread_.join();
	}
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
	try {
	outputEntries = cl::Buffer(context, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR,
	    streams*entriesPerStream*sizeof(I3CLSimTableEntry));
	} catch (cl::Error &err) {
		log_error_stream(err.what() << " " << err.errstr());
		throw;
	}
	numEntries = cl::Buffer(context, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR,
	    streams*sizeof(uint32_t));

}

struct KernelStatistics {
	
	KernelStatistics() : last_timestamp_(boost::posix_time::microsec_clock::universal_time()),
	    photons_(0), steps_(0), misses_(0), hostTime_(0), deviceTime_(0)
	{}
	~KernelStatistics()
	{
		double util = 100.*double(deviceTime_)/hostTime_;
		log_info_stream("Tracked "<<photons_<<" photons in "<<(double(hostTime_)*1e-9)<<"s: "<<1e-6*double(hostTime_)/photons_ << " millisec/photon ("<<util<<"% in kernel)");
		if (misses_*100 > steps_)
			log_warn_stream(misses_<<" of "<<steps_
			    <<" photon bunches ("<<std::setprecision(2)<<((100.*misses_)/steps_)
			    <<"%) ran out of storage space. "
			    "You might want to increase the storage space per photon bunch.");
	}

	void Record(const cl::Event &kernelFinished, uint64_t n_photons, size_t steps, size_t misses)
	{
		boost::posix_time::ptime this_timestamp(boost::posix_time::microsec_clock::universal_time());
	
		uint64_t timeStart, timeEnd;
		kernelFinished.getProfilingInfo(CL_PROFILING_COMMAND_START, &timeStart);
		kernelFinished.getProfilingInfo(CL_PROFILING_COMMAND_END, &timeEnd);
	
		photons_ += n_photons;
		deviceTime_ += (timeEnd-timeStart);
		hostTime_ += (this_timestamp - last_timestamp_).total_nanoseconds();
		misses_ += misses;
		steps_ += steps;
		
		last_timestamp_ = this_timestamp;
	}
	
	void Start()
	{
		last_timestamp_ = boost::posix_time::microsec_clock::universal_time();
	}

	boost::posix_time::ptime last_timestamp_;
	size_t photons_, steps_, misses_, hostTime_, deviceTime_;
	
	SET_LOGGER("I3CLSimStepToTableConverter");
};

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
	
	KernelStatistics stats;
	
	bunch_t bunch;
	while (1) {
	
		// Work on either the remainder of the last bunch or the next bunch
		// in the queue
		if (!bunch.first && !stepQueue_.GetNonBlocking(bunch)) {
			if (run_) {
				continue;
			} else {
				return;
			}
		}
		
		size_t n_photons = 0;
		size_t real_steps = 0;
		BOOST_FOREACH(const I3CLSimStep &step, *bunch.first) {
			if (step.GetNumPhotons() > 0) {
				n_photons += step.GetNumPhotons();
				real_steps++;
			}
		}
		
		VECTOR_CLASS<cl::Event> buffersFilled(3);
		VECTOR_CLASS<cl::Event> kernelFinished(1);
		VECTOR_CLASS<cl::Event> buffersRead(3);
		
		const size_t items = bunch.first->size();
		assert(items <= maxNumWorkitems_);
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
	
		commandQueue_.enqueueReadBuffer(buffers.inputSteps, CL_FALSE, 0,
		    items*sizeof(I3CLSimStep), &osteps[0], &kernelFinished, &buffersRead[0]);
		commandQueue_.enqueueReadBuffer(buffers.numEntries, CL_FALSE, 0,
		    items*sizeof(uint32_t), &numEntries[0], &kernelFinished, &buffersRead[1]);
		commandQueue_.enqueueReadBuffer(buffers.outputEntries, CL_FALSE, 0,
		    items*entriesPerStream_*sizeof(I3CLSimTableEntry), &tableEntries[0], &kernelFinished, &buffersRead[2]);

		commandQueue_.flush();
	
		cl::Event::waitForEvents(buffersRead);
		
		I3CLSimStepSeriesPtr rsteps;
		size_t misses = 0;
		for (size_t i = 0; i < items; i++) {
			// If any steps ran out of space, return them to the work list
			if (osteps[i].GetNumPhotons() > 0) {
				log_trace_stream(osteps[i].GetNumPhotons() << " left");
				if (!rsteps)
					rsteps = I3CLSimStepSeriesPtr(new I3CLSimStepSeries);
				rsteps->push_back(osteps[i]);
				n_photons -= osteps[i].GetNumPhotons();
				misses++;
			}
		}
		
		if (rsteps && rsteps->size() > 0) {
			assert(rsteps->size() <= maxNumWorkitems_);
			// Pad out to workgroup size
			I3CLSimStep dummy = rsteps->front();
			dummy.SetNumPhotons(0);
			while (rsteps->size() < maxNumWorkitems_)
				rsteps->push_back(dummy);
			bunch.first = rsteps;
		} else {
			// No steps ran out of space; reset the work list.
			bunch.first.reset();
		}
		

		for (size_t i = 0; i < items; i++) {
			size_t size = numEntries[i];
			size_t offset = i*entriesPerStream_;
			for (size_t j = 0; j < size; j++) {
				binContent_[tableEntries[offset+j].index] += tableEntries[offset+j].weight;
			}
			if (squaredWeights_.size() > 0) {
				for (size_t j = 0; j < size; j++) {
					squaredWeights_[tableEntries[offset+j].index] += std::pow(tableEntries[offset+j].weight, 2);
				}
			}
		}
		
		stats.Record(kernelFinished[0], n_photons, real_steps, misses);
	} // while (1)
}

void
I3CLSimStepToTableConverter::Normalize()
{
	const unsigned ndim = axes_->GetNDim();
	const std::vector<size_t> shape = axes_->GetShape();
	const std::vector<size_t> strides = axes_->GetStrides();
	std::vector<size_t> idxs(ndim, 0);

	// NB: assume that the first 3 dimensions are spatial
	const size_t spatial_stride = strides[2];
	for (size_t offset = 0; offset < binContent_.size(); offset += spatial_stride) {
		// unravel index
		for (unsigned j=0; j < ndim; j++) {
			// each dimension has an under- and an overflow bin.
			idxs[j] = std::min(std::max(int(offset/strides[j] % shape[j])-1, 0), int(shape[j]-3));
			assert(idxs[j] < shape[j]-2);
		}
		assert(idxs[ndim-1] == 0);
		
		// apply volume normalization to each spatial cell
		double norm = axes_->GetBinVolume(idxs)/(stepLength_*domArea_);
		for (size_t i=0; i < spatial_stride; i++)
			binContent_[i+offset] /= norm;
		// apply to squared weights as well if needed
		if (squaredWeights_.size() > 0) {
			norm *= norm;
			for (size_t i=0; i < spatial_stride; i++)
				squaredWeights_[i+offset] /= norm;
		}
	}
}

namespace {

std::string error_text(int error)
{
	std::string text(30, '\0');
	if (error != 0) {
		fits_get_errstatus(error, &text[0]);
	}
	
	return text;
}

/*
 * Create the bin content array with transposed axis
 * counts, like PyFITS does.
 */
void create_image(fitsfile *fits, std::vector<size_t> &shape, std::string name=std::string())
{
	int error(0);
	std::vector<long> naxes(shape.size());
	std::reverse_copy(shape.begin(), shape.end(), naxes.begin());
	fits_create_img(fits, FLOAT_IMG, shape.size(), &naxes[0], &error);
	if (error != 0) {
		log_fatal_stream("Could not create image: " << error_text(error));
	}
	if (name.size() > 0) {
		fits_write_key(fits, TSTRING, "EXTNAME", (void*)(name.c_str()),
		    NULL, &error);
	}
	if (error != 0) {
		log_fatal_stream("Could name HDU "<<name<<": " << error_text(error));
	}
}

void write_pixels(fitsfile *fits, std::vector<size_t> &shape,
    std::vector<float> &pixels)
{
	int error(0);
	std::vector<long> fpixel(shape.size(), 1);
	fits_write_pix(fits, TFLOAT, &fpixel[0], pixels.size(),
	    &pixels[0], &error);
	if (error != 0) {
		log_fatal_stream("Could not fill image: " << error_text(error));
	}
}

}


void I3CLSimStepToTableConverter::WriteFITSFile(const std::string &path, boost::python::dict tableHeader)
{
	fitsfile *fits;
	int error = 0;

	fits_create_diskfile(&fits, path.c_str(), &error);
	if (error != 0) {
		log_fatal_stream("Could not create " << path << ": " << error_text(error));
	}
	
	/*
	 * Write bin content
	 */
	this->Normalize();
	{
		std::vector<size_t> shape(axes_->GetShape());
		create_image(fits, shape);
		write_pixels(fits, shape, binContent_);
	}
	
	// Fill in things that only we know
	tableHeader["n_photons"] = spectralBiasFactor_*sumOfPhotonWeights_;
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
				log_fatal_stream("Could not write header keyword "<<name.str()<<": " << error_text(error));
			}
		}
	}
	
	/*
	 * Write squared weights
	 */
	if (squaredWeights_.size() > 0) {
		std::vector<size_t> shape(axes_->GetShape());
		create_image(fits, shape, "ERRORS");
		write_pixels(fits, shape, squaredWeights_);
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
			log_fatal_stream("Could not create edge array "<<i<<": " << error_text(error));
		}
		fits_write_key(fits, TSTRING, "EXTNAME", (void*)(name.str().c_str()),
		    NULL, &error);
		if (error != 0) {
			log_fatal_stream("Could not name HDU "<<name.str()<<": " << error_text(error));
		}
		
		fits_write_pix(fits, TDOUBLE, &fpixel, size,
		    &edges[0], &error);
		if (error != 0) {
			log_fatal_stream("Could not write edge array "<<i<<": " << error_text(error));
		}
	}
	
	fits_close_file(fits, &error);
	if (error != 0) {
		log_fatal_stream("Could not close " << path << ": " << error_text(error));
	}
}
