//
//   Copyright (c) 2011  Claudio Kopper
//   
//   $Id$
//
//   This file is part of the IceTray module g4sim-interface.
//
//   this file is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; either version 3 of the License, or
//   (at your option) any later version.
//
//   IceTray is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with this program.  If not, see <http://www.gnu.org/licenses/>.
//


#include "test/I3CLSimRandomDistributionTester.h"

#include <string>

#include "opencl/I3CLSimHelperLoadProgramSource.h"
#include "opencl/mwcrng_init.h"

I3CLSimRandomDistributionTester::I3CLSimRandomDistributionTester
(const std::pair<std::string, std::string> &platformAndDeviceName,
 uint64_t workgroupSize_,
 uint64_t workItemsPerIteration_,
 bool useNativeMath,
 I3RandomServicePtr randomService,
 I3CLSimRandomValueConstPtr randomDistribution):
I3CLSimTesterBase()
{
    std::vector<std::string> source;
    FillSource(source, randomDistribution);
    
    DoSetup(platformAndDeviceName,
            useNativeMath,
            workgroupSize_,
            workItemsPerIteration_,
            source);
    
    InitBuffers(randomService);
}

I3CLSimRandomDistributionTester::I3CLSimRandomDistributionTester
(boost::python::tuple platformAndDeviceName,
 uint64_t workgroupSize_,
 uint64_t workItemsPerIteration_,
 bool useNativeMath,
 I3RandomServicePtr randomService,
 I3CLSimRandomValueConstPtr randomDistribution):
I3CLSimTesterBase()
{
    std::string platformName = boost::python::extract<std::string>(platformAndDeviceName[0]);
    std::string deviceName = boost::python::extract<std::string>(platformAndDeviceName[1]);
    
    std::pair<std::string, std::string> argument(platformName, deviceName);
    
    std::vector<std::string> source;
    FillSource(source, randomDistribution);

    DoSetup(argument,
            useNativeMath,
            workgroupSize_,
            workItemsPerIteration_,
            source);
    
    InitBuffers(randomService);
}

void I3CLSimRandomDistributionTester::FillSource(std::vector<std::string> &source,
                                                 I3CLSimRandomValueConstPtr randomDistribution)
{
    source.clear();
    
    // load program source from files
    const std::string I3_SRC(getenv("I3_SRC"));
    const std::string kernelBaseDir = I3_SRC+"/clsim/resources/kernels";
    
    std::string mwcrngSource = I3CLSimHelper::LoadProgramSource(kernelBaseDir+"/mwcrng_kernel.cl");

    std::string randomDistSource = randomDistribution->GetOpenCLFunction("generateRandomNumberAccordingToDistribution",
                                                                         // these are all defined as macros by the rng code:
                                                                         "RNG_ARGS",               // function arguments for rng
                                                                         "RNG_ARGS_TO_CALL",       // if we call anothor function, this is how we pass on the rng state
                                                                         "RNG_CALL_UNIFORM_CO",    // the call to the rng for creating a uniform number [0;1[
                                                                         "RNG_CALL_UNIFORM_OC"     // the call to the rng for creating a uniform number ]0;1]
                                                                         );
    std::string rngDistTestKernelSource = I3CLSimHelper::LoadProgramSource(kernelBaseDir+"/rng_dist_test_kernel.cl");
    
    // collect the program sources
    source.push_back(mwcrngSource);
    source.push_back(randomDistSource);
    source.push_back(rngDistTestKernelSource);
}

void I3CLSimRandomDistributionTester::InitBuffers(I3RandomServicePtr randomService)
{
    // set up rng
    log_info("Setting up RNG for %" PRIu64 " workitems. (requested workgroupSize=%" PRIu64 ", maxWorkgroupSize=%" PRIu64 ")", workItemsPerIteration, workgroupSize, GetMaxWorkgroupSize());
    
    std::vector<uint64_t> MWC_RNG_x;
    std::vector<uint32_t> MWC_RNG_a;
    
    MWC_RNG_x.resize(workItemsPerIteration);
    MWC_RNG_a.resize(workItemsPerIteration);
    
    const std::string I3_SRC(getenv("I3_SRC"));
    if (init_MWC_RNG(&(MWC_RNG_x[0]), &(MWC_RNG_a[0]), workItemsPerIteration, (I3_SRC+"/clsim/resources/safeprimes_base32.txt").c_str(), randomService)!=0) 
        throw std::runtime_error("I3CLSimStepToPhotonConverterOpenCL already initialized!");
    log_info("RNG is set up..");
    
    log_debug("Setting up device buffers.");
    // set up device buffers from existing host buffers
    deviceBuffer_MWC_RNG_x = shared_ptr<cl::Buffer>(new cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, MWC_RNG_x.size() * sizeof(uint64_t), &(MWC_RNG_x[0])));
    deviceBuffer_MWC_RNG_a = shared_ptr<cl::Buffer>(new cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, MWC_RNG_a.size() * sizeof(uint32_t), &(MWC_RNG_a[0])));

    // allocate empty buffers on the device
    deviceBuffer_results = shared_ptr<cl::Buffer>(new cl::Buffer(*context, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, workItemsPerIteration*sizeof(float), NULL));
    log_debug("Device buffers are set up.");

    log_debug("Configuring kernel.");
    {
        kernel->setArg(0, *deviceBuffer_results);          // hit counter
        kernel->setArg(1, *deviceBuffer_MWC_RNG_x);        // rng state
        kernel->setArg(2, *deviceBuffer_MWC_RNG_a);        // rng state
    }
    log_debug("Kernel configured.");
}


I3VectorFloatPtr I3CLSimRandomDistributionTester::GenerateRandomNumbers(uint64_t iterations)
{
    // allocate the output vector
    I3VectorFloatPtr results = I3VectorFloatPtr(new I3VectorFloat());
    results->reserve(workItemsPerIteration*iterations);

    log_info("Starting iterations..");
    
    for (uint64_t i=0;i<iterations;++i)
    {
        log_debug("Starting kernel..");
        
        // run the kernel
        try {
            queue->enqueueNDRangeKernel(*kernel, 
                                       cl::NullRange,	// current implementations force this to be NULL
                                       cl::NDRange(workItemsPerIteration),	// number of work items
                                       cl::NDRange(workgroupSize),
                                       NULL,
                                       NULL);
        } catch (cl::Error &err) {
            log_fatal("OpenCL ERROR (running kernel): %s (%i)", err.what(), err.err());
        }
        
        try {
            // wait for the kernel
            queue->finish();
        } catch (cl::Error &err) {
            log_fatal("OpenCL ERROR (running kernel): %s (%i)", err.what(), err.err());
        }
        
        log_debug("kernel finished!");
        
        
        
        cl::Event mappingComplete;
        float *mapped_results = (float *)queue->enqueueMapBuffer(*deviceBuffer_results, CL_FALSE, CL_MAP_READ, 0, workItemsPerIteration*sizeof(float), NULL, &mappingComplete);
        
        // wait for the buffer to be mapped
        mappingComplete.wait();
        
        // add the values to the output vector
        for (uint64_t j=0;j<workItemsPerIteration;++j)
        {
            results->push_back(mapped_results[j]);
        }
        
        queue->enqueueUnmapMemObject(*deviceBuffer_results, mapped_results);
    }
    
    log_info("iterations complete.");
    
    return results;
}
