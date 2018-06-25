/**
 * Copyright (c) 2011, 2012
 * Claudio Kopper <claudio.kopper@icecube.wisc.edu>
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
 * @file I3CLSimRandomDistributionTester.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>

#include "test/I3CLSimRandomDistributionTester.h"

#include <string>

#include "opencl/I3CLSimHelperLoadProgramSource.h"
#include "opencl/mwcrng_init.h"

#include "clsim/I3CLSimHelperToFloatString.h"
using namespace I3CLSimHelper;

I3CLSimRandomDistributionTester::I3CLSimRandomDistributionTester
(const I3CLSimOpenCLDevice &device,
 uint64_t workgroupSize_,
 uint64_t workItemsPerIteration_,
 I3RandomServicePtr randomService,
 I3CLSimRandomValueConstPtr randomDistribution,
 const std::vector<double> &runtimeParameters):
I3CLSimTesterBase()
{
    std::vector<std::string> source;
    FillSource(source, randomDistribution, runtimeParameters);
    
    DoSetup(device,
            workgroupSize_,
            workItemsPerIteration_,
            source);
    
    InitBuffers(randomService);
}

void I3CLSimRandomDistributionTester::FillSource(std::vector<std::string> &source,
                                                 I3CLSimRandomValueConstPtr randomDistribution,
                                                 const std::vector<double> &runtimeParameters)
{
    source.clear();
    
    // load program source from files
    const std::string I3_BUILD(getenv("I3_BUILD"));
    const std::string kernelBaseDir = I3_BUILD+"/clsim/resources/kernels";
    
    std::string mwcrngSource = I3CLSimHelper::LoadProgramSource(kernelBaseDir+"/mwcrng_kernel.cl");

    std::string randomDistSource;
    
    if (runtimeParameters.size() > 0) {
        randomDistSource = randomDistSource +
        "\n"
        "#define DISTRIBUTION_ARGS ";
        
        for (std::size_t i=0;i<runtimeParameters.size();++i)
        {
            if (i>0) randomDistSource = randomDistSource + ", ";
            
            randomDistSource = randomDistSource + ToFloatString(runtimeParameters[i]);
        }

        randomDistSource = randomDistSource +
        "\n";
    }
    
    randomDistSource = randomDistSource +
    "\n" +
    randomDistribution->GetOpenCLFunction
    ("generateRandomNumberAccordingToDistribution",
     // these are all defined as macros by the rng code:
     "RNG_ARGS",               // function arguments for rng
     "RNG_ARGS_TO_CALL",       // if we call anothor function, this is how we pass on the rng state
     "RNG_CALL_UNIFORM_CO",    // the call to the rng for creating a uniform number [0;1[
     "RNG_CALL_UNIFORM_OC"     // the call to the rng for creating a uniform number ]0;1]
    );
    std::string rngDistTestKernelHeader = I3CLSimHelper::LoadProgramSource(kernelBaseDir+"/rng_dist_test_kernel.h.cl");
    std::string rngDistTestKernelSource = I3CLSimHelper::LoadProgramSource(kernelBaseDir+"/rng_dist_test_kernel.c.cl");
    
    // collect the program sources
    source.push_back(rngDistTestKernelHeader);
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
    
    if (init_MWC_RNG(&(MWC_RNG_x[0]), &(MWC_RNG_a[0]), workItemsPerIteration, randomService)!=0) 
        throw std::runtime_error("I3CLSimStepToPhotonConverterOpenCL already initialized!");
    log_info("RNG is set up..");
    
    log_debug("Setting up device buffers.");
    // set up device buffers from existing host buffers
    deviceBuffer_MWC_RNG_x = boost::shared_ptr<cl::Buffer>(new cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, MWC_RNG_x.size() * sizeof(uint64_t), &(MWC_RNG_x[0])));
    deviceBuffer_MWC_RNG_a = boost::shared_ptr<cl::Buffer>(new cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, MWC_RNG_a.size() * sizeof(uint32_t), &(MWC_RNG_a[0])));

    // allocate empty buffers on the device
    deviceBuffer_results = boost::shared_ptr<cl::Buffer>(new cl::Buffer(*context, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, workItemsPerIteration*sizeof(float), NULL));
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
                                       cl::NullRange,    // current implementations force this to be NULL
                                       cl::NDRange(workItemsPerIteration),    // number of work items
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
