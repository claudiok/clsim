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
 * @file I3CLSimFunctionTester.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>

#include "test/I3CLSimFunctionTester.h"

#include <string>

#include "opencl/I3CLSimHelperLoadProgramSource.h"

I3CLSimFunctionTester::I3CLSimFunctionTester
(const I3CLSimOpenCLDevice &device,
 uint64_t workgroupSize_,
 uint64_t workItemsPerIteration_,
 I3CLSimFunctionConstPtr wlenDependentValue):
I3CLSimTesterBase(),
wlenDependentValue_(wlenDependentValue)
{
    std::vector<std::string> source;
    FillSource(source, wlenDependentValue);
    
    DoSetup(device,
            workgroupSize_,
            workItemsPerIteration_,
            source);
    
    InitBuffers();
}

void I3CLSimFunctionTester::FillSource(std::vector<std::string> &source,
                                                 I3CLSimFunctionConstPtr wlenDependentValue)
{
    source.clear();
    
    // load program source from files
    const std::string I3_BUILD(getenv("I3_BUILD"));
    const std::string kernelBaseDir = I3_BUILD+"/clsim/resources/kernels";
    
    std::string functionSource = wlenDependentValue->GetOpenCLFunction("evaluateFunction");
    
    std::string derivativeSource;
    if (wlenDependentValue->HasDerivative()) {
        derivativeSource = wlenDependentValue->GetOpenCLFunctionDerivative("evaluateDerivative");
    } else {
        derivativeSource = "\ninline float evaluateDerivative(float x) {return 888888.f;}\n";
    }

    std::string testKernelHeader = I3CLSimHelper::LoadProgramSource(kernelBaseDir+"/wlen_dep_val_test_kernel.h.cl");
    std::string testKernelSource = I3CLSimHelper::LoadProgramSource(kernelBaseDir+"/wlen_dep_val_test_kernel.c.cl");
    
    // collect the program sources
    source.push_back(testKernelHeader);
    source.push_back(functionSource);
    source.push_back(derivativeSource);
    source.push_back(testKernelSource);
}

void I3CLSimFunctionTester::InitBuffers()
{
    log_debug("Setting up device buffers.");
    // allocate empty buffers on the device
    deviceBuffer_results = boost::shared_ptr<cl::Buffer>(new cl::Buffer(*context, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, workItemsPerIteration*sizeof(float), NULL));
    deviceBuffer_inputs = boost::shared_ptr<cl::Buffer>(new cl::Buffer(*context, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, workItemsPerIteration*sizeof(float), NULL));
    log_debug("Device buffers are set up.");
    
    log_debug("Configuring kernel.");
    {
        kernel->setArg(0, *deviceBuffer_inputs);          // input data
        kernel->setArg(1, *deviceBuffer_results);          // output data
    }
    log_debug("Kernel configured.");
}


I3VectorFloatPtr I3CLSimFunctionTester::EvaluateFunction(I3VectorFloatConstPtr xValues)
{
    return EvaluateIt(xValues, false);
}

I3VectorFloatPtr I3CLSimFunctionTester::EvaluateDerivative(I3VectorFloatConstPtr xValues)
{
    if (!wlenDependentValue_) log_fatal("Internal error: wlenDependentValue_ is NULL");
    
    if (!wlenDependentValue_->HasDerivative())
        log_fatal("No function derivative available!");
    
    return EvaluateIt(xValues, true);
}


I3VectorFloatPtr I3CLSimFunctionTester::EvaluateReferenceFunction(I3VectorFloatConstPtr xValues)
{
    if (!xValues) log_fatal("NULL pointer passed to EvaluateReferenceFunction.");
    if (!wlenDependentValue_) log_fatal("Internal error: wlenDependentValue_ is NULL");
    if (!wlenDependentValue_->HasNativeImplementation()) log_fatal("No native/reference implementation available!");
    
    // allocate the output vector
    I3VectorFloatPtr results = I3VectorFloatPtr(new I3VectorFloat(xValues->size(), NAN));
    
    for (std::size_t i=0;i<xValues->size();++i)
    {
        (*results)[i] = wlenDependentValue_->GetValue((*xValues)[i]);
    }
    
    return results;
}

I3VectorFloatPtr I3CLSimFunctionTester::EvaluateReferenceDerivative(I3VectorFloatConstPtr xValues)
{
    if (!xValues) log_fatal("NULL pointer passed to EvaluateReferenceDerivative.");
    if (!wlenDependentValue_) log_fatal("Internal error: wlenDependentValue_ is NULL");
    if (!wlenDependentValue_->HasNativeImplementation()) log_fatal("No native/reference implementation available!");
    if (!wlenDependentValue_->HasDerivative())
        log_fatal("No function derivative available!");

    // allocate the output vector
    I3VectorFloatPtr results = I3VectorFloatPtr(new I3VectorFloat(xValues->size(), NAN));

    for (std::size_t i=0;i<xValues->size();++i)
    {
        (*results)[i] = wlenDependentValue_->GetDerivative((*xValues)[i]);
    }

    return results;
}


I3VectorFloatPtr I3CLSimFunctionTester::EvaluateIt(I3VectorFloatConstPtr xValues, bool derivative)
{
    if (!xValues) log_fatal("NULL pointer passed to EvaluateIt.");

    // tell the kernel whether the function or its derivative is requested
    kernel->setArg(2, derivative?1:0);

    
    // allocate the output vector
    I3VectorFloatPtr results = I3VectorFloatPtr(new I3VectorFloat());

    // if nothing to do, exit here
    if (xValues->empty()) return results;

    // reserve space for output data
    results->reserve(xValues->size());

    // determine the number of iterations
    uint64_t iterations = xValues->size()/workItemsPerIteration;
    uint64_t numEntriesInLastIteration = xValues->size()%workItemsPerIteration;
    if (numEntriesInLastIteration == 0) {
        numEntriesInLastIteration=workItemsPerIteration;
    } else {
        iterations++;
    }
    
    std::size_t inputVecPos=0;
    
    
    log_debug("Starting iterations..");
    
    for (uint64_t i=0;i<iterations;++i)
    {
        log_debug("Filling input buffer..");
        
        {
            // mapped the buffer and wait for completion
            cl::Event mappingComplete;
            float *mapped_input = (float *)queue->enqueueMapBuffer(*deviceBuffer_inputs, CL_FALSE, CL_MAP_WRITE, 0, workItemsPerIteration*sizeof(float), NULL, &mappingComplete);
            mappingComplete.wait();
            
            // add the values to the input buffer
            if (i==iterations-1)
            {
                log_debug("(in)  iteration=%" PRIu64 " (last iteration)", i);
                // last iteration
                for (uint64_t j=0;j<numEntriesInLastIteration;++j)
                {
                    mapped_input[j] = (*xValues)[inputVecPos];
                    ++inputVecPos;
                    
                }
                for (uint64_t j=numEntriesInLastIteration+1;j<workItemsPerIteration;++j)
                {
                    mapped_input[j] = 0.f;
                }
            } 
            else
            {
                log_debug("(in)  iteration=%" PRIu64, i);
                // any other iteration
                for (uint64_t j=0;j<workItemsPerIteration;++j)
                {
                    mapped_input[j] = (*xValues)[inputVecPos];
                    ++inputVecPos;
                }
            }
            
            queue->enqueueUnmapMemObject(*deviceBuffer_inputs, mapped_input);
        }
        
        log_debug("Input buffer filled.");

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
        
        log_debug("Reading output buffer..");
        {
            // mapped the buffer and wait for completion
            cl::Event mappingComplete;
            float *mapped_results = (float *)queue->enqueueMapBuffer(*deviceBuffer_results, CL_FALSE, CL_MAP_READ, 0, workItemsPerIteration*sizeof(float), NULL, &mappingComplete);
            mappingComplete.wait();
            
            // add the values to the output vector
            uint64_t entriesInThisIteration;
            if (i==iterations-1) {
                // last iteration
                entriesInThisIteration = numEntriesInLastIteration;
            } else {
                // any other iteration
                entriesInThisIteration = workItemsPerIteration;
            }

            for (uint64_t j=0;j<entriesInThisIteration;++j)
            {
                results->push_back(mapped_results[j]);
            }

            queue->enqueueUnmapMemObject(*deviceBuffer_results, mapped_results);
        }
        log_debug("Output buffer read.");

    }
    
    log_debug("iterations complete.");
    
    return results;
}
