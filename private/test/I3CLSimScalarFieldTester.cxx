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
 * @file I3CLSimScalarFieldTester.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>

#include "test/I3CLSimScalarFieldTester.h"

#include <string>

#include "opencl/I3CLSimHelperLoadProgramSource.h"

I3CLSimScalarFieldTester::I3CLSimScalarFieldTester
(const I3CLSimOpenCLDevice &device,
 uint64_t workgroupSize_,
 uint64_t workItemsPerIteration_,
 I3CLSimScalarFieldConstPtr theField):
I3CLSimTesterBase(),
theField_(theField)
{
    std::vector<std::string> source;
    FillSource(source, theField);
    
    DoSetup(device,
            workgroupSize_,
            workItemsPerIteration_,
            source);
    
    InitBuffers();
}

void I3CLSimScalarFieldTester::FillSource(
    std::vector<std::string> &source,
    I3CLSimScalarFieldConstPtr theField)
{
    source.clear();
    
    // load program source from files
    const std::string I3_BUILD(getenv("I3_BUILD"));
    const std::string kernelBaseDir = I3_BUILD+"/clsim/resources/kernels";
    
    std::string functionSource = theField->GetOpenCLFunction("evaluateScalarField");
    
    std::string testKernelHeader = I3CLSimHelper::LoadProgramSource(kernelBaseDir+"/scalar_field_test_kernel.h.cl");
    std::string testKernelSource = I3CLSimHelper::LoadProgramSource(kernelBaseDir+"/scalar_field_test_kernel.c.cl");
    
    // collect the program sources
    source.push_back(testKernelHeader);
    source.push_back(functionSource);
    source.push_back(testKernelSource);
}

void I3CLSimScalarFieldTester::InitBuffers()
{
    log_debug("Setting up device buffers.");
    // allocate empty buffers on the device
    deviceBuffer_results  = boost::shared_ptr<cl::Buffer>(new cl::Buffer(*context, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, workItemsPerIteration*sizeof(float), NULL));
    deviceBuffer_inputs_x = boost::shared_ptr<cl::Buffer>(new cl::Buffer(*context, CL_MEM_READ_ONLY  | CL_MEM_ALLOC_HOST_PTR, workItemsPerIteration*sizeof(float), NULL));
    deviceBuffer_inputs_y = boost::shared_ptr<cl::Buffer>(new cl::Buffer(*context, CL_MEM_READ_ONLY  | CL_MEM_ALLOC_HOST_PTR, workItemsPerIteration*sizeof(float), NULL));
    deviceBuffer_inputs_z = boost::shared_ptr<cl::Buffer>(new cl::Buffer(*context, CL_MEM_READ_ONLY  | CL_MEM_ALLOC_HOST_PTR, workItemsPerIteration*sizeof(float), NULL));
    log_debug("Device buffers are set up.");
    
    log_debug("Configuring kernel.");
    {
        kernel->setArg(0, *deviceBuffer_inputs_x);        // input data
        kernel->setArg(1, *deviceBuffer_inputs_y);        // input data
        kernel->setArg(2, *deviceBuffer_inputs_z);        // input data
        kernel->setArg(3, *deviceBuffer_results);         // output data
    }
    log_debug("Kernel configured.");
}


I3VectorFloatPtr I3CLSimScalarFieldTester::EvaluateReferenceFunction(
    I3VectorFloatConstPtr xValues,
    I3VectorFloatConstPtr yValues,
    I3VectorFloatConstPtr zValues)
{
    if (!xValues) log_fatal("NULL pointer passed to EvaluateReferenceFunction. (xValues)");
    if (!yValues) log_fatal("NULL pointer passed to EvaluateReferenceFunction. (yValues)");
    if (!zValues) log_fatal("NULL pointer passed to EvaluateReferenceFunction. (zValues)");
    if (xValues->size() != yValues->size()) log_fatal("the input vectors have to have the same size!");
    if (xValues->size() != zValues->size()) log_fatal("the input vectors have to have the same size!");

    if (!theField_) log_fatal("Internal error: theField_ is NULL");
    if (!theField_->HasNativeImplementation()) log_fatal("No native/reference implementation available!");
    
    // allocate the output vector
    I3VectorFloatPtr results = I3VectorFloatPtr(new I3VectorFloat(xValues->size(), NAN));
    
    for (std::size_t i=0;i<xValues->size();++i)
    {
        (*results)[i] = theField_->GetValue((*xValues)[i], (*yValues)[i], (*zValues)[i]);
    }
    
    return results;
}

I3VectorFloatPtr I3CLSimScalarFieldTester::EvaluateFunction(
    I3VectorFloatConstPtr xValues,
    I3VectorFloatConstPtr yValues,
    I3VectorFloatConstPtr zValues)
{
    if (!xValues) log_fatal("NULL pointer passed to EvaluateFunction. (xValues)");
    if (!yValues) log_fatal("NULL pointer passed to EvaluateFunction. (yValues)");
    if (!zValues) log_fatal("NULL pointer passed to EvaluateFunction. (zValues)");
    if (xValues->size() != yValues->size()) log_fatal("the input vectors have to have the same size!");
    if (xValues->size() != zValues->size()) log_fatal("the input vectors have to have the same size!");

    // allocate the output vector
    I3VectorFloatPtr results = I3VectorFloatPtr(new I3VectorFloat());

    // if nothing to do, exit here
    if (xValues->empty()) return results;

    // reserve space for output data
    results->resize(xValues->size());

    // determine the number of iterations
    uint64_t iterations = xValues->size()/workItemsPerIteration;
    uint64_t numEntriesInLastIteration = xValues->size()%workItemsPerIteration;
    if (numEntriesInLastIteration == 0) {
        numEntriesInLastIteration=workItemsPerIteration;
    } else {
        iterations++;
    }
    
    std::size_t inputVecPos=0;
    std::size_t outputVecPos=0;
    
    
    log_debug("Starting iterations..");
    
    for (uint64_t i=0;i<iterations;++i)
    {
        log_debug("Filling input buffer..");
        
        {
            if (i==iterations-1)
            {
                queue->enqueueWriteBuffer(
                    *deviceBuffer_inputs_x,
                    CL_FALSE,
                    0,
                    numEntriesInLastIteration*sizeof(float),
                    &((*xValues)[inputVecPos]),
                    NULL, NULL);
                queue->enqueueWriteBuffer(
                    *deviceBuffer_inputs_y,
                    CL_FALSE,
                    0,
                    numEntriesInLastIteration*sizeof(float),
                    &((*yValues)[inputVecPos]),
                    NULL, NULL);
                queue->enqueueWriteBuffer(
                    *deviceBuffer_inputs_z,
                    CL_FALSE,
                    0,
                    numEntriesInLastIteration*sizeof(float),
                    &((*zValues)[inputVecPos]),
                    NULL, NULL);

                inputVecPos += numEntriesInLastIteration;
            } else {
                queue->enqueueWriteBuffer(
                    *deviceBuffer_inputs_x,
                    CL_FALSE,
                    0,
                    workItemsPerIteration*sizeof(float),
                    &((*xValues)[inputVecPos]),
                    NULL, NULL);
                queue->enqueueWriteBuffer(
                    *deviceBuffer_inputs_y,
                    CL_FALSE,
                    0,
                    workItemsPerIteration*sizeof(float),
                    &((*yValues)[inputVecPos]),
                    NULL, NULL);
                queue->enqueueWriteBuffer(
                    *deviceBuffer_inputs_z,
                    CL_FALSE,
                    0,
                    workItemsPerIteration*sizeof(float),
                    &((*zValues)[inputVecPos]),
                    NULL, NULL);

                inputVecPos += workItemsPerIteration;
            }
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
            queue->flush();
            queue->finish();
        } catch (cl::Error &err) {
            log_fatal("OpenCL ERROR (running kernel): %s (%i)", err.what(), err.err());
        }
        
        log_debug("kernel finished!");
        
        log_debug("Reading output buffer..");
        {
            if (i==iterations-1)
            {
                queue->enqueueReadBuffer(
                    *deviceBuffer_results,
                    CL_FALSE,
                    0,
                    numEntriesInLastIteration*sizeof(float),
                    &((*results)[outputVecPos]),
                    NULL, NULL);
                outputVecPos += numEntriesInLastIteration;
            } else {
                queue->enqueueReadBuffer(
                    *deviceBuffer_results,
                    CL_FALSE,
                    0,
                    workItemsPerIteration*sizeof(float),
                    &((*results)[outputVecPos]),
                    NULL, NULL);
                outputVecPos += workItemsPerIteration;
            }
        }
        log_debug("Output buffer read.");

    }
    
    log_debug("iterations complete.");
    
    queue->flush();
    queue->finish();

    return results;
}
