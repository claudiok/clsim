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
 * @file I3CLSimVectorTransformTester.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>

#include "test/I3CLSimVectorTransformTester.h"

#include <string>

#include "opencl/I3CLSimHelperLoadProgramSource.h"

I3CLSimVectorTransformTester::I3CLSimVectorTransformTester
(const I3CLSimOpenCLDevice &device,
 uint64_t workgroupSize_,
 uint64_t workItemsPerIteration_,
 I3CLSimVectorTransformConstPtr theTransform):
I3CLSimTesterBase(),
theTransform_(theTransform)
{
    std::vector<std::string> source;
    FillSource(source, theTransform);
    
    DoSetup(device,
            workgroupSize_,
            workItemsPerIteration_,
            source);
    
    InitBuffers();
}

void I3CLSimVectorTransformTester::FillSource(
    std::vector<std::string> &source,
    I3CLSimVectorTransformConstPtr theTransform)
{
    source.clear();
    
    // load program source from files
    const std::string I3_BUILD(getenv("I3_BUILD"));
    const std::string kernelBaseDir = I3_BUILD+"/clsim/resources/kernels";
    
    std::string functionSource = theTransform->GetOpenCLFunction("evaluateVectorTransform");
    
    std::string testKernelHeader = I3CLSimHelper::LoadProgramSource(kernelBaseDir+"/vector_transform_test_kernel.h.cl");
    std::string testKernelSource = I3CLSimHelper::LoadProgramSource(kernelBaseDir+"/vector_transform_test_kernel.c.cl");
    
    // collect the program sources
    source.push_back(testKernelHeader);
    source.push_back(functionSource);
    source.push_back(testKernelSource);
}

void I3CLSimVectorTransformTester::InitBuffers()
{
    log_debug("Setting up device buffers.");
    // allocate empty buffers on the device
    deviceBuffer_results_x = boost::shared_ptr<cl::Buffer>(new cl::Buffer(*context, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, workItemsPerIteration*sizeof(float), NULL));
    deviceBuffer_results_y = boost::shared_ptr<cl::Buffer>(new cl::Buffer(*context, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, workItemsPerIteration*sizeof(float), NULL));
    deviceBuffer_results_z = boost::shared_ptr<cl::Buffer>(new cl::Buffer(*context, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, workItemsPerIteration*sizeof(float), NULL));
    deviceBuffer_inputs_x  = boost::shared_ptr<cl::Buffer>(new cl::Buffer(*context, CL_MEM_READ_ONLY  | CL_MEM_ALLOC_HOST_PTR, workItemsPerIteration*sizeof(float), NULL));
    deviceBuffer_inputs_y  = boost::shared_ptr<cl::Buffer>(new cl::Buffer(*context, CL_MEM_READ_ONLY  | CL_MEM_ALLOC_HOST_PTR, workItemsPerIteration*sizeof(float), NULL));
    deviceBuffer_inputs_z  = boost::shared_ptr<cl::Buffer>(new cl::Buffer(*context, CL_MEM_READ_ONLY  | CL_MEM_ALLOC_HOST_PTR, workItemsPerIteration*sizeof(float), NULL));
    log_debug("Device buffers are set up.");
    
    log_debug("Configuring kernel.");
    {
        kernel->setArg(0, *deviceBuffer_inputs_x);        // input data
        kernel->setArg(1, *deviceBuffer_inputs_y);        // input data
        kernel->setArg(2, *deviceBuffer_inputs_z);        // input data
        kernel->setArg(3, *deviceBuffer_results_x);       // output data
        kernel->setArg(4, *deviceBuffer_results_y);       // output data
        kernel->setArg(5, *deviceBuffer_results_z);       // output data
    }
    log_debug("Kernel configured.");
}


std::vector<I3VectorFloatPtr> I3CLSimVectorTransformTester::EvaluateReferenceFunction(
        I3VectorFloatConstPtr xValues,
        I3VectorFloatConstPtr yValues,
        I3VectorFloatConstPtr zValues)
{
    if (!xValues) log_fatal("NULL pointer passed to EvaluateReferenceFunction. (xValues)");
    if (!yValues) log_fatal("NULL pointer passed to EvaluateReferenceFunction. (yValues)");
    if (!zValues) log_fatal("NULL pointer passed to EvaluateReferenceFunction. (zValues)");
    if (xValues->size() != yValues->size()) log_fatal("the input vectors have to have the same size!");
    if (xValues->size() != zValues->size()) log_fatal("the input vectors have to have the same size!");

    if (!theTransform_) log_fatal("Internal error: theTransform_ is NULL");
    if (!theTransform_->HasNativeImplementation()) log_fatal("No native/reference implementation available!");
    
    // allocate the output vector
    std::vector<I3VectorFloatPtr> result;
    for (std::size_t i=0;i<3;++i) {
        result.push_back(I3VectorFloatPtr(new I3VectorFloat(xValues->size(), NAN)));
    }
    
    for (std::size_t i=0;i<xValues->size();++i)
    {
        std::vector<double> oldVec(3,NAN);
        oldVec[0] = (*xValues)[i];
        oldVec[1] = (*yValues)[i];
        oldVec[2] = (*zValues)[i];

        const std::vector<double> newVec = theTransform_->ApplyTransform(oldVec);

        (*(result[0]))[i] = newVec[0];
        (*(result[1]))[i] = newVec[1];
        (*(result[2]))[i] = newVec[2];
    }
    
    return result;
}

std::vector<I3VectorFloatPtr> I3CLSimVectorTransformTester::EvaluateFunction(
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
    std::vector<I3VectorFloatPtr> result;
    for (std::size_t i=0;i<3;++i) {
        result.push_back(I3VectorFloatPtr(new I3VectorFloat(xValues->size(), NAN)));
        // reserve space for output data
    }

    // if nothing to do, exit here
    if (xValues->empty()) return result;

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
                    *deviceBuffer_results_x,
                    CL_FALSE,
                    0,
                    numEntriesInLastIteration*sizeof(float),
                    &((*(result[0]))[outputVecPos]),
                    NULL, NULL);
                queue->enqueueReadBuffer(
                    *deviceBuffer_results_y,
                    CL_FALSE,
                    0,
                    numEntriesInLastIteration*sizeof(float),
                    &((*(result[1]))[outputVecPos]),
                    NULL, NULL);
                queue->enqueueReadBuffer(
                    *deviceBuffer_results_z,
                    CL_FALSE,
                    0,
                    numEntriesInLastIteration*sizeof(float),
                    &((*(result[2]))[outputVecPos]),
                    NULL, NULL);
                outputVecPos += numEntriesInLastIteration;
            } else {
                queue->enqueueReadBuffer(
                    *deviceBuffer_results_x,
                    CL_FALSE,
                    0,
                    workItemsPerIteration*sizeof(float),
                    &((*(result[0]))[outputVecPos]),
                    NULL, NULL);
                queue->enqueueReadBuffer(
                    *deviceBuffer_results_y,
                    CL_FALSE,
                    0,
                    workItemsPerIteration*sizeof(float),
                    &((*(result[1]))[outputVecPos]),
                    NULL, NULL);
                queue->enqueueReadBuffer(
                    *deviceBuffer_results_z,
                    CL_FALSE,
                    0,
                    workItemsPerIteration*sizeof(float),
                    &((*(result[2]))[outputVecPos]),
                    NULL, NULL);
                outputVecPos += workItemsPerIteration;
            }
        }
        log_debug("Output buffer read.");

    }
    
    log_debug("iterations complete.");

    queue->flush();
    queue->finish();
    
    return result;
}
