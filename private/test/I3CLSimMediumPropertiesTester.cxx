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

#define __STDC_FORMAT_MACROS 
#include <inttypes.h>

#include "test/I3CLSimMediumPropertiesTester.h"

#include <string>

#include "opencl/I3CLSimHelperLoadProgramSource.h"
#include "opencl/I3CLSimHelperGenerateMediumPropertiesSource.h"
#include "opencl/mwcrng_init.h"

I3CLSimMediumPropertiesTester::I3CLSimMediumPropertiesTester
(const std::pair<std::string, std::string> &platformAndDeviceName,
 uint64_t workgroupSize_,
 uint64_t workItemsPerIteration_,
 bool useNativeMath,
 I3CLSimMediumPropertiesConstPtr mediumProperties,
 I3RandomServicePtr randomService)
:
I3CLSimTesterBase(),
mediumProperties_(mediumProperties),
randomService_(randomService)
{
    std::vector<std::string> source;
    FillSource(source, mediumProperties);
    
    const bool hasDispersion = mediumProperties->GetPhaseRefractiveIndices()[0]->HasDerivative();

    DoSetup(platformAndDeviceName,
            useNativeMath,
            workgroupSize_,
            workItemsPerIteration_,
            source,
            (!hasDispersion)?"-DNO_DISPERSION ":"");

    InitBuffers(randomService);
}

I3CLSimMediumPropertiesTester::I3CLSimMediumPropertiesTester
(boost::python::tuple platformAndDeviceName,
 uint64_t workgroupSize_,
 uint64_t workItemsPerIteration_,
 bool useNativeMath,
 I3CLSimMediumPropertiesConstPtr mediumProperties,
 I3RandomServicePtr randomService)
:
I3CLSimTesterBase(),
mediumProperties_(mediumProperties),
randomService_(randomService)
{
    std::string platformName = boost::python::extract<std::string>(platformAndDeviceName[0]);
    std::string deviceName = boost::python::extract<std::string>(platformAndDeviceName[1]);
    
    std::pair<std::string, std::string> argument(platformName, deviceName);
    
    std::vector<std::string> source;
    FillSource(source, mediumProperties);

    const bool hasDispersion = mediumProperties->GetPhaseRefractiveIndices()[0]->HasDerivative();
    
    DoSetup(argument,
            useNativeMath,
            workgroupSize_,
            workItemsPerIteration_,
            source,
            (!hasDispersion)?"-DNO_DISPERSION ":"");
    
    InitBuffers(randomService);
}

void I3CLSimMediumPropertiesTester::FillSource(std::vector<std::string> &source,
                                               I3CLSimMediumPropertiesConstPtr mediumProperties)
{
    source.clear();
    
    // load program source from files
    const std::string I3_SRC(getenv("I3_SRC"));
    const std::string kernelBaseDir = I3_SRC+"/clsim/resources/kernels";
    
    std::string mwcrngSource = I3CLSimHelper::LoadProgramSource(kernelBaseDir+"/mwcrng_kernel.cl");

    std::string mediumPropertiesSource = I3CLSimHelper::GenerateMediumPropertiesSource(*mediumProperties);

    std::string testKernelSource = I3CLSimHelper::LoadProgramSource(kernelBaseDir+"/medium_properties_test_kernel.cl");
    
    // collect the program sources
    source.push_back(mwcrngSource);
    source.push_back(mediumPropertiesSource);
    source.push_back(testKernelSource);
}

void I3CLSimMediumPropertiesTester::InitBuffers(I3RandomServicePtr randomService)
{
    if (randomService) {
        // set up rng
        log_debug("Setting up RNG for %" PRIu64 " workitems. (requested workgroupSize=%" PRIu64 ", maxWorkgroupSize=%" PRIu64 ")", workItemsPerIteration, workgroupSize, GetMaxWorkgroupSize());

        std::vector<uint64_t> MWC_RNG_x;
        std::vector<uint32_t> MWC_RNG_a;
    
        MWC_RNG_x.resize(workItemsPerIteration);
        MWC_RNG_a.resize(workItemsPerIteration);
    
        const std::string I3_SRC(getenv("I3_SRC"));
        if (init_MWC_RNG(&(MWC_RNG_x[0]), &(MWC_RNG_a[0]), workItemsPerIteration, (I3_SRC+"/clsim/resources/safeprimes_base32.txt").c_str(), randomService)!=0) 
            throw std::runtime_error("I3CLSimStepToPhotonConverterOpenCL already initialized!");
        log_debug("RNG is set up..");

        // set up device buffers from existing host buffers
        deviceBuffer_MWC_RNG_x = shared_ptr<cl::Buffer>(new cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, MWC_RNG_x.size() * sizeof(uint64_t), &(MWC_RNG_x[0])));
        deviceBuffer_MWC_RNG_a = shared_ptr<cl::Buffer>(new cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, MWC_RNG_a.size() * sizeof(uint32_t), &(MWC_RNG_a[0])));
    }
    else
    {
        // set up dummy buffers
        deviceBuffer_MWC_RNG_x = shared_ptr<cl::Buffer>(new cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, workItemsPerIteration * sizeof(uint64_t), NULL ));
        deviceBuffer_MWC_RNG_a = shared_ptr<cl::Buffer>(new cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, workItemsPerIteration * sizeof(uint32_t), NULL ));
    }
    
    log_debug("Setting up device buffers.");
    // allocate empty buffers on the device
    deviceBuffer_results = shared_ptr<cl::Buffer>(new cl::Buffer(*context, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, workItemsPerIteration*sizeof(float), NULL));
    deviceBuffer_inputs = shared_ptr<cl::Buffer>(new cl::Buffer(*context, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, workItemsPerIteration*sizeof(float), NULL));
    log_debug("Device buffers are set up.");
    
    log_debug("Configuring kernel.");
    {
        kernel->setArg(0, *deviceBuffer_MWC_RNG_x);        // rng state
        kernel->setArg(1, *deviceBuffer_MWC_RNG_a);        // rng state
        kernel->setArg(2, *deviceBuffer_inputs);            // input data
        kernel->setArg(3, *deviceBuffer_results);          // output data
    }
    log_debug("Kernel configured.");
}


I3VectorFloatPtr I3CLSimMediumPropertiesTester::EvaluatePhaseRefIndex(I3VectorFloatConstPtr xValues, uint32_t layer)
{
    return EvaluateIt(xValues, layer, 0);
}

I3VectorFloatPtr I3CLSimMediumPropertiesTester::EvaluateDispersion(I3VectorFloatConstPtr xValues, uint32_t layer)
{
    return EvaluateIt(xValues, layer, 1);
}

I3VectorFloatPtr I3CLSimMediumPropertiesTester::EvaluateGroupVelocity(I3VectorFloatConstPtr xValues, uint32_t layer)
{
    return EvaluateIt(xValues, layer, 2);
}

I3VectorFloatPtr I3CLSimMediumPropertiesTester::EvaluateAbsorptionLength(I3VectorFloatConstPtr xValues, uint32_t layer)
{
    return EvaluateIt(xValues, layer, 3);
}

I3VectorFloatPtr I3CLSimMediumPropertiesTester::EvaluateScatteringLength(I3VectorFloatConstPtr xValues, uint32_t layer)
{
    return EvaluateIt(xValues, layer, 4);
}



I3VectorFloatPtr I3CLSimMediumPropertiesTester::EvaluateIt(I3VectorFloatConstPtr xValues, uint32_t layer, uint32_t mode)
{
    if (!xValues) log_fatal("NULL pointer passed to EvaluateFunction.");

    if (layer >= mediumProperties_->GetLayersNum())
        log_fatal("Medium only has layers 0 to %" PRIu32 ". You requested layer %" PRIu32,
                  mediumProperties_->GetLayersNum()-1, layer);
    
    // tell the kernel whether the function or its derivative is requested
    kernel->setArg(4, static_cast<cl_uint>(layer));
    kernel->setArg(5, static_cast<cl_uint>(mode));

    
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
