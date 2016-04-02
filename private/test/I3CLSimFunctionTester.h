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
 * @file I3CLSimFunctionTester.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMWLENDEPENDENTVALUETESTER_H_INCLUDED
#define I3CLSIMWLENDEPENDENTVALUETESTER_H_INCLUDED

#include "test/I3CLSimTesterBase.h"

#include "dataclasses/I3Vector.h"
#include "clsim/function/I3CLSimFunction.h"

class I3CLSimFunctionTester : public I3CLSimTesterBase
{
public:
    I3CLSimFunctionTester(const I3CLSimOpenCLDevice &device,
                                    uint64_t workgroupSize_,
                                    uint64_t workItemsPerIteration_,
                                    I3CLSimFunctionConstPtr wlenDependentValue);

    // evaluates the function using an OpenCL kernel
    I3VectorFloatPtr EvaluateFunction(I3VectorFloatConstPtr xValues);
    I3VectorFloatPtr EvaluateDerivative(I3VectorFloatConstPtr xValues);

    // evaluates the function using compiled code (i.e. using the 
    // I3CLSimFunction object)
    I3VectorFloatPtr EvaluateReferenceFunction(I3VectorFloatConstPtr xValues);
    I3VectorFloatPtr EvaluateReferenceDerivative(I3VectorFloatConstPtr xValues);

private:
    I3VectorFloatPtr EvaluateIt(I3VectorFloatConstPtr xValues, bool derivative);

    
    void FillSource(std::vector<std::string> &source,
                    I3CLSimFunctionConstPtr wlenDependentValue);

    void InitBuffers();

    boost::shared_ptr<cl::Buffer> deviceBuffer_results;
    boost::shared_ptr<cl::Buffer> deviceBuffer_inputs;

    I3CLSimFunctionConstPtr wlenDependentValue_;
    
    SET_LOGGER("I3CLSimFunctionTester");
};



#endif //I3CLSIMWLENDEPENDENTVALUETESTER_H_INCLUDED
