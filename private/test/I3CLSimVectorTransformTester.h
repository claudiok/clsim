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
 * @file I3CLSimVectorTransformTester.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSimVectorTransformTESTER_H_INCLUDED
#define I3CLSimVectorTransformTESTER_H_INCLUDED

#include "test/I3CLSimTesterBase.h"

#include "dataclasses/I3Vector.h"
#include "clsim/function/I3CLSimVectorTransform.h"

#include <vector>

class I3CLSimVectorTransformTester : public I3CLSimTesterBase
{
public:
    I3CLSimVectorTransformTester(const I3CLSimOpenCLDevice &device,
                                    uint64_t workgroupSize_,
                                    uint64_t workItemsPerIteration_,
                                    I3CLSimVectorTransformConstPtr theTransform);

    // evaluates the function using an OpenCL kernel
    std::vector<I3VectorFloatPtr> EvaluateFunction(
        I3VectorFloatConstPtr xValues,
        I3VectorFloatConstPtr yValues,
        I3VectorFloatConstPtr zValues);

    // evaluates the function using compiled code (i.e. using the 
    // I3CLSimVectorTransform object)
    std::vector<I3VectorFloatPtr> EvaluateReferenceFunction(
        I3VectorFloatConstPtr xValues,
        I3VectorFloatConstPtr yValues,
        I3VectorFloatConstPtr zValues);

private:
    void FillSource(std::vector<std::string> &source,
                    I3CLSimVectorTransformConstPtr theTransform);

    void InitBuffers();

    boost::shared_ptr<cl::Buffer> deviceBuffer_results_x;
    boost::shared_ptr<cl::Buffer> deviceBuffer_results_y;
    boost::shared_ptr<cl::Buffer> deviceBuffer_results_z;
    boost::shared_ptr<cl::Buffer> deviceBuffer_inputs_x;
    boost::shared_ptr<cl::Buffer> deviceBuffer_inputs_y;
    boost::shared_ptr<cl::Buffer> deviceBuffer_inputs_z;

    I3CLSimVectorTransformConstPtr theTransform_;
    
    SET_LOGGER("I3CLSimVectorTransformTester");
};



#endif //I3CLSimVectorTransformTESTER_H_INCLUDED
