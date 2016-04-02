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
 * @file I3CLSimScalarFieldTester.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMSCALARFIELDTESTER_H_INCLUDED
#define I3CLSIMSCALARFIELDTESTER_H_INCLUDED

#include "test/I3CLSimTesterBase.h"

#include "dataclasses/I3Vector.h"
#include "clsim/function/I3CLSimScalarField.h"

class I3CLSimScalarFieldTester : public I3CLSimTesterBase
{
public:
    I3CLSimScalarFieldTester(const I3CLSimOpenCLDevice &device,
                                    uint64_t workgroupSize_,
                                    uint64_t workItemsPerIteration_,
                                    I3CLSimScalarFieldConstPtr theField);

    // evaluates the function using an OpenCL kernel
    I3VectorFloatPtr EvaluateFunction(
        I3VectorFloatConstPtr xValues,
        I3VectorFloatConstPtr yValues,
        I3VectorFloatConstPtr zValues);

    // evaluates the function using compiled code (i.e. using the 
    // I3CLSimScalarField object)
    I3VectorFloatPtr EvaluateReferenceFunction(
        I3VectorFloatConstPtr xValues,
        I3VectorFloatConstPtr yValues,
        I3VectorFloatConstPtr zValues);

private:
    void FillSource(std::vector<std::string> &source,
                    I3CLSimScalarFieldConstPtr theField);

    void InitBuffers();

    boost::shared_ptr<cl::Buffer> deviceBuffer_results;
    boost::shared_ptr<cl::Buffer> deviceBuffer_inputs_x;
    boost::shared_ptr<cl::Buffer> deviceBuffer_inputs_y;
    boost::shared_ptr<cl::Buffer> deviceBuffer_inputs_z;

    I3CLSimScalarFieldConstPtr theField_;
    
    SET_LOGGER("I3CLSimScalarFieldTester");
};



#endif //I3CLSIMSCALARFIELDTESTER_H_INCLUDED
