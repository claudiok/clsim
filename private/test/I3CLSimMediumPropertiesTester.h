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
 * @file I3CLSimMediumPropertiesTester.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSimMediumPropertiesTester_H_INCLUDED
#define I3CLSimMediumPropertiesTester_H_INCLUDED

#include "test/I3CLSimTesterBase.h"

#include "phys-services/I3RandomService.h"

#include "dataclasses/I3Vector.h"
#include "clsim/I3CLSimMediumProperties.h"

class I3CLSimMediumPropertiesTester : public I3CLSimTesterBase
{
public:
    I3CLSimMediumPropertiesTester(const I3CLSimOpenCLDevice &device,
                                  uint64_t workgroupSize_,
                                  uint64_t workItemsPerIteration_,
                                  I3CLSimMediumPropertiesConstPtr mediumProperties,
                                  I3RandomServicePtr randomService = I3RandomServicePtr());

    // evaluates the function using an OpenCL kernel
    I3VectorFloatPtr EvaluatePhaseRefIndex(I3VectorFloatConstPtr xValues, uint32_t layer);
    I3VectorFloatPtr EvaluateDispersion(I3VectorFloatConstPtr xValues, uint32_t layer);
    I3VectorFloatPtr EvaluateGroupVelocity(I3VectorFloatConstPtr xValues, uint32_t layer);
    I3VectorFloatPtr EvaluateAbsorptionLength(I3VectorFloatConstPtr xValues, uint32_t layer);
    I3VectorFloatPtr EvaluateScatteringLength(I3VectorFloatConstPtr xValues, uint32_t layer);

private:
    I3VectorFloatPtr EvaluateIt(I3VectorFloatConstPtr xValues, uint32_t layer, uint32_t mode);

    
    void FillSource(std::vector<std::string> &source,
                    I3CLSimMediumPropertiesConstPtr wlenDependentValue);

    void InitBuffers(I3RandomServicePtr randomService);

    boost::shared_ptr<cl::Buffer> deviceBuffer_MWC_RNG_x;
    boost::shared_ptr<cl::Buffer> deviceBuffer_MWC_RNG_a;

    boost::shared_ptr<cl::Buffer> deviceBuffer_results;
    boost::shared_ptr<cl::Buffer> deviceBuffer_inputs;

    I3CLSimMediumPropertiesConstPtr mediumProperties_;
    I3RandomServicePtr randomService_;
};



#endif //I3CLSimMediumPropertiesTester_H_INCLUDED
