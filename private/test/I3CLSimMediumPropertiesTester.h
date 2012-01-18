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

    shared_ptr<cl::Buffer> deviceBuffer_MWC_RNG_x;
    shared_ptr<cl::Buffer> deviceBuffer_MWC_RNG_a;

    shared_ptr<cl::Buffer> deviceBuffer_results;
    shared_ptr<cl::Buffer> deviceBuffer_inputs;

    I3CLSimMediumPropertiesConstPtr mediumProperties_;
    I3RandomServicePtr randomService_;
};



#endif //I3CLSimMediumPropertiesTester_H_INCLUDED
