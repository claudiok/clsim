//
//   Copyright (c) 2011  Claudio Kopper
//   
//   $Id$
//
//   This file is part of the IceTray module clsim.
//
//   clsim is free software; you can redistribute it and/or modify
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

    shared_ptr<cl::Buffer> deviceBuffer_results;
    shared_ptr<cl::Buffer> deviceBuffer_inputs;

    I3CLSimFunctionConstPtr wlenDependentValue_;
    
    SET_LOGGER("I3CLSimFunctionTester");
};



#endif //I3CLSIMWLENDEPENDENTVALUETESTER_H_INCLUDED
