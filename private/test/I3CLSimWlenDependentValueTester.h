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

#ifndef I3CLSIMWLENDEPENDENTVALUETESTER_H_INCLUDED
#define I3CLSIMWLENDEPENDENTVALUETESTER_H_INCLUDED

#include "test/I3CLSimTesterBase.h"

#include "dataclasses/I3Vector.h"
#include "clsim/I3CLSimWlenDependentValue.h"

class I3CLSimWlenDependentValueTester : public I3CLSimTesterBase
{
public:
    // version to be called from c++
    I3CLSimWlenDependentValueTester(const std::pair<std::string, std::string> &platformAndDeviceName,
                                    uint64_t workgroupSize_,
                                    uint64_t workItemsPerIteration_,
                                    bool useNativeMath,
                                    I3CLSimWlenDependentValueConstPtr wlenDependentValue);

    // version to be called from boost::python
    I3CLSimWlenDependentValueTester(boost::python::tuple platformAndDeviceName,
                                    uint64_t workgroupSize_,
                                    uint64_t workItemsPerIteration_,
                                    bool useNativeMath,
                                    I3CLSimWlenDependentValueConstPtr wlenDependentValue);

    // evaluates the function using an OpenCL kernel
    I3VectorFloatPtr EvaluateFunction(I3VectorFloatConstPtr xValues);
    I3VectorFloatPtr EvaluateDerivative(I3VectorFloatConstPtr xValues);

    // evaluates the function using compiled code (i.e. using the 
    // I3CLSimWlenDependentValue object)
    I3VectorFloatPtr EvaluateReferenceFunction(I3VectorFloatConstPtr xValues);
    I3VectorFloatPtr EvaluateReferenceDerivative(I3VectorFloatConstPtr xValues);

private:
    I3VectorFloatPtr EvaluateIt(I3VectorFloatConstPtr xValues, bool derivative);

    
    void FillSource(std::vector<std::string> &source,
                    I3CLSimWlenDependentValueConstPtr wlenDependentValue);

    void InitBuffers();

    shared_ptr<cl::Buffer> deviceBuffer_results;
    shared_ptr<cl::Buffer> deviceBuffer_inputs;

    I3CLSimWlenDependentValueConstPtr wlenDependentValue_;
};



#endif //I3CLSIMWLENDEPENDENTVALUETESTER_H_INCLUDED
