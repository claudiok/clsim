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

#ifndef I3CLSIMRANDOMDISTRIBUTIONTESTER_H_INCLUDED
#define I3CLSIMRANDOMDISTRIBUTIONTESTER_H_INCLUDED

#include "test/I3CLSimTesterBase.h"

#include "phys-services/I3RandomService.h"

#include "dataclasses/I3Vector.h"
#include "clsim/I3CLSimRandomValue.h"

class I3CLSimRandomDistributionTester : public I3CLSimTesterBase
{
public:
    I3CLSimRandomDistributionTester(const I3CLSimOpenCLDevice &device,
                                    uint64_t workgroupSize_,
                                    uint64_t workItemsPerIteration_,
                                    I3RandomServicePtr randomService,
                                    I3CLSimRandomValueConstPtr randomDistribution);

    I3VectorFloatPtr GenerateRandomNumbers(uint64_t iterations);
    
private:
    void FillSource(std::vector<std::string> &source,
                    I3CLSimRandomValueConstPtr randomDistribution);
    
    void InitBuffers(I3RandomServicePtr randomService);
    
    shared_ptr<cl::Buffer> deviceBuffer_MWC_RNG_x;
    shared_ptr<cl::Buffer> deviceBuffer_MWC_RNG_a;

    shared_ptr<cl::Buffer> deviceBuffer_results;

};



#endif //I3CLSIMRANDOMDISTRIBUTIONTESTER_H_INCLUDED
