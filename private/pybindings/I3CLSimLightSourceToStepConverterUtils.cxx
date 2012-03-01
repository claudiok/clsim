//
//   Copyright (c) 2011  Claudio Kopper
//   
//   $Id$
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

#include <sstream>

#include <clsim/I3CLSimLightSourceToStepConverterUtils.h>

#include <boost/preprocessor/seq.hpp>

#include <boost/foreach.hpp>

using namespace boost::python;
namespace bp = boost::python;


void register_I3CLSimLightSourceToStepConverterUtils()
{
    double (* gammaDistributedNumber_smartPtr)(double, I3RandomServicePtr)
    = &I3CLSimLightSourceToStepConverterUtils::gammaDistributedNumber;

    // this can be used for testing purposes
    bp::def("NumberOfPhotonsPerMeter", &I3CLSimLightSourceToStepConverterUtils::NumberOfPhotonsPerMeter);
    bp::def("PhotonNumberCorrectionFactorAfterBias", &I3CLSimLightSourceToStepConverterUtils::PhotonNumberCorrectionFactorAfterBias);
    bp::def("gammaDistributedNumber", gammaDistributedNumber_smartPtr);

    //bp::def("scatterDirectionByAngle", &I3CLSimLightSourceToStepConverterUtils::scatterDirectionByAngle);
    //bp::def("GenerateStep", &I3CLSimLightSourceToStepConverterUtils::GenerateStep);
    //bp::def("GenerateStepForMuon", &I3CLSimLightSourceToStepConverterUtils::GenerateStepForMuon);
    
    
}
