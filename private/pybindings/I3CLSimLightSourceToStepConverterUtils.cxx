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
 * @file I3CLSimLightSourceToStepConverterUtils.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

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
