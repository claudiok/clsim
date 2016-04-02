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
 * @file I3CLSimModuleHelper.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <sstream>

#include <clsim/I3CLSimModuleHelper.h>
#include <clsim/random_value/I3CLSimRandomValueInterpolatedDistribution.h>
#include <clsim/random_value/I3CLSimRandomValueWlenCherenkovNoDispersion.h>

#include <boost/preprocessor/seq.hpp>

#include <boost/foreach.hpp>

using namespace boost::python;
namespace bp = boost::python;


void register_I3ModuleHelper()
{
    // this can be used for testing purposes
    bp::def("makeCherenkovWavelengthGenerator", &I3CLSimModuleHelper::makeCherenkovWavelengthGenerator);
    bp::def("makeWavelengthGenerator", &I3CLSimModuleHelper::makeWavelengthGenerator);
    bp::def("initializeOpenCL", &I3CLSimModuleHelper::initializeOpenCL,
        (bp::arg("openCLDevice"), "randomService", "geometry", "mediumProperties",
	"wavelengthGenerationBias", "wavelengthGenerators",
	bp::arg("enableDoubleBuffering")=false, bp::arg("doublePrecision")=false,
	bp::arg("stopDetectedPhotons")=true, bp::arg("saveAllPhotons")=false,
	bp::arg("saveAllPhotonsPrescale")=0.01, bp::arg("fixedNumberOfAbsorptionLengths")=NAN,
	bp::arg("pancakeFactor")=1., bp::arg("photonHistoryEntries")=0,
	bp::arg("limitWorkgroupSize")=0));
    
}
