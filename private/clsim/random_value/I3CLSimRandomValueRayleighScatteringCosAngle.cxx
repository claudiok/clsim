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
 * @file I3CLSimRandomValueRayleighScatteringCosAngle.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <icetray/serialization.h>
#include <icetray/I3Logging.h>
#include <clsim/random_value/I3CLSimRandomValueRayleighScatteringCosAngle.h>

#include "clsim/I3CLSimHelperToFloatString.h"
using namespace I3CLSimHelper;

I3CLSimRandomValueRayleighScatteringCosAngle::I3CLSimRandomValueRayleighScatteringCosAngle()
{ 
}

I3CLSimRandomValueRayleighScatteringCosAngle::~I3CLSimRandomValueRayleighScatteringCosAngle() 
{ 

}

std::size_t I3CLSimRandomValueRayleighScatteringCosAngle::NumberOfParameters() const {return 0;}

double I3CLSimRandomValueRayleighScatteringCosAngle::SampleFromDistribution(const I3RandomServicePtr &random,
                                                                            const std::vector<double> &parameters) const
{
    if (!random) log_fatal("random service is NULL!");
    if (parameters.size() != 0) log_fatal("This distribution expects 0 parameters. Got %zu.", parameters.size());

    const double b = 0.835;
    const double p = 1./0.835;

    const double q = (b+3.)*((random->Uniform())-0.5)/b;
    const double d = q*q + p*p*p;

    const double u1 = -q+std::sqrt(d);
    const double u = std::pow((std::abs(u1)),(1./3.)) * ((u1<0.)?-1.:1.);

    const double v1 = -q-std::sqrt(d);
    const double v = std::pow((std::abs(v1)),(1./3.)) * ((v1<0.)?-1.:1.);

    return std::min(std::max(u+v, -1.), 1.);
}

std::string I3CLSimRandomValueRayleighScatteringCosAngle::GetOpenCLFunction
(const std::string &functionName,
 const std::string &functionArgs,
 const std::string &functionArgsToCall,
 const std::string &uniformRandomCall_co,
 const std::string &uniformRandomCall_oc
 ) const
{
    const std::string functionDecl = std::string("inline float ") + functionName + "(" + functionArgs + ")";
    
    return functionDecl + ";\n\n" + functionDecl + "\n"
    "{\n"
    "    const float b = 0.835f;\n"
    "    const float p = 1.f/0.835f;\n"
    "    \n"
    "    const float q = (b+3.f)*((" + uniformRandomCall_co + ")-0.5f)/b;\n"
    "    const float d = q*q + p*p*p;\n"
    "    \n"
    "#ifdef USE_NATIVE_MATH\n"
    "    const float u1 = -q+native_sqrt(d);\n"
    "    const float u = native_powr((fabs(u1)),(1.f/3.f)) * sign(u1);\n"
    "    \n"
    "    const float v1 = -q-native_sqrt(d);\n"
    "    const float v = native_powr((fabs(v1)),(1.f/3.f)) * sign(v1);\n"
    "#else\n"
    "    const float u1 = -q+sqrt(d);\n"
    "    const float u = powr((fabs(u1)),(1.f/3.f)) * sign(u1);\n"
    "    \n"
    "    const float v1 = -q-sqrt(d);\n"
    "    const float v = powr((fabs(v1)),(1.f/3.f)) * sign(v1);\n"
    "#endif\n"
    "    \n"
    "    return clamp(u+v, -1.f, 1.f);\n"
    "}\n"
    ;
}

bool I3CLSimRandomValueRayleighScatteringCosAngle::CompareTo(const I3CLSimRandomValue &other) const
{
    try
    {
        //const I3CLSimRandomValueRayleighScatteringCosAngle &other_ = dynamic_cast<const I3CLSimRandomValueRayleighScatteringCosAngle &>(other);
        return true; // does not have internal state
    }
    catch (const std::bad_cast& e)
    {
        // not of the same type, treat it as non-equal
        return false;
    }
    
}

template <class Archive>
void I3CLSimRandomValueRayleighScatteringCosAngle::serialize(Archive &ar, unsigned version)
{
    if (version>i3clsimrandomvaluerayleighscatteringcosangle_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimRandomValueRayleighScatteringCosAngle class.",version,i3clsimrandomvaluerayleighscatteringcosangle_version_);

    ar & make_nvp("I3CLSimRandomValue", base_object<I3CLSimRandomValue>(*this));
}     


I3_SERIALIZABLE(I3CLSimRandomValueRayleighScatteringCosAngle);
