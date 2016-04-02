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
 * @file I3CLSimRandomValueSimplifiedLiu.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <icetray/serialization.h>
#include <icetray/I3Logging.h>
#include <clsim/random_value/I3CLSimRandomValueSimplifiedLiu.h>

#include "clsim/I3CLSimHelperToFloatString.h"
using namespace I3CLSimHelper;

I3CLSimRandomValueSimplifiedLiu::I3CLSimRandomValueSimplifiedLiu(double meanCosine)
:
meanCosine_(meanCosine)
{ 
    if ((meanCosine_<=0.) || (meanCosine_>=1.)) 
        log_fatal("Invalid mean cosine (%f) in construction of I3CLSimRandomValueSimplifiedLiu.",
                  meanCosine_);
}

I3CLSimRandomValueSimplifiedLiu::~I3CLSimRandomValueSimplifiedLiu() 
{ 

}

I3CLSimRandomValueSimplifiedLiu::I3CLSimRandomValueSimplifiedLiu() {;}

std::size_t I3CLSimRandomValueSimplifiedLiu::NumberOfParameters() const {return 0;}

double I3CLSimRandomValueSimplifiedLiu::SampleFromDistribution(const I3RandomServicePtr &random,
                                                               const std::vector<double> &parameters) const
{
    if (!random) log_fatal("random service is NULL!");
    if (parameters.size() != 0) log_fatal("This distribution expects 0 parameters. Got %zu.", parameters.size());

    const double beta = (1.-meanCosine_)/(1.+meanCosine_);

    return std::min(std::max(2.*std::pow((random->Uniform()),beta)-1., -1.), 1.);
}


std::string I3CLSimRandomValueSimplifiedLiu::GetOpenCLFunction
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
    "    //const float g = " + ToFloatString(meanCosine_) + ";\n"
    "    //const float beta = (1.f-g)/(1.f+g);\n"
    "    const float beta = " + ToFloatString((1.-meanCosine_)/(1.+meanCosine_)) + ";\n"
    "    \n"
    "#ifdef USE_NATIVE_MATH\n"
    "    return clamp(2.f * native_powr((" + uniformRandomCall_co + "), beta) - 1.f, -1.f, 1.f);\n"
    "#else\n"
    "    return clamp(2.f * powr((" + uniformRandomCall_co + "), beta) - 1.f, -1.f, 1.f);\n"
    "#endif\n"
    "}\n"
    ;
}

bool I3CLSimRandomValueSimplifiedLiu::CompareTo(const I3CLSimRandomValue &other) const
{
    try
    {
        const I3CLSimRandomValueSimplifiedLiu &other_ = dynamic_cast<const I3CLSimRandomValueSimplifiedLiu &>(other);
        return ((other_.meanCosine_ == meanCosine_));
    }
    catch (const std::bad_cast& e)
    {
        // not of the same type, treat it as non-equal
        return false;
    }
    
}

template <class Archive>
void I3CLSimRandomValueSimplifiedLiu::serialize(Archive &ar, unsigned version)
{
    if (version>i3clsimrandomvaluesimplifiedliu_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimRandomValueSimplifiedLiu class.",version,i3clsimrandomvaluesimplifiedliu_version_);

    ar & make_nvp("I3CLSimRandomValue", base_object<I3CLSimRandomValue>(*this));
    ar & make_nvp("meanCosine", meanCosine_);
}     


I3_SERIALIZABLE(I3CLSimRandomValueSimplifiedLiu);
