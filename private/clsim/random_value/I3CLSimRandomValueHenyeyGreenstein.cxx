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
 * @file I3CLSimRandomValueHenyeyGreenstein.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <icetray/serialization.h>
#include <icetray/I3Logging.h>
#include <clsim/random_value/I3CLSimRandomValueHenyeyGreenstein.h>

#include "clsim/I3CLSimHelperToFloatString.h"
using namespace I3CLSimHelper;

I3CLSimRandomValueHenyeyGreenstein::I3CLSimRandomValueHenyeyGreenstein(double meanCosine)
:
meanCosine_(meanCosine)
{ 
    if ((meanCosine_<=0.) || (meanCosine_>=1.)) 
        log_fatal("Invalid mean cosine (%f) in construction of I3CLSimRandomValueHenyeyGreenstein.",
                  meanCosine_);
}

I3CLSimRandomValueHenyeyGreenstein::~I3CLSimRandomValueHenyeyGreenstein() 
{ 

}

I3CLSimRandomValueHenyeyGreenstein::I3CLSimRandomValueHenyeyGreenstein() {;}

std::size_t I3CLSimRandomValueHenyeyGreenstein::NumberOfParameters() const {return 0;}

double I3CLSimRandomValueHenyeyGreenstein::SampleFromDistribution(const I3RandomServicePtr &random,
                                                                  const std::vector<double> &parameters) const
{
    if (!random) log_fatal("random service is NULL!");
    if (parameters.size() != 0) log_fatal("This distribution expects 0 parameters. Got %zu.", parameters.size());

    const double g = meanCosine_;
    const double g2 = meanCosine_*meanCosine_;
    
    // a random number [-1;+1]
    const double s = 2.*(random->Uniform())-1.;
    
    const double ii = ((1. - g2)/(1. + g*s));
    
    return std::min(std::max((1. + g2 - ii*ii) / (2.*g), -1.), 1.);
}

std::string I3CLSimRandomValueHenyeyGreenstein::GetOpenCLFunction
(const std::string &functionName,
 const std::string &functionArgs,
 const std::string &functionArgsToCall,
 const std::string &uniformRandomCall_co,
 const std::string &uniformRandomCall_oc
 ) const
{
    const std::string functionDecl = std::string("float ") + functionName + "(" + functionArgs + ")";
    
    return functionDecl + ";\n\n" + "inline " + functionDecl + "\n"
    "{\n"
    "    const float g = " + ToFloatString(meanCosine_) + ";\n"
    "    const float g2 = " + ToFloatString(meanCosine_*meanCosine_) + ";\n"
    "    \n"
    "    // a random number [-1;+1]\n"
    "    const float s = 2.f*(" + uniformRandomCall_co + ")-1.f;\n"
    "    \n"
    "    const float ii = ((1.f - g2)/(1.f + g*s));\n"
    "    return clamp((1.f + g2 - ii*ii) / (2.f*g), -1.f, 1.f);\n"
    "}\n"
    ;
}

bool I3CLSimRandomValueHenyeyGreenstein::CompareTo(const I3CLSimRandomValue &other) const
{
    try
    {
        const I3CLSimRandomValueHenyeyGreenstein &other_ = dynamic_cast<const I3CLSimRandomValueHenyeyGreenstein &>(other);
        return ((other_.meanCosine_ == meanCosine_));
    }
    catch (const std::bad_cast& e)
    {
        // not of the same type, treat it as non-equal
        return false;
    }
    
}

template <class Archive>
void I3CLSimRandomValueHenyeyGreenstein::serialize(Archive &ar, unsigned version)
{
    if (version>i3clsimrandomvaluehenyeygreenstein_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimRandomValueHenyeyGreenstein class.",version,i3clsimrandomvaluehenyeygreenstein_version_);

    ar & make_nvp("I3CLSimRandomValue", base_object<I3CLSimRandomValue>(*this));
    ar & make_nvp("meanCosine", meanCosine_);
}     


I3_SERIALIZABLE(I3CLSimRandomValueHenyeyGreenstein);
