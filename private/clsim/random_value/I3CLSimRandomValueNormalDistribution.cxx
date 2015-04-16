/**
 * Copyright (c) 2012
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
 * @file I3CLSimRandomValueNormalDistribution.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <icetray/serialization.h>
#include <icetray/I3Logging.h>
#include <clsim/random_value/I3CLSimRandomValueNormalDistribution.h>

#include "clsim/I3CLSimHelperToFloatString.h"
using namespace I3CLSimHelper;

//I3CLSimRandomValueNormalDistribution::I3CLSimRandomValueNormalDistribution()
//{ 
//}

I3CLSimRandomValueNormalDistribution::~I3CLSimRandomValueNormalDistribution() 
{ 

}

I3CLSimRandomValueNormalDistribution::I3CLSimRandomValueNormalDistribution() {;}

std::size_t I3CLSimRandomValueNormalDistribution::NumberOfParameters() const {return 2;}

double I3CLSimRandomValueNormalDistribution::SampleFromDistribution(const I3RandomServicePtr &random,
                                                                    const std::vector<double> &parameters) const
{
    if (!random) log_fatal("random service is NULL!");
    if (parameters.size() != 2) log_fatal("This distribution expects 2 parameters. Got %zu.", parameters.size());

    const double mean=parameters[0];
    const double sigma=parameters[1];
    
    // sample a value using the Box-Muller transform (needs two random samples from uniform (0,1])
    return (std::sqrt(-2.*std::log(random->Uniform()))*std::sin(2.*M_PI*random->Uniform()))*sigma + mean;
}


std::string I3CLSimRandomValueNormalDistribution::GetOpenCLFunction
(const std::string &functionName,
 const std::string &functionArgs,
 const std::string &functionArgsToCall,
 const std::string &uniformRandomCall_co,
 const std::string &uniformRandomCall_oc
 ) const
{
    const std::string functionDecl = std::string("inline float ") + functionName + "(" + functionArgs + ", float mean, float sigma)";

    return functionDecl + ";\n\n" + functionDecl + "\n"
    "{\n"
    "    const float rnd1 = " + uniformRandomCall_oc + ";\n"
    "    const float rnd2 = " + uniformRandomCall_oc + ";\n"
    "#ifdef USE_NATIVE_MATH\n"
    "    return (native_sqrt(-2.f*native_log(rnd1))*native_sin(2.f*M_PI_F*rnd2))*sigma + mean;\n"
    "#else\n"
    "    return (sqrt(-2.f*log(rnd1))*sin(2.f*M_PI_F*rnd2))*sigma + mean;\n"
    "#endif\n"
    "}\n"
    ;
}

bool I3CLSimRandomValueNormalDistribution::CompareTo(const I3CLSimRandomValue &other) const
{
    try
    {
        // this would be used to compare data members if this class had any:
        // const I3CLSimRandomValueNormalDistribution &other_ = dynamic_cast<const I3CLSimRandomValueNormalDistribution &>(other);
        return true;
    }
    catch (const std::bad_cast& e)
    {
        // not of the same type, treat it as non-equal
        return false;
    }
    
}

template <class Archive>
void I3CLSimRandomValueNormalDistribution::serialize(Archive &ar, unsigned version)
{
    if (version>i3clsimrandomvaluenormaldistribution_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimRandomValueNormalDistribution class.",version,i3clsimrandomvaluenormaldistribution_version_);

    ar & make_nvp("I3CLSimRandomValue", base_object<I3CLSimRandomValue>(*this));
    //ar & make_nvp("meanCosine", meanCosine_);
}     


I3_SERIALIZABLE(I3CLSimRandomValueNormalDistribution);
