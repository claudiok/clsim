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
 * @file I3CLSimRandomValueConstant.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <icetray/serialization.h>
#include <icetray/I3Logging.h>
#include <clsim/random_value/I3CLSimRandomValueConstant.h>
#include <cmath>

#include "clsim/I3CLSimHelperToFloatString.h"
using namespace I3CLSimHelper;

I3CLSimRandomValueConstant::I3CLSimRandomValueConstant()
:
value_(NAN)
{ 
}

I3CLSimRandomValueConstant::I3CLSimRandomValueConstant(double value)
:
value_(value)
{ 
}

I3CLSimRandomValueConstant::~I3CLSimRandomValueConstant() 
{ 

}

std::size_t I3CLSimRandomValueConstant::NumberOfParameters() const
{
    if (std::isnan(value_)) {
        return 1;
    } else {
        return 0;
    }
}

double I3CLSimRandomValueConstant::SampleFromDistribution(const I3RandomServicePtr &random,
                                                          const std::vector<double> &parameters) const
{
    if (!random) log_fatal("random service is NULL!");
    
    if (std::isnan(value_)) {
        if (parameters.size() != 1) log_fatal("This distribution expects 1 parameter. Got %zu.", parameters.size());
        return parameters[0];
    } else {
        if (parameters.size() != 0) log_fatal("This distribution expects 0 parameters. Got %zu.", parameters.size());
        return value_;
    }
}


std::string I3CLSimRandomValueConstant::GetOpenCLFunction
(const std::string &functionName,
 const std::string &functionArgs,
 const std::string &functionArgsToCall,
 const std::string &uniformRandomCall_co,
 const std::string &uniformRandomCall_oc
 ) const
{
    if (std::isnan(value_)) {
        const std::string functionDecl = std::string("inline float ") + functionName + "(" + functionArgs + ", float value)";
        
        return functionDecl + ";\n\n" + functionDecl + "\n"
        "{\n"
        "    return value;\n"
        "}\n"
        ;
    } else {
        const std::string functionDecl = std::string("inline float ") + functionName + "(" + functionArgs + ")";

        return functionDecl + ";\n\n" + functionDecl + "\n"
        "{\n"
        "    return " + ToFloatString(value_) + ";\n"
        "}\n"
        ;
    }
}

bool I3CLSimRandomValueConstant::CompareTo(const I3CLSimRandomValue &other) const
{
    try
    {
        const I3CLSimRandomValueConstant &other_ = dynamic_cast<const I3CLSimRandomValueConstant &>(other);
        
        if (std::isnan(other_.value_) && std::isnan(value_)) return true;
        if (other_.value_ == value_) return true;
        
        return false;
    }
    catch (const std::bad_cast& e)
    {
        // not of the same type, treat it as non-equal
        return false;
    }
    
}

template <class Archive>
void I3CLSimRandomValueConstant::serialize(Archive &ar, unsigned version)
{
    if (version>i3clsimrandomvalueconstant_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimRandomValueConstant class.",version,i3clsimrandomvalueconstant_version_);

    ar & make_nvp("I3CLSimRandomValue", base_object<I3CLSimRandomValue>(*this));
    ar & make_nvp("value", value_);
}     


I3_SERIALIZABLE(I3CLSimRandomValueConstant);
