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
 * @file I3CLSimRandomValueUniform.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <icetray/serialization.h>
#include <icetray/I3Logging.h>
#include <clsim/random_value/I3CLSimRandomValueUniform.h>
#include <cmath>

#include "clsim/I3CLSimHelperToFloatString.h"
using namespace I3CLSimHelper;

I3CLSimRandomValueUniform::I3CLSimRandomValueUniform()
:
from_(NAN),
to_(NAN)
{ 
}

I3CLSimRandomValueUniform::I3CLSimRandomValueUniform(double from, double to)
:
from_(from),
to_(to)
{ 
}

I3CLSimRandomValueUniform::~I3CLSimRandomValueUniform() 
{ 

}

std::size_t I3CLSimRandomValueUniform::NumberOfParameters() const
{
    if (std::isnan(from_)) {
        if (std::isnan(to_))
        {
            return 2;
        }
        else
        {
            return 1;
        }
    } else {
        if (std::isnan(to_))
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }
}

double I3CLSimRandomValueUniform::SampleFromDistribution(const I3RandomServicePtr &random,
                                                         const std::vector<double> &parameters) const
{
    if (!random) log_fatal("random service is NULL!");
    if (parameters.size() != NumberOfParameters())
        log_fatal("This distribution expects %zu parameters. Got %zu.",
                  NumberOfParameters(),
                  parameters.size());

    double from=from_, to=to_;
    
    {
        std::size_t vec_index=0;
        if (std::isnan(from)) {
            from = parameters[vec_index];
            ++vec_index;
        }

        if (std::isnan(to)) {
            to = parameters[vec_index];
            ++vec_index;
        }
    }
    
    const double width = to-from;
    
    return random->Uniform()*width + from;
}


std::string I3CLSimRandomValueUniform::GetOpenCLFunction
(const std::string &functionName,
 const std::string &functionArgs,
 const std::string &functionArgsToCall,
 const std::string &uniformRandomCall_co,
 const std::string &uniformRandomCall_oc
 ) const
{
    std::string parametersInDecl;
    std::string fromSubstitution;
    std::string toSubstitution;
    
    if (std::isnan(from_)) {
        parametersInDecl += ", float from";
        fromSubstitution = "from";
    } else {
        fromSubstitution = ToFloatString(from_);
    }

    if (std::isnan(to_)) {
        parametersInDecl += ", float to";
        fromSubstitution = "to";
    } else {
        fromSubstitution = ToFloatString(to_);
    }

    const std::string functionDecl = std::string("inline float ") + functionName + "(" + functionArgs + parametersInDecl + ")";

    return functionDecl + ";\n\n" + functionDecl + "\n"
    "{\n"
    "    const float width = " + toSubstitution + " - " + fromSubstitution + ";\n"
    "    \n"
    "    return " + uniformRandomCall_co + "*width + " + fromSubstitution + ";\n"
    "}\n"
    ;
}

bool I3CLSimRandomValueUniform::CompareTo(const I3CLSimRandomValue &other) const
{
    try
    {
        const I3CLSimRandomValueUniform &other_ = dynamic_cast<const I3CLSimRandomValueUniform &>(other);
        if (!(std::isnan(other_.from_) && std::isnan(from_))) {
            if (other_.from_ != from_) return false;
        }

        if (!(std::isnan(other_.to_) && std::isnan(to_))) {
            if (other_.to_ != to_) return false;
        }

        return true;
    }
    catch (const std::bad_cast& e)
    {
        // not of the same type, treat it as non-equal
        return false;
    }
    
}

template <class Archive>
void I3CLSimRandomValueUniform::serialize(Archive &ar, unsigned version)
{
    if (version>i3clsimrandomvalueuniform_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimRandomValueUniform class.",version,i3clsimrandomvalueuniform_version_);

    ar & make_nvp("I3CLSimRandomValue", base_object<I3CLSimRandomValue>(*this));
    ar & make_nvp("from", from_);
    ar & make_nvp("to", to_);
}     


I3_SERIALIZABLE(I3CLSimRandomValueUniform);
