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
 * @file I3CLSimRandomValueFixParameter.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <icetray/serialization.h>
#include <icetray/I3Logging.h>
#include <clsim/random_value/I3CLSimRandomValueFixParameter.h>

#include <boost/lexical_cast.hpp>

#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <limits>

#include "clsim/I3CLSimHelperToFloatString.h"
using namespace I3CLSimHelper;

I3CLSimRandomValueFixParameter::
I3CLSimRandomValueFixParameter(const I3CLSimRandomValuePtr &randomDistUsed,
                               std::size_t parameterIndex,
                               double parameterValue)
:
randomDistUsed_(randomDistUsed),
parameterIndex_(parameterIndex),
parameterValue_(parameterValue)
{
    if (!randomDistUsed_) {
        log_fatal("The \"randomDistUsed\" parameter must not be NULL!");
    }
    
    if (parameterIndex_ >= randomDistUsed_->NumberOfParameters()) {
        log_fatal("The input distibution only has %zu parameters, cannot fix parameter index %zu!",
                  randomDistUsed_->NumberOfParameters(), parameterIndex_);
    }
}

I3CLSimRandomValueFixParameter::~I3CLSimRandomValueFixParameter() 
{ 
}

I3CLSimRandomValueFixParameter::I3CLSimRandomValueFixParameter() {;}

std::size_t I3CLSimRandomValueFixParameter::NumberOfParameters() const
{
    return randomDistUsed_->NumberOfParameters()-1;
}

double I3CLSimRandomValueFixParameter::SampleFromDistribution(const I3RandomServicePtr &random,
                                                              const std::vector<double> &parameters) const
{
    if (!random) log_fatal("random service is NULL!");
    if (parameters.size() != randomDistUsed_->NumberOfParameters()-1)
        log_fatal("This distribution expects %zu parameters. Got %zu.",
                  randomDistUsed_->NumberOfParameters()-1,
                  parameters.size());
    
    // build the new parameter list (with our new value inserted)
    std::vector<double> new_parameters(randomDistUsed_->NumberOfParameters()); 
    for (std::size_t i=0;i<parameterIndex_;++i){
        new_parameters[i] = parameters[i];
    }
    new_parameters[parameterIndex_] = parameterValue_;
    for (std::size_t i=parameterIndex_+1;i<new_parameters.size();++i){
        new_parameters[i] = parameters[i-1];
    }
    
    // call the original distribution
    return randomDistUsed_->SampleFromDistribution(random, new_parameters);
}


std::string I3CLSimRandomValueFixParameter::GetOpenCLFunction
(const std::string &functionName,
 const std::string &functionArgs,
 const std::string &functionArgsToCall,
 const std::string &uniformRandomCall_co,
 const std::string &uniformRandomCall_oc
 ) const
{
    // build the argument and parameter lists
    std::string parametersInDecl, parametersToCall;
    for (std::size_t i=0;i<parameterIndex_;++i){
        parametersInDecl += ", float p" + boost::lexical_cast<std::string>(i);
        parametersToCall += ", p" + boost::lexical_cast<std::string>(i);
    }
    parametersToCall += ", " + ToFloatString(parameterValue_);
    for (std::size_t i=parameterIndex_+1;i<randomDistUsed_->NumberOfParameters();++i){
        parametersInDecl += ", float p" + boost::lexical_cast<std::string>(i);
        parametersToCall += ", p" + boost::lexical_cast<std::string>(i);
    }

    const std::string functionDecl = std::string("inline float ") + functionName + "(" + functionArgs + parametersInDecl +")";
    const std::string functionNameUsed = std::string("") + functionName + "_used";

    const std::string usedFunctionDecl =
    randomDistUsed_->GetOpenCLFunction
    (
     functionNameUsed,
     functionArgs,
     functionArgsToCall,
     uniformRandomCall_co,
     uniformRandomCall_oc
    );
    
    return usedFunctionDecl + "\n\n" + functionDecl + ";\n\n" + functionDecl + "\n"
    "{\n"
    "    return " + functionNameUsed + "(" + functionArgsToCall + parametersToCall + ");\n"
    "}\n"
    ;
}

bool I3CLSimRandomValueFixParameter::CompareTo(const I3CLSimRandomValue &other) const
{
    try
    {
        const I3CLSimRandomValueFixParameter &other_ = dynamic_cast<const I3CLSimRandomValueFixParameter &>(other);

        if (other_.parameterIndex_ != parameterIndex_) return false;
        if (!(std::isnan(other_.parameterValue_) && (std::isnan(parameterValue_)))) {
            if (other_.parameterValue_ != parameterValue_) return false;
        }
        if (!(*(other_.randomDistUsed_) == *randomDistUsed_)) return false;
        
        return true;
    }
    catch (const std::bad_cast& e)
    {
        // not of the same type, treat it as non-equal
        return false;
    }
    
}

template <class Archive>
void I3CLSimRandomValueFixParameter::serialize(Archive &ar, unsigned version)
{
    if (version>i3clsimrandomvaluefixparameter_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimRandomValueFixParameter class.",
                  version,
                  i3clsimrandomvaluefixparameter_version_);

    ar & make_nvp("I3CLSimRandomValue", base_object<I3CLSimRandomValue>(*this));
    ar & make_nvp("parameterIndex", parameterIndex_);
    ar & make_nvp("parameterValue", parameterValue_);
    ar & make_nvp("randomDistUsed", randomDistUsed_);
}     


I3_SERIALIZABLE(I3CLSimRandomValueFixParameter);
