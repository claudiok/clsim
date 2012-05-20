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
 * @file I3CLSimRandomValueApplyFunction.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <icetray/serialization.h>
#include <icetray/I3Logging.h>
#include <clsim/random_value/I3CLSimRandomValueApplyFunction.h>

#include <boost/lexical_cast.hpp>

#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <limits>

I3CLSimRandomValueApplyFunction::
I3CLSimRandomValueApplyFunction(const I3CLSimRandomValuePtr &randomDistUsed,
                                const std::string &applyFunctionName)
:
randomDistUsed_(randomDistUsed),
applyFunctionName_(applyFunctionName)
{
    if (!randomDistUsed_) {
        log_fatal("The \"randomDistUsed\" parameter must not be NULL!");
    }
    
    if (applyFunctionName_=="") {
        log_fatal("The \"applyFunctionName\" parameter must be set!");
    }
}

I3CLSimRandomValueApplyFunction::~I3CLSimRandomValueApplyFunction() 
{ 
}

I3CLSimRandomValueApplyFunction::I3CLSimRandomValueApplyFunction() {;}

std::size_t I3CLSimRandomValueApplyFunction::NumberOfParameters() const
{
    return randomDistUsed_->NumberOfParameters();
}

double I3CLSimRandomValueApplyFunction::SampleFromDistribution(const I3RandomServicePtr &random,
                                                               const std::vector<double> &parameters) const
{
    if (!random) log_fatal("random service is NULL!");
    
    const double sampledValue = randomDistUsed_->SampleFromDistribution(random, parameters);
    
    if      (applyFunctionName_=="cos")  return std::cos (sampledValue); 
    else if (applyFunctionName_=="sin")  return std::sin (sampledValue); 
    else if (applyFunctionName_=="tan")  return std::tan (sampledValue); 
    else if (applyFunctionName_=="acos") return std::acos(sampledValue); 
    else if (applyFunctionName_=="asin") return std::asin(sampledValue); 
    else if (applyFunctionName_=="atan") return std::atan(sampledValue); 
    else if (applyFunctionName_=="exp")  return std::exp (sampledValue); 
    else {
        log_fatal("The function \"%s\" is currently not implemented in host code.", applyFunctionName_.c_str());
        return NAN;
    }
}


std::string I3CLSimRandomValueApplyFunction::GetOpenCLFunction
(const std::string &functionName,
 const std::string &functionArgs,
 const std::string &functionArgsToCall,
 const std::string &uniformRandomCall_co,
 const std::string &uniformRandomCall_oc
 ) const
{
    // some simple substitutions
    std::string applyFunctionNameNativeMath;
    if (applyFunctionName_=="cos") applyFunctionNameNativeMath="native_cos";
    else if (applyFunctionName_=="sin") applyFunctionNameNativeMath="native_sin";
    else if (applyFunctionName_=="tan") applyFunctionNameNativeMath="native_tan";
    else if (applyFunctionName_=="acos") applyFunctionNameNativeMath="native_acos";
    else if (applyFunctionName_=="asin") applyFunctionNameNativeMath="native_asin";
    else if (applyFunctionName_=="atan") applyFunctionNameNativeMath="native_atan";
    else if (applyFunctionName_=="exp") applyFunctionNameNativeMath="native_exp";
    else applyFunctionNameNativeMath=applyFunctionName_;

    std::string parametersInDecl, parametersToCall;
    for (std::size_t i=0;i<randomDistUsed_->NumberOfParameters();++i)
    {
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
    "    const float number = " + functionNameUsed + "(" + functionArgsToCall + parametersToCall + ");\n"
    "    \n"
    "#ifdef USE_NATIVE_MATH\n"
    "    return " + applyFunctionNameNativeMath + "(number);\n"
    "#else\n"
    "    return " + applyFunctionName_ + "(number);\n"
    "#endif\n"
    "}\n"
    ;
}

bool I3CLSimRandomValueApplyFunction::CompareTo(const I3CLSimRandomValue &other) const
{
    try
    {
        const I3CLSimRandomValueApplyFunction &other_ = dynamic_cast<const I3CLSimRandomValueApplyFunction &>(other);

        if (other_.applyFunctionName_ != applyFunctionName_) return false;
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
void I3CLSimRandomValueApplyFunction::serialize(Archive &ar, unsigned version)
{
    if (version>i3clsimrandomvalueapplyfunction_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimRandomValueApplyFunction class.",
                  version,
                  i3clsimrandomvalueapplyfunction_version_);

    ar & make_nvp("I3CLSimRandomValue", base_object<I3CLSimRandomValue>(*this));
    ar & make_nvp("applyFunctionName", applyFunctionName_);
    ar & make_nvp("randomDistUsed", randomDistUsed_);
}     


I3_SERIALIZABLE(I3CLSimRandomValueApplyFunction);
