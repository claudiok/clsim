#include <icetray/serialization.h>
#include <icetray/I3Logging.h>
#include <clsim/I3CLSimRandomValueHenyeyGreenstein.h>

#include "clsim/to_float_string.h"
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

std::string I3CLSimRandomValueHenyeyGreenstein::GetOpenCLFunction
(const std::string &functionName,
 const std::string &functionArgs,
 const std::string &functionArgsToCall,
 const std::string &uniformRandomCall_co,
 const std::string &uniformRandomCall_oc
 ) const
{
    const std::string functionDecl = std::string("inline float ") + functionName + "(" + functionArgs + ")";
    
    return functionDecl + "\n"
    "{\n"
    "    const float g = " + to_float_string(meanCosine_) + ";\n"
    "    const float g2 = " + to_float_string(meanCosine_*meanCosine_) + ";\n"
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
