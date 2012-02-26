#include <icetray/serialization.h>
#include <icetray/I3Logging.h>
#include <clsim/I3CLSimRandomValueSimplifiedLiu.h>

#include "clsim/to_float_string.h"
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
    "    //const float g = " + to_float_string(meanCosine_) + ";\n"
    "    //const float beta = (1.f-g)/(1.f+g);\n"
    "    const float beta = " + to_float_string((1.-meanCosine_)/(1.+meanCosine_)) + ";\n"
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
