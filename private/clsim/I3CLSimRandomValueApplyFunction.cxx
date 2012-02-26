#include <icetray/serialization.h>
#include <icetray/I3Logging.h>
#include <clsim/I3CLSimRandomValueApplyFunction.h>

#include <boost/lexical_cast.hpp>

#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <limits>

I3CLSimRandomValueApplyFunction::
I3CLSimRandomValueApplyFunction(I3CLSimRandomValuePtr randomDistUsed,
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

    const std::string functionDecl = std::string("inline float ") + functionName + "(" + functionArgs + ")";
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
    "    const float number = " + functionNameUsed + "(" + functionArgsToCall + ");\n"
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
