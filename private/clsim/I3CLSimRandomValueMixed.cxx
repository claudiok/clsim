#include <icetray/serialization.h>
#include <icetray/I3Logging.h>
#include <clsim/I3CLSimRandomValueMixed.h>

#include "clsim/to_float_string.h"
using namespace I3CLSimHelper;

I3CLSimRandomValueMixed::I3CLSimRandomValueMixed(double fractionOfFirstDistribution,
                                                 I3CLSimRandomValueConstPtr firstDistribution,
                                                 I3CLSimRandomValueConstPtr secondDistribution)
:
fractionOfFirstDistribution_(fractionOfFirstDistribution),
firstDistribution_(firstDistribution),
secondDistribution_(secondDistribution)
{ 
    if ((fractionOfFirstDistribution_ < 0.) || (fractionOfFirstDistribution_ > 1.))
        log_fatal("The \"fractionOfFirstDistribution\" parameter must not be <0 or >1. It is %f.",
                  fractionOfFirstDistribution_);
    
    if (!firstDistribution_)
        log_fatal("The \"firstDistribution\" parameter is NULL!");

    if (!secondDistribution_)
        log_fatal("The \"secondDistribution\" parameter is NULL!");
}

I3CLSimRandomValueMixed::~I3CLSimRandomValueMixed() 
{ 

}

I3CLSimRandomValueMixed::I3CLSimRandomValueMixed() {;}

bool I3CLSimRandomValueMixed::OpenCLFunctionWillOnlyUseASingleRandomNumber() const
{
    return ((firstDistribution_->OpenCLFunctionWillOnlyUseASingleRandomNumber()) && 
            (secondDistribution_->OpenCLFunctionWillOnlyUseASingleRandomNumber()));
}

std::string I3CLSimRandomValueMixed::GetOpenCLFunction
(const std::string &functionName,
 const std::string &functionArgs,
 const std::string &functionArgsToCall,
 const std::string &uniformRandomCall_co,
 const std::string &uniformRandomCall_oc
 ) const
{
    std::string firstFunc, secondFunc, mainFunc;
    
    if (OpenCLFunctionWillOnlyUseASingleRandomNumber())
    {
        firstFunc = 
        firstDistribution_->GetOpenCLFunction(functionName+"_mix1",
                                              "float rrrr__",
                                              "rrrr__",
                                              "rrrr__",
                                              "(1.f-rrrr__)");
        secondFunc = 
        secondDistribution_->GetOpenCLFunction(functionName+"_mix2",
                                              "float rrrr__",
                                              "rrrr__",
                                              "rrrr__",
                                              "(1.f-rrrr__)");
        
        const std::string functionDecl = std::string("inline float ") + functionName + "(" + functionArgs + ")";
        
        mainFunc = functionDecl + ";\n\n" + functionDecl + "\n"
        "{\n"    
        "    const float rr = " + uniformRandomCall_co + ";\n"
        "    if (rr < " + to_float_string(fractionOfFirstDistribution_) + ")\n"
        "    {\n"
        "        return " + functionName + "_mix1(rr/" + to_float_string(fractionOfFirstDistribution_) + ");\n"
        "    }\n"
        "    else\n"
        "    {\n"
        "        return " + functionName + "_mix2((1.f-rr)/" + to_float_string(1.-fractionOfFirstDistribution_) + ");\n"
        "    }\n"
        "}\n"
        ;
    }
    else        
    {
        firstFunc = 
        firstDistribution_->GetOpenCLFunction(functionName+"_mix1",
                                              functionArgs,
                                              functionArgsToCall,
                                              uniformRandomCall_co,
                                              uniformRandomCall_oc);
        secondFunc = 
        secondDistribution_->GetOpenCLFunction(functionName+"_mix2",
                                              functionArgs,
                                              functionArgsToCall,
                                              uniformRandomCall_co,
                                              uniformRandomCall_oc);
        
        const std::string functionDecl = std::string("inline float ") + functionName + "(" + functionArgs + ")";
        
        mainFunc = functionDecl + ";\n\n" + functionDecl + "\n"
        "{\n"    
        "    const float rr = " + uniformRandomCall_co + ";\n"
        "    if (rr < " + to_float_string(fractionOfFirstDistribution_) + ")\n"
        "    {\n"
        "        return " + functionName + "_mix1(" + functionArgsToCall + ");\n"
        "    }\n"
        "    else\n"
        "    {\n"
        "        return " + functionName + "_mix2(" + functionArgsToCall + ");\n"
        "    }\n"
        "}\n"
        ;
    }
    
    return std::string() + 
    "\n"
    "/////// BEGIN mix scattering angle generator \"" + functionName + "\" ////////\n" +
    firstFunc + "\n" +
    secondFunc + "\n" +
    mainFunc + "\n"
    "///////   END mix scattering angle generator \"" + functionName + "\" ////////\n"
    ;
}

bool I3CLSimRandomValueMixed::CompareTo(const I3CLSimRandomValue &other) const
{
    try
    {
        const I3CLSimRandomValueMixed &other_ = dynamic_cast<const I3CLSimRandomValueMixed &>(other);
        if (other_.fractionOfFirstDistribution_ != fractionOfFirstDistribution_) return false;
        
        if (!(firstDistribution_->CompareTo(*(other_.firstDistribution_)))) return false;
        if (!(secondDistribution_->CompareTo(*(other_.secondDistribution_)))) return false;
        
        return true;
    }
    catch (const std::bad_cast& e)
    {
        // not of the same type, treat it as non-equal
        return false;
    }
    
}

namespace {
    template <typename T, class Archive>
    void LoadFromArchiveIntoConstPtr(Archive &ar, const std::string &name, boost::shared_ptr<const T> &arg)
    {
        shared_ptr<T> tmp;
        ar >> make_nvp(name.c_str(), tmp);
        arg = tmp;
    }
}


template <class Archive>
void I3CLSimRandomValueMixed::load(Archive &ar, unsigned version)
{
    if (version>i3clsimrandomvaluemixed_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimRandomValueMixed class.",version,i3clsimrandomvaluemixed_version_);

    ar >> make_nvp("I3CLSimRandomValue", base_object<I3CLSimRandomValue>(*this));
    ar >> make_nvp("fractionOfFirstDistribution", fractionOfFirstDistribution_);
    LoadFromArchiveIntoConstPtr(ar, "firstDistribution", firstDistribution_ );
    LoadFromArchiveIntoConstPtr(ar, "secondDistribution", secondDistribution_ );
}     


template <class Archive>
void I3CLSimRandomValueMixed::save(Archive &ar, unsigned version) const
{
    ar << make_nvp("I3CLSimRandomValue", base_object<I3CLSimRandomValue>(*this));
    ar << make_nvp("fractionOfFirstDistribution", fractionOfFirstDistribution_);
    ar << make_nvp("firstDistribution", firstDistribution_);
    ar << make_nvp("secondDistribution", secondDistribution_);
}     


I3_SERIALIZABLE(I3CLSimRandomValueMixed);
