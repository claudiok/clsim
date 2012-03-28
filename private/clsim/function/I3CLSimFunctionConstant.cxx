#include <icetray/serialization.h>
#include <clsim/function/I3CLSimFunctionConstant.h>

#include <typeinfo>
#include <cmath>
#include <math.h>
#include <sstream>
#include <stdexcept>

#include "clsim/I3CLSimHelperToFloatString.h"
using namespace I3CLSimHelper;

I3CLSimFunctionConstant::
I3CLSimFunctionConstant(double value)
:
value_(value)
{
}

I3CLSimFunctionConstant::I3CLSimFunctionConstant() {;}

I3CLSimFunctionConstant::~I3CLSimFunctionConstant() 
{;}

bool I3CLSimFunctionConstant::HasNativeImplementation() const 
{
    return true;
}

bool I3CLSimFunctionConstant::HasDerivative() const 
{
    return true;
}

double I3CLSimFunctionConstant::GetValue(double wlen) const
{
    return value_;
}

double I3CLSimFunctionConstant::GetDerivative(double wlen) const
{
    return 0.;
}

double I3CLSimFunctionConstant::GetMinWlen() const
{
    return -std::numeric_limits<double>::infinity();
}

double I3CLSimFunctionConstant::GetMaxWlen() const
{
    return std::numeric_limits<double>::infinity();
}

std::string I3CLSimFunctionConstant::GetOpenCLFunction(const std::string &functionName) const
{
    std::string funcDef = 
    std::string("inline float ") + functionName + std::string("(float wavelength)\n");
    
    std::ostringstream output(std::ostringstream::out);
    output.setf(std::ios::scientific,std::ios::floatfield);
    output.precision(std::numeric_limits<float>::digits10+4); // maximum precision for a float
    
    output << value_ << "f";
    std::string valueString = output.str();
    
    std::string funcBody = std::string() + 
    "{\n"
    "    return " + valueString + ";\n"
    "}\n"
    ;
    
    return funcDef + ";\n\n" + funcDef + funcBody;
}

std::string I3CLSimFunctionConstant::GetOpenCLFunctionDerivative(const std::string &functionName) const
{
    std::string funcDef = 
    std::string("inline float ") + functionName + std::string("(float wavelength)\n");
    
    std::string funcBody = std::string() + 
    "{\n"
    "    return 0.f;\n"
    "}\n"
    ;
    
    return funcDef + ";\n\n" + funcDef + funcBody;
}

bool I3CLSimFunctionConstant::CompareTo(const I3CLSimFunction &other) const
{
    try
    {
        const I3CLSimFunctionConstant &other_ = dynamic_cast<const I3CLSimFunctionConstant &>(other);
        if ((other_.value_ != value_))
            return false;
        
        return true;
    }
    catch (const std::bad_cast& e)
    {
        // not of the same type, treat it as non-equal
        return false;
    }
    
}



template <class Archive>
void I3CLSimFunctionConstant::serialize(Archive &ar, unsigned version)
{
    if (version>i3clsimfunctionconstant_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimFunctionConstant class.",version,i3clsimfunctionconstant_version_);

    ar & make_nvp("I3CLSimFunction", base_object<I3CLSimFunction>(*this));
    ar & make_nvp("value", value_);
}     


I3_SERIALIZABLE(I3CLSimFunctionConstant);
