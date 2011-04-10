#include <icetray/serialization.h>
#include <clsim/I3CLSimWlenDependentValueConstant.h>

#include <typeinfo>
#include <cmath>
#include <math.h>
#include <sstream>
#include <stdexcept>

#include "clsim/to_float_string.h"
using namespace I3CLSimHelper;

I3CLSimWlenDependentValueConstant::
I3CLSimWlenDependentValueConstant(double value)
:
value_(value)
{
}

I3CLSimWlenDependentValueConstant::I3CLSimWlenDependentValueConstant() {;}

I3CLSimWlenDependentValueConstant::~I3CLSimWlenDependentValueConstant() 
{;}

bool I3CLSimWlenDependentValueConstant::HasNativeImplementation() const 
{
    return true;
}

bool I3CLSimWlenDependentValueConstant::HasDerivative() const 
{
    return true;
}

double I3CLSimWlenDependentValueConstant::GetValue(double wlen) const
{
    return value_;
}

double I3CLSimWlenDependentValueConstant::GetDerivative(double wlen) const
{
    return 0.;
}

double I3CLSimWlenDependentValueConstant::GetMinWlen() const
{
    return -std::numeric_limits<double>::infinity();
}

double I3CLSimWlenDependentValueConstant::GetMaxWlen() const
{
    return std::numeric_limits<double>::infinity();
}

std::string I3CLSimWlenDependentValueConstant::GetOpenCLFunction(const std::string &functionName) const
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
    
    return funcDef + funcBody;
}

std::string I3CLSimWlenDependentValueConstant::GetOpenCLFunctionDerivative(const std::string &functionName) const
{
    std::string funcDef = 
    std::string("inline float ") + functionName + std::string("(float wavelength)\n");
    
    std::string funcBody = std::string() + 
    "{\n"
    "    return 0.f;\n"
    "}\n"
    ;
    
    return funcDef + funcBody;
}

bool I3CLSimWlenDependentValueConstant::CompareTo(const I3CLSimWlenDependentValue &other) const
{
    try
    {
        const I3CLSimWlenDependentValueConstant &other_ = dynamic_cast<const I3CLSimWlenDependentValueConstant &>(other);
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
void I3CLSimWlenDependentValueConstant::serialize(Archive &ar, unsigned version)
{
    if (version>i3clsimwlendependentvalueconstant_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimWlenDependentValueConstant class.",version,i3clsimwlendependentvalueconstant_version_);

    ar & make_nvp("I3CLSimWlenDependentValue", base_object<I3CLSimWlenDependentValue>(*this));
    ar & make_nvp("value", value_);
}     


I3_SERIALIZABLE(I3CLSimWlenDependentValueConstant);
