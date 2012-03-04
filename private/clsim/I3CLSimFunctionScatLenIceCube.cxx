#include <icetray/serialization.h>
#include <icetray/I3Units.h>
#include <clsim/I3CLSimFunctionScatLenIceCube.h>

#include <typeinfo>
#include <cmath>

#include "clsim/to_float_string.h"
using namespace I3CLSimHelper;

I3CLSimFunctionScatLenIceCube::
I3CLSimFunctionScatLenIceCube(double alpha,
                                        double b400
                                        )
:
alpha_(alpha),
b400_(b400)
{ 
}

I3CLSimFunctionScatLenIceCube::I3CLSimFunctionScatLenIceCube() {;}

I3CLSimFunctionScatLenIceCube::~I3CLSimFunctionScatLenIceCube() 
{;}


double I3CLSimFunctionScatLenIceCube::GetValue(double wlen) const
{
    const double x = wlen/I3Units::nanometer;
    return I3Units::m/( b400_ * std::pow(x/400., -alpha_) );
}


std::string I3CLSimFunctionScatLenIceCube::GetOpenCLFunction(const std::string &functionName) const
{
    std::string funcDef = 
    std::string("inline float ") + functionName + std::string("(float wlen)\n");

    const std::string meterAsString = to_float_string(I3Units::m);
    const std::string refWlenAsString = to_float_string(1./(400.*I3Units::nanometer));
    
    
    std::string funcBody = std::string() + 
    "{\n"
    "    const float alpha = " + to_float_string(alpha_) + ";\n"
    "    const float b400 = " + to_float_string(b400_) + ";\n"
    "    \n"
    "#ifdef USE_NATIVE_MATH\n"
    "    return " + to_float_string(I3Units::m) + "*native_recip( b400 * native_powr(wlen*" + refWlenAsString + ", -alpha) );\n"
    "#else\n"
    "    return " + to_float_string(I3Units::m) + "/( b400 * powr(wlen*" + refWlenAsString + ", -alpha) );\n"
    "#endif\n"
    "}\n"
    ;
    
    return funcDef + ";\n\n" + funcDef + funcBody;
}

bool I3CLSimFunctionScatLenIceCube::CompareTo(const I3CLSimFunction &other) const
{
    try
    {
        const I3CLSimFunctionScatLenIceCube &other_ = dynamic_cast<const I3CLSimFunctionScatLenIceCube &>(other);
        return ((other_.alpha_ == alpha_) &&
                (other_.b400_ == b400_));
    }
    catch (const std::bad_cast& e)
    {
        // not of the same type, treat it as non-equal
        return false;
    }
    
}


template <class Archive>
void I3CLSimFunctionScatLenIceCube::serialize(Archive &ar, unsigned version)
{
    if (version>i3clsimfunctionscatlenicecube_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimFunctionScatLenIceCube class.",version,i3clsimfunctionscatlenicecube_version_);

    ar & make_nvp("I3CLSimFunction", base_object<I3CLSimFunction>(*this));
    ar & make_nvp("alpha", alpha_);
    ar & make_nvp("b400", b400_);
}


I3_SERIALIZABLE(I3CLSimFunctionScatLenIceCube);
