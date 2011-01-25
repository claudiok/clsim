#include <icetray/serialization.h>
#include <icetray/I3Units.h>
#include <clsim/I3CLSimWlenDependentValueScatLenIceCube.h>

#include <typeinfo>
#include <cmath>

#include "clsim/to_float_string.h"
using namespace I3CLSimHelper;

I3CLSimWlenDependentValueScatLenIceCube::
I3CLSimWlenDependentValueScatLenIceCube(double alpha,
                                        double be400
                                        )
:
alpha_(alpha),
be400_(be400)
{ 
}

I3CLSimWlenDependentValueScatLenIceCube::I3CLSimWlenDependentValueScatLenIceCube() {;}

I3CLSimWlenDependentValueScatLenIceCube::~I3CLSimWlenDependentValueScatLenIceCube() 
{;}


double I3CLSimWlenDependentValueScatLenIceCube::GetValue(double wlen) const
{
    const double x = wlen/I3Units::nanometer;
    return I3Units::m/( be400_ * std::pow(x/400., -alpha_) );
}


std::string I3CLSimWlenDependentValueScatLenIceCube::GetOpenCLFunction(const std::string &functionName) const
{
    std::string funcDef = 
    std::string("inline float ") + functionName + std::string("(float wlen)\n");

    const std::string meterAsString = to_float_string(I3Units::m);
    const std::string refWlenAsString = to_float_string(1./(400.*I3Units::nanometer));
    
    
    std::string funcBody = std::string() + 
    "{\n"
    "    const float alpha = " + to_float_string(alpha_) + ";\n"
    "    const float be400 = " + to_float_string(be400_) + ";\n"
    "    \n"
    "#ifdef USE_NATIVE_MATH\n"
    "    return " + to_float_string(I3Units::m) + "*native_recip( be400 * native_powr(wlen*" + refWlenAsString + ", -alpha) );\n"
    "#else\n"
    "    return " + to_float_string(I3Units::m) + "/( be400 * powr(wlen*" + refWlenAsString + ", -alpha) );\n"
    "#endif\n"
    "}\n"
    ;
    
    return funcDef + funcBody;
}

bool I3CLSimWlenDependentValueScatLenIceCube::CompareTo(const I3CLSimWlenDependentValue &other) const
{
    try
    {
        const I3CLSimWlenDependentValueScatLenIceCube &other_ = dynamic_cast<const I3CLSimWlenDependentValueScatLenIceCube &>(other);
        return ((other_.alpha_ == alpha_) &&
                (other_.be400_ == be400_));
    }
    catch (const std::bad_cast& e)
    {
        // not of the same type, treat it as non-equal
        return false;
    }
    
}


template <class Archive>
void I3CLSimWlenDependentValueScatLenIceCube::serialize(Archive &ar, unsigned version)
{
    if (version>i3clsimwlendependentvaluescatlenicecube_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimWlenDependentValueScatLenIceCube class.",version,i3clsimwlendependentvaluescatlenicecube_version_);

    ar & make_nvp("I3CLSimWlenDependentValue", base_object<I3CLSimWlenDependentValue>(*this));
    ar & make_nvp("alpha", alpha_);
    ar & make_nvp("be400", be400_);
}


I3_SERIALIZABLE(I3CLSimWlenDependentValueScatLenIceCube);
