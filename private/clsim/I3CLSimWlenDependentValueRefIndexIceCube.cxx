#include <icetray/serialization.h>
#include <icetray/I3Units.h>
#include <clsim/I3CLSimWlenDependentValueRefIndexIceCube.h>

#include <typeinfo>
#include <cmath>

#include "clsim/to_float_string.h"
using namespace I3CLSimHelper;

const double I3CLSimWlenDependentValueRefIndexIceCube::default_n0= 1.55749;
const double I3CLSimWlenDependentValueRefIndexIceCube::default_n1=-1.57988;
const double I3CLSimWlenDependentValueRefIndexIceCube::default_n2= 3.99993;
const double I3CLSimWlenDependentValueRefIndexIceCube::default_n3=-4.68271;
const double I3CLSimWlenDependentValueRefIndexIceCube::default_n4= 2.09354;

I3CLSimWlenDependentValueRefIndexIceCube::
I3CLSimWlenDependentValueRefIndexIceCube(double n0,
                                         double n1,
                                         double n2,
                                         double n3,
                                         double n4
                                         )
:
n0_(n0),
n1_(n1),
n2_(n2),
n3_(n3),
n4_(n4)
{ 
}

I3CLSimWlenDependentValueRefIndexIceCube::~I3CLSimWlenDependentValueRefIndexIceCube() 
{;}


double I3CLSimWlenDependentValueRefIndexIceCube::GetValue(double wlen) const
{
    const double x = wlen/I3Units::micrometer;
    return n0_ + x*(n1_ + x*(n2_ + x*(n3_ + x*n4_)));
}

double I3CLSimWlenDependentValueRefIndexIceCube::GetDerivative(double wlen) const
{
    const double x = wlen/I3Units::micrometer;
    return n1_ + x*(2.*n2_ + x*(3.*n3_ + x*4.*n4_));
}

std::string I3CLSimWlenDependentValueRefIndexIceCube::GetOpenCLFunction(const std::string &functionName) const
{
    std::string funcDef = 
    std::string("inline float ") + functionName + std::string("(float wlen)\n");

    std::string funcBody = std::string() + 
    "{\n"
    "    const float n0 = " + to_float_string(n0_) + ";\n"
    "    const float n1 = " + to_float_string(n1_) + ";\n"
    "    const float n2 = " + to_float_string(n2_) + ";\n"
    "    const float n3 = " + to_float_string(n3_) + ";\n"
    "    const float n4 = " + to_float_string(n4_) + ";\n"
    "    \n"
    "    const float x = wlen/" + to_float_string(I3Units::micrometer) + ";\n"
    "    return n0 + x*(n1 + x*(n2 + x*(n3 + x*n4)));\n"
    "}\n"
    ;
    
    return funcDef + funcBody;
}

std::string I3CLSimWlenDependentValueRefIndexIceCube::GetOpenCLFunctionDerivative(const std::string &functionName) const
{
    std::string funcDef = 
    std::string("inline float ") + functionName + std::string("(float wlen)\n");
    
    std::string funcBody = std::string() + 
    "{\n"
    "    const float n1 = " + to_float_string(n1_) + ";\n"
    "    const float n2 = " + to_float_string(n2_) + ";\n"
    "    const float n3 = " + to_float_string(n3_) + ";\n"
    "    const float n4 = " + to_float_string(n4_) + ";\n"
    "    \n"
    "    const float x = wlen/" + to_float_string(I3Units::micrometer) + ";\n"
    "    return n1 + x*(2.f*n2 + x*(3.f*n3 + x*4.f*n4));\n"
    "}\n"
    ;
    
    return funcDef + funcBody;
}

bool I3CLSimWlenDependentValueRefIndexIceCube::CompareTo(const I3CLSimWlenDependentValue &other) const
{
    try
    {
        const I3CLSimWlenDependentValueRefIndexIceCube &other_ = dynamic_cast<const I3CLSimWlenDependentValueRefIndexIceCube &>(other);
        return ((other_.n0_ == n0_) &&
                (other_.n1_ == n1_) &&
                (other_.n2_ == n2_) &&
                (other_.n3_ == n3_) &&
                (other_.n4_ == n4_));
    }
    catch (const std::bad_cast& e)
    {
        // not of the same type, treat it as non-equal
        return false;
    }
    
}


template <class Archive>
void I3CLSimWlenDependentValueRefIndexIceCube::serialize(Archive &ar, unsigned version)
{
    if (version>i3clsimwlendependentvaluerefindexicecube_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimWlenDependentValueRefIndexIceCube class.",version,i3clsimwlendependentvaluerefindexicecube_version_);

    ar & make_nvp("I3CLSimWlenDependentValue", base_object<I3CLSimWlenDependentValue>(*this));
    ar & make_nvp("n0", n0_);
    ar & make_nvp("n1", n1_);
    ar & make_nvp("n2", n2_);
    ar & make_nvp("n3", n3_);
    ar & make_nvp("n4", n4_);
}


I3_SERIALIZABLE(I3CLSimWlenDependentValueRefIndexIceCube);
