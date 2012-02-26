#include <icetray/serialization.h>
#include <icetray/I3Units.h>
#include <clsim/I3CLSimWlenDependentValueAbsLenIceCube.h>

#include <typeinfo>
#include <cmath>

#include "clsim/to_float_string.h"
using namespace I3CLSimHelper;

I3CLSimWlenDependentValueAbsLenIceCube::
I3CLSimWlenDependentValueAbsLenIceCube(double kappa,
                                       double A,
                                       double B,
                                       double D,
                                       double E,
                                       double aDust400,
                                       double deltaTau
                                       )
:
kappa_(kappa),
A_(A),
B_(B),
D_(D),
E_(E),
aDust400_(aDust400),
deltaTau_(deltaTau)
{ 
}

I3CLSimWlenDependentValueAbsLenIceCube::I3CLSimWlenDependentValueAbsLenIceCube() {;}

I3CLSimWlenDependentValueAbsLenIceCube::~I3CLSimWlenDependentValueAbsLenIceCube() 
{;}


double I3CLSimWlenDependentValueAbsLenIceCube::GetValue(double wlen) const
{
    const double x = wlen/I3Units::nanometer;
    return I3Units::m/( (D_*aDust400_+E_) * std::pow(x, -kappa_)  +  A_*std::exp(-B_/x) * (1.f + 0.01f*deltaTau_) );
}


std::string I3CLSimWlenDependentValueAbsLenIceCube::GetOpenCLFunction(const std::string &functionName) const
{
    std::string funcDef = 
    std::string("inline float ") + functionName + std::string("(float wlen)\n");

    std::string funcBody = std::string() + 
    "{\n"
    "    const float kappa = " + to_float_string(kappa_) + ";\n"
    "    const float A = " + to_float_string(A_) + ";\n"
    "    const float B = " + to_float_string(B_) + ";\n"
    "    const float D = " + to_float_string(D_) + ";\n"
    "    const float E = " + to_float_string(E_) + ";\n"
    "    const float aDust400 = " + to_float_string(aDust400_) + ";\n"
    "    const float deltaTau = " + to_float_string(deltaTau_) + ";\n"
    "    \n"
    "    const float x = wlen/" + to_float_string(I3Units::nanometer) + ";\n"
    "    \n"
    "#ifdef USE_NATIVE_MATH\n"
    "    return " + to_float_string(I3Units::m) + "*native_recip( (D*aDust400+E) * native_powr(x, -kappa)  +  A*native_exp(-B/x) * (1.f + 0.01f*deltaTau) );\n"
    "#else\n"
    "    return " + to_float_string(I3Units::m) + "/( (D*aDust400+E) * powr(x, -kappa)  +  A*exp(-B/x) * (1.f + 0.01f*deltaTau) );\n"
    "#endif\n"
    "}\n"
    ;
    
    return funcDef + ";\n\n" + funcDef + funcBody;
}

bool I3CLSimWlenDependentValueAbsLenIceCube::CompareTo(const I3CLSimWlenDependentValue &other) const
{
    try
    {
        const I3CLSimWlenDependentValueAbsLenIceCube &other_ = dynamic_cast<const I3CLSimWlenDependentValueAbsLenIceCube &>(other);
        return ((other_.kappa_ == kappa_) &&
                (other_.A_ == A_) &&
                (other_.B_ == B_) &&
                (other_.D_ == D_) &&
                (other_.E_ == E_) &&
                (other_.aDust400_ == aDust400_) &&
                (other_.deltaTau_ == deltaTau_));
    }
    catch (const std::bad_cast& e)
    {
        // not of the same type, treat it as non-equal
        return false;
    }
    
}


template <class Archive>
void I3CLSimWlenDependentValueAbsLenIceCube::serialize(Archive &ar, unsigned version)
{
    if (version>i3clsimwlendependentvalueabslenicecube_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimWlenDependentValueAbsLenIceCube class.",version,i3clsimwlendependentvalueabslenicecube_version_);

    ar & make_nvp("I3CLSimWlenDependentValue", base_object<I3CLSimWlenDependentValue>(*this));
    ar & make_nvp("kappa", kappa_);
    ar & make_nvp("A", A_);
    ar & make_nvp("B", B_);
    ar & make_nvp("D", D_);
    ar & make_nvp("E", E_);
    ar & make_nvp("aDust400", aDust400_);
    ar & make_nvp("deltaTau", deltaTau_);
}


I3_SERIALIZABLE(I3CLSimWlenDependentValueAbsLenIceCube);
