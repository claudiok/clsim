#include <icetray/serialization.h>
#include <icetray/I3Units.h>
#include <clsim/I3CLSimWlenDependentValuePolynomial.h>

#include <typeinfo>
#include <cmath>

#include "clsim/to_float_string.h"
using namespace I3CLSimHelper;

I3CLSimWlenDependentValuePolynomial::
I3CLSimWlenDependentValuePolynomial(const std::vector<double> &coeffs)
:
coefficients_(coeffs)
{ 
}

I3CLSimWlenDependentValuePolynomial::I3CLSimWlenDependentValuePolynomial() {;}

I3CLSimWlenDependentValuePolynomial::~I3CLSimWlenDependentValuePolynomial() 
{;}


double I3CLSimWlenDependentValuePolynomial::GetValue(double wlen) const
{
    if (coefficients_.size()==0) return 0.;
    
    double sum=coefficients_[0];
    double multiplier=1.;
    
    for (std::size_t i=1;i<coefficients_.size();++i)
    {
        multiplier *= wlen;
        sum += coefficients_[i]*multiplier;
    }
    
    return sum;
}


std::string I3CLSimWlenDependentValuePolynomial::GetOpenCLFunction(const std::string &functionName) const
{
    std::ostringstream output(std::ostringstream::out);
    output.setf(std::ios::scientific,std::ios::floatfield);
    output.precision(std::numeric_limits<float>::digits10+4); // maximum precision for a float

    output << "inline float " << functionName << "(float x)" << std::endl;
    output << "{" << std::endl;

    if (coefficients_.size()==0)
    {
        output << "return 0.f;" << std::endl;
    }
    else if (coefficients_.size()==1)
    {
        output << "return " << to_float_string(coefficients_[0]) << "f;" << std::endl;
    }
    else // 2 or more coefficients
    {
        output << "return ";
        
        for (std::size_t i=0;i<coefficients_.size()-1;++i)
        {
            output << to_float_string(coefficients_[i]) << " + x*(";
        }

        output << to_float_string(coefficients_[coefficients_.size()-1]);
        
        for (std::size_t i=0;i<coefficients_.size()-1;++i)
        {
            output << ")";
        }
        output << ";" << std::endl;
    }
    
    output << "}" << std::endl;

    return output.str();
}

bool I3CLSimWlenDependentValuePolynomial::CompareTo(const I3CLSimWlenDependentValue &other) const
{
    try
    {
        const I3CLSimWlenDependentValuePolynomial &other_ = dynamic_cast<const I3CLSimWlenDependentValuePolynomial &>(other);
        
        if (other_.coefficients_.size() != coefficients_.size()) return false;
        
        for (std::size_t i=0;i<coefficients_.size();++i)
        {
            if (other_.coefficients_[i] != coefficients_[i]) return false;
        }
        
        return true;
    }
    catch (const std::bad_cast& e)
    {
        // not of the same type, treat it as non-equal
        return false;
    }
    
}


template <class Archive>
void I3CLSimWlenDependentValuePolynomial::serialize(Archive &ar, unsigned version)
{
    if (version>i3clsimwlendependentvaluepolynomial_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimWlenDependentValuePolynomial class.",
                  version,i3clsimwlendependentvaluepolynomial_version_);

    ar & make_nvp("I3CLSimWlenDependentValue", base_object<I3CLSimWlenDependentValue>(*this));
    ar & make_nvp("coefficients", coefficients_);
}


I3_SERIALIZABLE(I3CLSimWlenDependentValuePolynomial);
