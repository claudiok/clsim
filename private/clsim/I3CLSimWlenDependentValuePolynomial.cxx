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
coefficients_(coeffs),
rangemin_(-std::numeric_limits<double>::infinity()),
rangemax_(std::numeric_limits<double>::infinity()),
underflow_(NAN),
overflow_(NAN)
{ 
}

I3CLSimWlenDependentValuePolynomial::
I3CLSimWlenDependentValuePolynomial(const std::vector<double> &coeffs,
                                    double rangemin, double rangemax)
:
coefficients_(coeffs),
rangemin_(rangemin),
rangemax_(rangemax),
underflow_(GetValue(rangemin)),
overflow_(GetValue(rangemax))
{ 
    if (rangemax_ <= rangemin_) log_fatal("Trying to initalize a polynomial with a invalid range! [%f,%f]", rangemin_, rangemax_);
}



I3CLSimWlenDependentValuePolynomial::
I3CLSimWlenDependentValuePolynomial(const std::vector<double> &coeffs,
                                    double rangemin, double rangemax,
                                    double underflow, double overflow)
:
coefficients_(coeffs),
rangemin_(rangemin),
rangemax_(rangemax),
underflow_(underflow),
overflow_(overflow)
{ 
    if (rangemax_ <= rangemin_) log_fatal("Trying to initalize a polynomial with a invalid range! [%f,%f]", rangemin_, rangemax_);
}


I3CLSimWlenDependentValuePolynomial::I3CLSimWlenDependentValuePolynomial()
{;}

I3CLSimWlenDependentValuePolynomial::~I3CLSimWlenDependentValuePolynomial() 
{;}


double I3CLSimWlenDependentValuePolynomial::GetValue(double wlen) const
{
    if (coefficients_.size()==0) return 0.;
    if (wlen < rangemin_) return underflow_;
    else if (wlen > rangemax_) return overflow_;
    
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

    output << "inline float " << functionName << "(float x);" << std::endl;
    output << "inline float " << functionName << "(float x)" << std::endl;
    output << "{" << std::endl;
    
    // check the range
    if (!std::isinf(rangemin_)) {
        output << "if (x < " << to_float_string(rangemin_) << ")" << std::endl;
        output << "{" << std::endl;
        output << "    return " << to_float_string(underflow_) << ";" << std::endl;
        output << "}" << std::endl;
    }
    
    if (!std::isinf(rangemax_)) {
        output << "if (x > " << to_float_string(rangemax_) << ")" << std::endl;
        output << "{" << std::endl;
        output << "    return " << to_float_string(overflow_) << ";" << std::endl;
        output << "}" << std::endl;
    }

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
        if (other_.rangemin_ != rangemin_) return false;
        if (other_.rangemax_ != rangemax_) return false;
        if (other_.underflow_ != underflow_) return false;
        if (other_.overflow_ != overflow_) return false;
        
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
    ar & make_nvp("rangemin", rangemin_);
    ar & make_nvp("rangemax", rangemax_);
    ar & make_nvp("underflow", underflow_);
    ar & make_nvp("overflow", overflow_);
}


I3_SERIALIZABLE(I3CLSimWlenDependentValuePolynomial);
