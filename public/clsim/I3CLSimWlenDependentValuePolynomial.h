#ifndef I3CLSIMWLENDEPENDENTVALUEPOLYNOMIAL_H_INCLUDED
#define I3CLSIMWLENDEPENDENTVALUEPOLYNOMIAL_H_INCLUDED

#include "clsim/I3CLSimWlenDependentValue.h"

#include <vector>

/**
 * @brief A simple polynomial
 */
static const unsigned i3clsimwlendependentvaluepolynomial_version_ = 0;

struct I3CLSimWlenDependentValuePolynomial : public I3CLSimWlenDependentValue
{
public:
    I3CLSimWlenDependentValuePolynomial(const std::vector<double> &coeffs);
    virtual ~I3CLSimWlenDependentValuePolynomial();
    
    /**
     * If this is true, it is assumed that GetValue() and GetDerivative() return
     * meaningful values. If not, GetValue will not be called;
     * only the OpenCL implementation will be used.
     */
    virtual bool HasNativeImplementation() const {return true;};
    
    /**
     * If this is true, derivatives can be used.
     */
    virtual bool HasDerivative() const {return false;};
    
    /**
     * Shall return the value at a requested wavelength (n)
     */
    virtual double GetValue(double wlen) const;
    
    /**
     * Shall return the minimal supported wavelength (possibly -inf)
     */
    virtual double GetMinWlen() const {return -std::numeric_limits<double>::infinity();}
    
    /**
     * Shall return the maximal supported wavelength (possibly +inf)
     */
    virtual double GetMaxWlen() const {return std::numeric_limits<double>::infinity();}
    
    /**
     * Shall return an OpenCL-compatible function named
     * functionName with a single float argument (float wlen)
     */
    virtual std::string GetOpenCLFunction(const std::string &functionName) const;
    
    /**
     * Shall compare to another I3CLSimWlenDependentValue object
     */
    virtual bool CompareTo(const I3CLSimWlenDependentValue &other) const;
    
    
    // access to the internal state
    inline const std::vector<double> &GetCoefficients() const {return coefficients_;};
    
private:
    I3CLSimWlenDependentValuePolynomial();
    
    std::vector<double> coefficients_;
    
    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive & ar, unsigned version);
};


BOOST_CLASS_VERSION(I3CLSimWlenDependentValuePolynomial, i3clsimwlendependentvaluepolynomial_version_);

I3_POINTER_TYPEDEFS(I3CLSimWlenDependentValuePolynomial);

#endif //I3CLSIMWLENDEPENDENTVALUEPOLYNOMIAL_H_INCLUDED
