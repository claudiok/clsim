#ifndef I3CLSimWlenDependentValueConstant_H_INCLUDED
#define I3CLSimWlenDependentValueConstant_H_INCLUDED

#include "clsim/I3CLSimWlenDependentValue.h"

#include <vector>

/**
 * @brief A constant value.
 */
static const unsigned i3clsimwlendependentvalueconstant_version_ = 0;

struct I3CLSimWlenDependentValueConstant : public I3CLSimWlenDependentValue
{
public:
    
    I3CLSimWlenDependentValueConstant(double value);
    virtual ~I3CLSimWlenDependentValueConstant();
    
    /**
     * If this is true, it is assumed that GetValue() and GetDerivative() return
     * meaningful values. If not, GetValue will not be called;
     * only the OpenCL implementation will be used.
     */
    virtual bool HasNativeImplementation() const; //{return true;};
    
    /**
     * If this is true, derivatives can be used.
     */
    virtual bool HasDerivative() const; //{return true;};
    
    /**
     * Shall return the value at a requested wavelength (n)
     */
    virtual double GetValue(double wlen) const;
    
    /**
     * Shall return the derivative at a requested wavelength (dn/dlambda)
     */
    virtual double GetDerivative(double wlen) const;

    /**
     * Shall return the minimal supported wavelength (possibly -inf)
     */
    virtual double GetMinWlen() const;
    
    /**
     * Shall return the maximal supported wavelength (possibly +inf)
     */
    virtual double GetMaxWlen() const;
    
    /**
     * Shall return an OpenCL-compatible function named
     * functionName with a single float argument (float wlen)
     */
    virtual std::string GetOpenCLFunction(const std::string &functionName) const;
    
    /**
     * Shall return an OpenCL-compatible function named
     * functionName with a single float argument (float wlen)
     */
    virtual std::string GetOpenCLFunctionDerivative(const std::string &functionName) const;

    /**
     * Shall compare to another I3CLSimWlenDependentValue object
     */
    virtual bool CompareTo(const I3CLSimWlenDependentValue &other) const;
    
private:
    I3CLSimWlenDependentValueConstant();
    
    double value_;
    
    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive & ar, unsigned version);
};


BOOST_CLASS_VERSION(I3CLSimWlenDependentValueConstant, i3clsimwlendependentvalueconstant_version_);

I3_POINTER_TYPEDEFS(I3CLSimWlenDependentValueConstant);

#endif //I3CLSimWlenDependentValueConstant_H_INCLUDED
