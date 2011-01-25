#ifndef I3CLSIMWLENDEPENDENTVALUE_H_INCLUDED
#define I3CLSIMWLENDEPENDENTVALUE_H_INCLUDED

#include "icetray/serialization.h"
#include "icetray/I3TrayHeaders.h"

#include <string>

/**
 * @brief A value dependent on photon wavelength
 */
static const unsigned i3clsimwlendependentvalue_version_ = 0;

struct I3CLSimWlenDependentValue 
{
public:
    
    I3CLSimWlenDependentValue();
    virtual ~I3CLSimWlenDependentValue();

    /**
     * If this is true, it is assumed that GetValue() and GetDerivative() return
     * meaningful values. If not, GetValue will not be called;
     * only the OpenCL implementation will be used.
     */
    virtual bool HasNativeImplementation() const = 0;

    /**
     * If this is true, derivatives can be used.
     */
    virtual bool HasDerivative() const = 0;
    
    /**
     * Shall return the value at a requested wavelength (n)
     */
    virtual double GetValue(double wlen) const = 0;

    /**
     * Shall return the derivative at a requested wavelength (dn/dlambda)
     */
    virtual double GetDerivative(double wlen) const {return NAN;}
    
    /**
     * Shall return the minimal supported wavelength (possibly -inf)
     */
    virtual double GetMinWlen() const = 0;

    /**
     * Shall return the maximal supported wavelength (possibly +inf)
     */
    virtual double GetMaxWlen() const = 0;
    
    /**
     * Shall return an OpenCL-compatible function named
     * functionName with a single float argument (float wlen)
     */
    virtual std::string GetOpenCLFunction(const std::string &functionName) const = 0;

    /**
     * Shall return an OpenCL-compatible function named
     * functionName with a single float argument (float wlen)
     */
    virtual std::string GetOpenCLFunctionDerivative(const std::string &functionName) const {return std::string();}
    
    /**
     * Shall compare to another I3CLSimWlenDependentValue object
     */
    virtual bool CompareTo(const I3CLSimWlenDependentValue &other) const = 0;
    
private:
    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive & ar, unsigned version);
};

inline bool operator==(const I3CLSimWlenDependentValue& a, const I3CLSimWlenDependentValue& b)
{
    return a.CompareTo(b);
}

inline bool operator!=(const I3CLSimWlenDependentValue& a, const I3CLSimWlenDependentValue& b)
{
    return (!a.CompareTo(b));
}


BOOST_CLASS_VERSION(I3CLSimWlenDependentValue, i3clsimwlendependentvalue_version_);

I3_POINTER_TYPEDEFS(I3CLSimWlenDependentValue);

#endif //I3CLSIMWLENDEPENDENTVALUE_H_INCLUDED
