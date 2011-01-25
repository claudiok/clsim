#ifndef I3CLSIMWLENDEPENDENTVALUEREFINDEXICECUBE_H_INCLUDED
#define I3CLSIMWLENDEPENDENTVALUEREFINDEXICECUBE_H_INCLUDED

#include "clsim/I3CLSimWlenDependentValue.h"

#include <limits>

/**
 * @brief The phase refractive index for IceCube glacial ice, taken from
 * "Role of Group and Phase Velocity in High-Energy Neutrino Observatories",
 * P.B. Price and K. Woschnagg, Astropart. Phys. 15 (2001) 97
 */
static const unsigned i3clsimwlendependentvaluerefindexicecube_version_ = 0;

struct I3CLSimWlenDependentValueRefIndexIceCube : public I3CLSimWlenDependentValue
{
public:
    static const double default_n0;
    static const double default_n1;
    static const double default_n2;
    static const double default_n3;
    static const double default_n4;
    
    
    I3CLSimWlenDependentValueRefIndexIceCube(double n0 = default_n0,  // coefficients 0-4
                                             double n1 = default_n1,
                                             double n2 = default_n2,
                                             double n3 = default_n3,
                                             double n4 = default_n4
                                             );
    virtual ~I3CLSimWlenDependentValueRefIndexIceCube();
    
    /**
     * If this is true, it is assumed that GetValue() and GetDerivative() return
     * meaningful values. If not, GetValue will not be called;
     * only the OpenCL implementation will be used.
     */
    virtual bool HasNativeImplementation() const {return true;};
    
    /**
     * If this is true, derivatives can be used.
     */
    virtual bool HasDerivative() const {return true;};
    
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
     * Shall return an OpenCL-compatible function named
     * functionName with a single float argument (float wlen)
     */
    virtual std::string GetOpenCLFunctionDerivative(const std::string &functionName) const;
    
    /**
     * Shall compare to another I3CLSimWlenDependentValue object
     */
    virtual bool CompareTo(const I3CLSimWlenDependentValue &other) const;
    
private:
    double n0_;
    double n1_;
    double n2_;
    double n3_;
    double n4_;
    
    
    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive & ar, unsigned version);
};


BOOST_CLASS_VERSION(I3CLSimWlenDependentValueRefIndexIceCube, i3clsimwlendependentvaluerefindexicecube_version_);

I3_POINTER_TYPEDEFS(I3CLSimWlenDependentValueRefIndexIceCube);

#endif //I3CLSIMWLENDEPENDENTVALUEREFINDEXICECUBE_H_INCLUDED
