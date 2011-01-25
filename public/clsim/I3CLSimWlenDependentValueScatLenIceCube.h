#ifndef I3CLSIMWLENDEPENDENTVALUESCATLENICECUBE_H_INCLUDED
#define I3CLSIMWLENDEPENDENTVALUESCATLENICECUBE_H_INCLUDED

#include "clsim/I3CLSimWlenDependentValue.h"

#include <limits>

/**
 * @brief Six-parameter ice model as fitted in the SPICE ice model.
 * (scattering length)
 */
static const unsigned i3clsimwlendependentvaluescatlenicecube_version_ = 0;

struct I3CLSimWlenDependentValueScatLenIceCube : public I3CLSimWlenDependentValue
{
public:
    I3CLSimWlenDependentValueScatLenIceCube(double alpha,
                                            double be400
                                            );
    virtual ~I3CLSimWlenDependentValueScatLenIceCube();
    
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
    inline double GetAlpha() const {return alpha_;}
    inline double GetBe400() const {return be400_;}
    
private:
    I3CLSimWlenDependentValueScatLenIceCube();
    
    double alpha_;
    double be400_;
    
    
    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive & ar, unsigned version);
};


BOOST_CLASS_VERSION(I3CLSimWlenDependentValueScatLenIceCube, i3clsimwlendependentvaluescatlenicecube_version_);

I3_POINTER_TYPEDEFS(I3CLSimWlenDependentValueScatLenIceCube);

#endif //I3CLSIMWLENDEPENDENTVALUESCATLENICECUBE_H_INCLUDED
