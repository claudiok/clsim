#ifndef I3CLSIMFUNCTIONSCATLENPARTIC_H_INCLUDED
#define I3CLSIMFUNCTIONSCATLENPARTIC_H_INCLUDED

#include "clsim/I3CLSimFunction.h"

#include <limits>

/**
 * @brief The scattering length as defined by the "partic"
 * scattering model.
 */
static const unsigned i3clsimfunctionscatlenpartic_version_ = 0;

struct I3CLSimFunctionScatLenPartic : public I3CLSimFunction
{
public:
    static const double default_volumeConcentrationSmallParticles;
    static const double default_volumeConcentrationLargeParticles;
    
    I3CLSimFunctionScatLenPartic(double volumeConcentrationSmallParticles=default_volumeConcentrationSmallParticles,   // fraction (e.g. 0.0075*I3Units::perMillion)
                                           double volumeConcentrationLargeParticles=default_volumeConcentrationLargeParticles    // fraction (e.g. 0.0075*I3Units::perMillion)
                                           );
    virtual ~I3CLSimFunctionScatLenPartic();
    
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
     * Shall compare to another I3CLSimFunction object
     */
    virtual bool CompareTo(const I3CLSimFunction &other) const;
    
private:
    double volumeConcentrationSmallParticles_;
    double volumeConcentrationLargeParticles_;
    
    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive & ar, unsigned version);
};


BOOST_CLASS_VERSION(I3CLSimFunctionScatLenPartic, i3clsimfunctionscatlenpartic_version_);

I3_POINTER_TYPEDEFS(I3CLSimFunctionScatLenPartic);

#endif //I3CLSIMFUNCTIONSCATLENPARTIC_H_INCLUDED
