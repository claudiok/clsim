#ifndef I3CLSIMWLENDEPENDENTVALUEFROMTABLE_H_INCLUDED
#define I3CLSIMWLENDEPENDENTVALUEFROMTABLE_H_INCLUDED

#include "clsim/I3CLSimWlenDependentValue.h"

#include <vector>

/**
 * @brief An arbitrary value interpolated from a uniformly distributed number
 * of values.
 */
static const unsigned i3clsimwlendependentvaluefromtable_version_ = 0;

struct I3CLSimWlenDependentValueFromTable : public I3CLSimWlenDependentValue
{
public:
    
    I3CLSimWlenDependentValueFromTable(double startWlen,
                                       double wlenStep,
                                       const std::vector<double> &values);
    virtual ~I3CLSimWlenDependentValueFromTable();
    
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
    virtual double GetMinWlen() const;
    
    /**
     * Shall return the maximal supported wavelength (possibly +inf)
     */
    virtual double GetMaxWlen() const;

    /**
     * Returns the internal state: first wavelength
     */
    inline double GetFirstWavelength() const {return startWlen_;}
    
    /**
     * Returns the internal state: wavelength stepping
     */
    inline double GetWavelengthStepping() const {return wlenStep_;}

    /**
     * Returns the internal state: wavelength stepping
     */
    inline std::size_t GetWavelengthNumValues() const {return values_.size();}

    /**
     * Returns the internal state: wavelength stepping
     */
    inline double GetWavelengthValue(std::size_t i) const {return values_[i];}

    /**
     * Shall return an OpenCL-compatible function named
     * functionName with a single float argument (float wlen)
     */
    virtual std::string GetOpenCLFunction(const std::string &functionName) const;
    
    /**
     * Shall compare to another I3CLSimWlenDependentValue object
     */
    virtual bool CompareTo(const I3CLSimWlenDependentValue &other) const;
    
private:
    I3CLSimWlenDependentValueFromTable();
    
    double startWlen_;
    double wlenStep_;
    std::vector<double> values_;
    
    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive & ar, unsigned version);
};


BOOST_CLASS_VERSION(I3CLSimWlenDependentValueFromTable, i3clsimwlendependentvaluefromtable_version_);

I3_POINTER_TYPEDEFS(I3CLSimWlenDependentValueFromTable);

#endif //I3CLSIMWLENDEPENDENTVALUEFROMTABLE_H_INCLUDED
