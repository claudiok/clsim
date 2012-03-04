#ifndef I3CLSIMFUNCTIONFROMTABLE_H_INCLUDED
#define I3CLSIMFUNCTIONFROMTABLE_H_INCLUDED

#include "clsim/I3CLSimFunction.h"

#include <vector>

/**
 * @brief An arbitrary value interpolated from a uniformly distributed number
 * of values.
 */
static const unsigned i3clsimfunctionfromtable_version_ = 0;

struct I3CLSimFunctionFromTable : public I3CLSimFunction
{
public:
    
    // arbitrary wavelength values
    I3CLSimFunctionFromTable(const std::vector<double> &wlens,
                                       const std::vector<double> &values);
    
    // wavelengths with constant spacing (more efficient)
    I3CLSimFunctionFromTable(double startWlen,
                                       double wlenStep,
                                       const std::vector<double> &values);
    virtual ~I3CLSimFunctionFromTable();
    
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
     * Returns the internal state: number of bins
     */
    inline std::size_t GetNumEntries() const {return values_.size();}

    /**
     * Returns the internal state: value at entry i
     */
    inline double GetEntryValue(std::size_t i) const {return values_[i];}

    /**
     * Returns the internal state: wavenelgth at entry i
     */
    inline double GetEntryWavelength(std::size_t i) const {return wlens_[i];}

    /**
     * Returns the internal state: has equally spaced bins?
     */
    inline bool GetInEqualSpacingMode() const {return equalSpacingMode_;}

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
    I3CLSimFunctionFromTable();
    
    double startWlen_;
    double wlenStep_;
    std::vector<double> wlens_;
    std::vector<double> values_;
    
    bool equalSpacingMode_;
    
    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive & ar, unsigned version);
};


BOOST_CLASS_VERSION(I3CLSimFunctionFromTable, i3clsimfunctionfromtable_version_);

I3_POINTER_TYPEDEFS(I3CLSimFunctionFromTable);

#endif //I3CLSIMFUNCTIONFROMTABLE_H_INCLUDED
