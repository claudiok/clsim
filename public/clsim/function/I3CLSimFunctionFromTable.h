/**
 * Copyright (c) 2011, 2012
 * Claudio Kopper <claudio.kopper@icecube.wisc.edu>
 * and the IceCube Collaboration <http://www.icecube.wisc.edu>
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
 * SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
 * OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
 * CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 *
 * $Id$
 *
 * @file I3CLSimFunctionFromTable.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMFUNCTIONFROMTABLE_H_INCLUDED
#define I3CLSIMFUNCTIONFROMTABLE_H_INCLUDED

#include "clsim/function/I3CLSimFunction.h"

#include <vector>

/**
 * @brief An arbitrary value interpolated from a uniformly distributed number
 * of values.
 */
static const unsigned i3clsimfunctionfromtable_version_ = 0;

struct I3CLSimFunctionFromTable : public I3CLSimFunction
{
public:
    static const bool default_storeDataAsHalfPrecision;
    
    // arbitrary wavelength values
    I3CLSimFunctionFromTable(const std::vector<double> &wlens,
                             const std::vector<double> &values,
                             bool storeDataAsHalfPrecision=default_storeDataAsHalfPrecision);
    
    // wavelengths with constant spacing (more efficient)
    I3CLSimFunctionFromTable(double startWlen,
                             double wlenStep,
                             const std::vector<double> &values,
                             bool storeDataAsHalfPrecision=default_storeDataAsHalfPrecision);
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
    
    bool storeDataAsHalfPrecision_;
    
    friend class icecube::serialization::access;
    template <class Archive> void serialize(Archive & ar, unsigned version);
};


I3_CLASS_VERSION(I3CLSimFunctionFromTable, i3clsimfunctionfromtable_version_);

I3_POINTER_TYPEDEFS(I3CLSimFunctionFromTable);

#endif //I3CLSIMFUNCTIONFROMTABLE_H_INCLUDED
