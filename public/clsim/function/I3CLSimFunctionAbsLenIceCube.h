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
 * @file I3CLSimFunctionAbsLenIceCube.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMFUNCTIONABSLENICECUBE_H_INCLUDED
#define I3CLSIMFUNCTIONABSLENICECUBE_H_INCLUDED

#include "clsim/function/I3CLSimFunction.h"

#include <limits>

/**
 * @brief Six-parameter ice model as fitted in the SPICE ice model.
 * (Absorption length)
 */
static const unsigned i3clsimfunctionabslenicecube_version_ = 0;

struct I3CLSimFunctionAbsLenIceCube : public I3CLSimFunction
{
public:
    I3CLSimFunctionAbsLenIceCube(double kappa,
                                           double A,
                                           double B,
                                           double D,
                                           double E,
                                           double aDust400,
                                           double deltaTau
                                           );
    virtual ~I3CLSimFunctionAbsLenIceCube();
    
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
    
    
    // access to the internal state
    inline double GetKappa() const {return kappa_;}
    inline double GetA() const {return A_;}
    inline double GetB() const {return B_;}
    inline double GetD() const {return D_;}
    inline double GetE() const {return E_;}
    inline double GetADust400() const {return aDust400_;}
    inline double GetDeltaTau() const {return deltaTau_;}
    
private:
    I3CLSimFunctionAbsLenIceCube();
    
    double kappa_;
    double A_;
    double B_;
    double D_;
    double E_;
    double aDust400_;
    double deltaTau_;
    
    
    friend class icecube::serialization::access;
    template <class Archive> void serialize(Archive & ar, unsigned version);
};


I3_CLASS_VERSION(I3CLSimFunctionAbsLenIceCube, i3clsimfunctionabslenicecube_version_);

I3_POINTER_TYPEDEFS(I3CLSimFunctionAbsLenIceCube);

#endif //I3CLSIMFUNCTIONABSLENICECUBE_H_INCLUDED
