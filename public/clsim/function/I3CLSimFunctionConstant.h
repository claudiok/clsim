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
 * @file I3CLSimFunctionConstant.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMFUNCTIONCONSTANT_H_INCLUDED
#define I3CLSIMFUNCTIONCONSTANT_H_INCLUDED

#include "clsim/function/I3CLSimFunction.h"

#include <vector>

/**
 * @brief A constant value.
 */
static const unsigned i3clsimfunctionconstant_version_ = 0;

struct I3CLSimFunctionConstant : public I3CLSimFunction
{
public:
    
    I3CLSimFunctionConstant(double value);
    virtual ~I3CLSimFunctionConstant();
    
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
     * Shall compare to another I3CLSimFunction object
     */
    virtual bool CompareTo(const I3CLSimFunction &other) const;
    
private:
    I3CLSimFunctionConstant();
    
    double value_;
    
    friend class icecube::serialization::access;
    template <class Archive> void serialize(Archive & ar, unsigned version);
};


I3_CLASS_VERSION(I3CLSimFunctionConstant, i3clsimfunctionconstant_version_);

I3_POINTER_TYPEDEFS(I3CLSimFunctionConstant);

#endif //I3CLSIMFUNCTIONCONSTANT_H_INCLUDED
