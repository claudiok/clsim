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
 * @file I3CLSimScalarFieldConstant.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMSCALARFIELDCONSTANT_H_INCLUDED
#define I3CLSIMSCALARFIELDCONSTANT_H_INCLUDED

#include "clsim/function/I3CLSimScalarField.h"

#include <vector>

/**
 * @brief A constant value.
 */
static const unsigned i3clsimscalarfieldconstant_version_ = 0;

struct I3CLSimScalarFieldConstant : public I3CLSimScalarField
{
public:
    
    I3CLSimScalarFieldConstant(double value);
    virtual ~I3CLSimScalarFieldConstant();
    
    /**
     * if this is true, it is assumed that GetValue() returns a
     * meaningful value. If not, GetValue will not be called;
     * only the OpenCL implementation will be used.
     */
    virtual bool HasNativeImplementation() const;
    
    /**
     * return the value at a requested 3-vector
     */
    virtual double GetValue(double x, double y, double z) const;

    /**
     * Shall return an OpenCL-compatible function named
     * functionName with a single float argument (float wlen)
     */
    virtual std::string GetOpenCLFunction(const std::string &functionName) const;

    /**
     * Shall compare to another I3CLSimScalarField object
     */
    virtual bool CompareTo(const I3CLSimScalarField &other) const;
    
private:
    I3CLSimScalarFieldConstant();
    
    double value_;
    
    friend class icecube::serialization::access;
    template <class Archive> void serialize(Archive & ar, unsigned version);
};


I3_CLASS_VERSION(I3CLSimScalarFieldConstant, i3clsimscalarfieldconstant_version_);

I3_POINTER_TYPEDEFS(I3CLSimScalarFieldConstant);

#endif //I3CLSIMSCALARFIELDCONSTANT_H_INCLUDED
