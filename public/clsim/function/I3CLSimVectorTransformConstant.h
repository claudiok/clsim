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
 * @file I3CLSimVectorTransformConstant.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMVECTORTRANSFORMCONSTANT_H_INCLUDED
#define I3CLSIMVECTORTRANSFORMCONSTANT_H_INCLUDED

#include "clsim/function/I3CLSimVectorTransform.h"

#include <vector>

/**
 * @brief Will keep the input vector un-changed.
 */
static const unsigned i3clsimvectortransformconstant_version_ = 0;

struct I3CLSimVectorTransformConstant : public I3CLSimVectorTransform
{
public:
    
    I3CLSimVectorTransformConstant();
    virtual ~I3CLSimVectorTransformConstant();
    
    virtual bool HasNativeImplementation() const;

    virtual std::vector<double> ApplyTransform(const std::vector<double> &vec) const;

    /**
     * Return an OpenCL-compatible function named
     * functionName with a single float4* argument.
     * The function will transform the vector in place.
     * The value of the fourth element will be kept
     * un-changed.
     */
    virtual std::string GetOpenCLFunction(const std::string &functionName) const;

    /**
     * Shall compare to another I3CLSimVectorTransform object
     */
    virtual bool CompareTo(const I3CLSimVectorTransform &other) const;
    
private:
    
    friend class icecube::serialization::access;
    template <class Archive> void serialize(Archive & ar, unsigned version);
};


I3_CLASS_VERSION(I3CLSimVectorTransformConstant, i3clsimvectortransformconstant_version_);

I3_POINTER_TYPEDEFS(I3CLSimVectorTransformConstant);

#endif //I3CLSIMVECTORTRANSFORMCONSTANT_H_INCLUDED
