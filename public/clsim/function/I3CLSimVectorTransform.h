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
 * @file I3CLSimVectorTransform.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMVECTORTRANSFORM_H_INCLUDED
#define I3CLSIMVECTORTRANSFORM_H_INCLUDED

#include "icetray/serialization.h"
#include "icetray/I3TrayHeaders.h"

#include <string>
#include <vector>
#include <stdexcept>

/**
 * @brief Transforms a vector into another vector.
 * (works on 3-vectors)
 */
static const unsigned i3clsimvectortransform_version_ = 0;

struct I3CLSimVectorTransform 
{
public:
    
    I3CLSimVectorTransform();
    virtual ~I3CLSimVectorTransform();

    /**
     * if this is true, it is assumed that GetValue() returns a
     * meaningful value. If not, GetValue will not be called;
     * only the OpenCL implementation will be used.
     */
    virtual bool HasNativeImplementation() const = 0;

    /**
     * apply the transform
     */
    virtual std::vector<double> ApplyTransform(const std::vector<double> &vec) const = 0;

    /**
     * Return an OpenCL-compatible function named
     * functionName with a single float4* argument.
     * The function will transform the vector in place.
     * The value of the fourth element will be kept
     * un-changed.
     */
    virtual std::string GetOpenCLFunction(const std::string &functionName) const = 0;

    /**
     * compare to another I3CLSimVectorTransform object
     */
    virtual bool CompareTo(const I3CLSimVectorTransform &other) const = 0;
    
private:
    friend class icecube::serialization::access;
    template <class Archive> void serialize(Archive & ar, unsigned version);
};

inline bool operator==(const I3CLSimVectorTransform& a, const I3CLSimVectorTransform& b)
{
    return a.CompareTo(b);
}

inline bool operator!=(const I3CLSimVectorTransform& a, const I3CLSimVectorTransform& b)
{
    return (!a.CompareTo(b));
}


I3_CLASS_VERSION(I3CLSimVectorTransform, i3clsimvectortransform_version_);

I3_POINTER_TYPEDEFS(I3CLSimVectorTransform);

#endif //I3CLSIMVECTORTRANSFORM_H_INCLUDED
