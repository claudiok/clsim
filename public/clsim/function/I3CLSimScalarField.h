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
 * @file I3CLSimScalarField.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMSCALARFIELD_H_INCLUDED
#define I3CLSIMSCALARFIELD_H_INCLUDED

#include "icetray/serialization.h"
#include "icetray/I3TrayHeaders.h"

#include <string>
#include <vector>
#include <stdexcept>

/**
 * @brief A scalar field (i.e. a function returning a scalar
 * for each point given by a 3-vector).
 */
static const unsigned i3clsimscalarfield_version_ = 0;

struct I3CLSimScalarField 
{
public:
    
    I3CLSimScalarField();
    virtual ~I3CLSimScalarField();

    /**
     * if this is true, it is assumed that GetValue() returns a
     * meaningful value. If not, GetValue will not be called;
     * only the OpenCL implementation will be used.
     */
    virtual bool HasNativeImplementation() const = 0;

    /**
     * return the value at a requested 3-vector
     */
    virtual double GetValue(double x, double y, double z) const = 0;

    /**
     * return the value at a requested 3-vector
     */
    inline double GetValue(const std::vector<double> &vec) const
    {
        if (vec.size() != 3)
            throw std::range_error("vector must contain excatly 3 elements!");
        return GetValue(vec[0], vec[1], vec[2]);
    }

    /**
     * return an OpenCL-compatible function named
     * functionName with a single float argument (float wlen)
     */
    virtual std::string GetOpenCLFunction(const std::string &functionName) const = 0;

    /**
     * compare to another I3CLSimScalarField object
     */
    virtual bool CompareTo(const I3CLSimScalarField &other) const = 0;
    
private:
    friend class icecube::serialization::access;
    template <class Archive> void serialize(Archive & ar, unsigned version);
};

inline bool operator==(const I3CLSimScalarField& a, const I3CLSimScalarField& b)
{
    return a.CompareTo(b);
}

inline bool operator!=(const I3CLSimScalarField& a, const I3CLSimScalarField& b)
{
    return (!a.CompareTo(b));
}


I3_CLASS_VERSION(I3CLSimScalarField, i3clsimscalarfield_version_);

I3_POINTER_TYPEDEFS(I3CLSimScalarField);

#endif //I3CLSIMSCALARFIELD_H_INCLUDED
