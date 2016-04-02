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
 * @file I3CLSimSimpleGeometry.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMSIMPLEGEOMETRY_H_INCLUDED
#define I3CLSIMSIMPLEGEOMETRY_H_INCLUDED

#include "icetray/I3TrayHeaders.h"

#include <vector>
#include <string>

/**
 * @brief Describes a detector geometry (abstract base class)
 */

class I3CLSimSimpleGeometry 
{
    
public:
    virtual std::size_t size() const = 0;
    
    /// This is the radius *with* oversizing applied!
    virtual double GetOMRadius() const = 0;

    virtual const std::vector<int32_t> &GetStringIDVector() const = 0;
    virtual const std::vector<uint32_t> &GetDomIDVector() const = 0;
    virtual const std::vector<double> &GetPosXVector() const = 0;
    virtual const std::vector<double> &GetPosYVector() const = 0;
    virtual const std::vector<double> &GetPosZVector() const = 0;
    virtual const std::vector<std::string> &GetSubdetectorVector() const = 0;
    
    virtual int32_t GetStringID(std::size_t pos) const = 0;
    virtual uint32_t GetDomID(std::size_t pos) const = 0;
    virtual double GetPosX(std::size_t pos) const = 0;
    virtual double GetPosY(std::size_t pos) const = 0;
    virtual double GetPosZ(std::size_t pos) const = 0;
    virtual std::string GetSubdetector(std::size_t pos) const = 0;
    
    
};

I3_POINTER_TYPEDEFS(I3CLSimSimpleGeometry);

#endif //I3CLSIMSIMPLEGEOMETRY_H_INCLUDED
