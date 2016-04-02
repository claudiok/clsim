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
 * @file I3CLSimSimpleGeometryUserConfigurable.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMSIMPLEGEOMETRYUSERCONFIGURABLE_H_INCLUDED
#define I3CLSIMSIMPLEGEOMETRYUSERCONFIGURABLE_H_INCLUDED

#include "clsim/I3CLSimSimpleGeometry.h"

/**
 * @brief Describes a detector geometry.
 *
 * The DOM properties have to be explicitly set by the user.
 */

class I3CLSimSimpleGeometryUserConfigurable : public I3CLSimSimpleGeometry
{
    
public:
    I3CLSimSimpleGeometryUserConfigurable(double OMRadius, std::size_t numOMs);
    virtual ~I3CLSimSimpleGeometryUserConfigurable();

    virtual std::size_t size() const {return numOMs_;}

    /// This is the radius *with* oversizing applied!
    virtual double GetOMRadius() const {return OMRadius_;}
    
    virtual const std::vector<int32_t> &GetStringIDVector() const {return stringIDs_;}
    virtual const std::vector<uint32_t> &GetDomIDVector() const {return domIDs_;}
    virtual const std::vector<double> &GetPosXVector() const {return posX_;}
    virtual const std::vector<double> &GetPosYVector() const {return posY_;}
    virtual const std::vector<double> &GetPosZVector() const {return posZ_;}
    virtual const std::vector<std::string> &GetSubdetectorVector() const {return subdetectors_;}

    virtual int32_t GetStringID(std::size_t pos) const {return stringIDs_.at(pos);}
    virtual uint32_t GetDomID(std::size_t pos) const {return domIDs_.at(pos);}
    virtual double GetPosX(std::size_t pos) const {return posX_.at(pos);}
    virtual double GetPosY(std::size_t pos) const {return posY_.at(pos);}
    virtual double GetPosZ(std::size_t pos) const {return posZ_.at(pos);}
    virtual std::string GetSubdetector(std::size_t pos) const {return subdetectors_.at(pos);}
    
    
    // additional methods
    virtual void SetStringID(std::size_t pos, int32_t val) {stringIDs_.at(pos) = val;}
    virtual void SetDomID(std::size_t pos, uint32_t val) {domIDs_.at(pos) = val;}
    virtual void SetPosX(std::size_t pos, double val) {posX_.at(pos) = val;}
    virtual void SetPosY(std::size_t pos, double val) {posY_.at(pos) = val;}
    virtual void SetPosZ(std::size_t pos, double val) {posZ_.at(pos) = val;}
    virtual void SetSubdetector(std::size_t pos, const std::string &val) {subdetectors_.at(pos) = val;}

    
private:
    double OMRadius_;
    std::size_t numOMs_;
    
    std::vector<int32_t> stringIDs_;
    std::vector<uint32_t> domIDs_;
    std::vector<double> posX_;
    std::vector<double> posY_;
    std::vector<double> posZ_;
    std::vector<std::string> subdetectors_;
};

I3_POINTER_TYPEDEFS(I3CLSimSimpleGeometryUserConfigurable);

#endif //I3CLSIMSIMPLEGEOMETRYUSERCONFIGURABLE_H_INCLUDED
