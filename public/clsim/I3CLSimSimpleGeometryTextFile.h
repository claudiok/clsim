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
 * @file I3CLSimSimpleGeometryTextFile.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMSIMPLEGEOMETRYTEXTFILE_H_INCLUDED
#define I3CLSIMSIMPLEGEOMETRYTEXTFILE_H_INCLUDED

#include "clsim/I3CLSimSimpleGeometry.h"

#include <string>

/**
 * @brief Describes a detector geometry.
 *
 * Reads from a simple text file.
 */
class I3CLSimSimpleGeometryTextFile : public I3CLSimSimpleGeometry
{
public:
    static const int32_t default_ignoreStringIDsSmallerThan;
    static const int32_t default_ignoreStringIDsLargerThan;
    static const uint32_t default_ignoreDomIDsSmallerThan;
    static const uint32_t default_ignoreDomIDsLargerThan;
    
    I3CLSimSimpleGeometryTextFile(double OMRadius, 
                                  const std::string &filename,
                                  int32_t ignoreStringIDsSmallerThan=default_ignoreStringIDsSmallerThan,
                                  int32_t ignoreStringIDsLargerThan=default_ignoreStringIDsLargerThan,
                                  uint32_t ignoreDomIDsSmallerThan=default_ignoreDomIDsSmallerThan,
                                  uint32_t ignoreDomIDsLargerThan=default_ignoreDomIDsLargerThan
                                  );
    virtual ~I3CLSimSimpleGeometryTextFile();

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
    
private:
    double OMRadius_;
    std::size_t numOMs_;
    
    std::vector<int32_t> stringIDs_;
    std::vector<uint32_t> domIDs_;
    std::vector<double> posX_;
    std::vector<double> posY_;
    std::vector<double> posZ_;
    std::vector<std::string> subdetectors_;
    
    int32_t ignoreStringIDsSmallerThan_;
    int32_t ignoreStringIDsLargerThan_;
    uint32_t ignoreDomIDsSmallerThan_;
    uint32_t ignoreDomIDsLargerThan_;
    
};

I3_POINTER_TYPEDEFS(I3CLSimSimpleGeometryTextFile);

#endif //I3CLSIMSIMPLEGEOMETRYTEXTFILE_H_INCLUDED
