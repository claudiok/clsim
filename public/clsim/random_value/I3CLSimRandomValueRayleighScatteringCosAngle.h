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
 * @file I3CLSimRandomValueRayleighScatteringCosAngle.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMRANDOMVALUERAYLEIGHSCATTERINGCOSANGLE_H_INCLUDED
#define I3CLSIMRANDOMVALUERAYLEIGHSCATTERINGCOSANGLE_H_INCLUDED

#include "clsim/random_value/I3CLSimRandomValue.h"

/**
 * @brief The cosine of the scattering angle chosen according to
 * Rayleigh scattering.
 *
 */
static const unsigned i3clsimrandomvaluerayleighscatteringcosangle_version_ = 0;

struct I3CLSimRandomValueRayleighScatteringCosAngle : public I3CLSimRandomValue
{
public:
    
    I3CLSimRandomValueRayleighScatteringCosAngle();
    virtual ~I3CLSimRandomValueRayleighScatteringCosAngle();

    virtual std::size_t NumberOfParameters() const;

    virtual double SampleFromDistribution(const I3RandomServicePtr &random,
                                          const std::vector<double> &parameters) const;

    virtual bool OpenCLFunctionWillOnlyUseASingleRandomNumber() const {return true;}

    virtual std::string GetOpenCLFunction(const std::string &functionName,
                                          const std::string &functionArgs,
                                          const std::string &functionArgsToCall,
                                          const std::string &uniformRandomCall_co,
                                          const std::string &uniformRandomCall_oc
                                          ) const;

    virtual bool CompareTo(const I3CLSimRandomValue &other) const;
    
private:
    friend class icecube::serialization::access;
    template <class Archive> void serialize(Archive & ar, unsigned version);
};


I3_CLASS_VERSION(I3CLSimRandomValueRayleighScatteringCosAngle, i3clsimrandomvaluerayleighscatteringcosangle_version_);

I3_POINTER_TYPEDEFS(I3CLSimRandomValueRayleighScatteringCosAngle);

#endif //I3CLSIMRANDOMVALUERAYLEIGHSCATTERINGCOSANGLE_H_INCLUDED
