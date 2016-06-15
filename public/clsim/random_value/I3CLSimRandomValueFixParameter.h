/**
 * Copyright (c) 2012
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
 * @file I3CLSimRandomValueFixParameter.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMRANDOMVALUEFIXPARAMETER_H_INCLUDED
#define I3CLSIMRANDOMVALUEFIXPARAMETER_H_INCLUDED

#include "clsim/random_value/I3CLSimRandomValue.h"

#include <vector>
#include <string>

/**
 * @brief Takes a random distribution that depends on a number
 * of parameters and fixes one of them to a constant value.
 *
 * You could use this to fix the mean of the two-parameter
 * normal distribution while keeping the width a parameter.
 */
static const unsigned i3clsimrandomvaluefixparameter_version_ = 0;

struct I3CLSimRandomValueFixParameter : public I3CLSimRandomValue
{
public:
    
    I3CLSimRandomValueFixParameter(const I3CLSimRandomValuePtr &randomDistUsed,
                                   std::size_t parameterIndex,
                                   double parameterValue);

    virtual ~I3CLSimRandomValueFixParameter();

    virtual std::size_t NumberOfParameters() const;

    virtual double SampleFromDistribution(const I3RandomServicePtr &random,
                                          const std::vector<double> &parameters) const;

    virtual bool OpenCLFunctionWillOnlyUseASingleRandomNumber() const {return randomDistUsed_->OpenCLFunctionWillOnlyUseASingleRandomNumber();}

    virtual std::string GetOpenCLFunction(const std::string &functionName,
                                          const std::string &functionArgs,
                                          const std::string &functionArgsToCall,
                                          const std::string &uniformRandomCall_co,
                                          const std::string &uniformRandomCall_oc
                                          ) const;

    virtual bool CompareTo(const I3CLSimRandomValue &other) const;
    
private:
    I3CLSimRandomValueFixParameter();

    I3CLSimRandomValuePtr randomDistUsed_;
    std::size_t parameterIndex_;
    double parameterValue_;
    
    friend class icecube::serialization::access;
    template <class Archive> void serialize(Archive & ar, unsigned version);
};


I3_CLASS_VERSION(I3CLSimRandomValueFixParameter, i3clsimrandomvaluefixparameter_version_);

I3_POINTER_TYPEDEFS(I3CLSimRandomValueFixParameter);

#endif //I3CLSIMRANDOMVALUEFIXPARAMETER_H_INCLUDED
