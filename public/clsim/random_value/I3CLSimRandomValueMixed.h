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
 * @file I3CLSimRandomValueMixed.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMRANDOMVALUEMIXED_H_INCLUDED
#define I3CLSIMRANDOMVALUEMIXED_H_INCLUDED

#include "clsim/random_value/I3CLSimRandomValue.h"

/**
 * @brief A mix of two random values.
 */
static const unsigned i3clsimrandomvaluemixed_version_ = 0;

struct I3CLSimRandomValueMixed : public I3CLSimRandomValue
{
public:
    
    I3CLSimRandomValueMixed(double fractionOfFirstDistribution,
                            I3CLSimRandomValueConstPtr firstDistribution,
                            I3CLSimRandomValueConstPtr secondDistribution);
    virtual ~I3CLSimRandomValueMixed();

    virtual std::size_t NumberOfParameters() const;

    virtual double SampleFromDistribution(const I3RandomServicePtr &random,
                                          const std::vector<double> &parameters) const;

    virtual bool OpenCLFunctionWillOnlyUseASingleRandomNumber() const;

    virtual std::string GetOpenCLFunction(const std::string &functionName,
                                          const std::string &functionArgs,
                                          const std::string &functionArgsToCall,
                                          const std::string &uniformRandomCall_co,
                                          const std::string &uniformRandomCall_oc
                                          ) const;

    virtual bool CompareTo(const I3CLSimRandomValue &other) const;
    
private:
    I3CLSimRandomValueMixed();

    double fractionOfFirstDistribution_;
    I3CLSimRandomValueConstPtr firstDistribution_;
    I3CLSimRandomValueConstPtr secondDistribution_;
    
    friend class icecube::serialization::access;
    template <class Archive> void load(Archive & ar, unsigned version);
    template <class Archive> void save(Archive & ar, unsigned version) const;
    I3_SERIALIZATION_SPLIT_MEMBER();
};


I3_CLASS_VERSION(I3CLSimRandomValueMixed, i3clsimrandomvaluemixed_version_);

I3_POINTER_TYPEDEFS(I3CLSimRandomValueMixed);

#endif //I3CLSIMRANDOMVALUEMIXED_H_INCLUDED
