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
 * @file I3CLSimRandomValueInterpolatedDistribution.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMRANDOMVALUEINTERPOLATEDDISTRIBUTION_H_INCLUDED
#define I3CLSIMRANDOMVALUEINTERPOLATEDDISTRIBUTION_H_INCLUDED

#include "clsim/random_value/I3CLSimRandomValue.h"

#include <vector>

/**
 * @brief A random value chosen according to a given distribution.
 * The distribution is linearly interpolated between the given data
 * points.
 */
static const unsigned i3clsimrandomvalueinterpolateddistribution_version_ = 0;

struct I3CLSimRandomValueInterpolatedDistribution : public I3CLSimRandomValue
{
public:
    
    // arbitrary x values
    I3CLSimRandomValueInterpolatedDistribution(const std::vector<double> &x,
                                               const std::vector<double> &y);

    // x values with constant spacing (more efficient)
    I3CLSimRandomValueInterpolatedDistribution(double xFirst, double xSpacing,
                                               const std::vector<double> &y);

    virtual ~I3CLSimRandomValueInterpolatedDistribution();

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
    void InitTables();
    std::string WriteTableCode(const std::string &prefix) const;
    
    I3CLSimRandomValueInterpolatedDistribution();

    std::vector<double> data_acu_;
    std::vector<double> data_beta_;

    std::vector<double> x_;
    std::vector<double> y_;
    double constantXSpacing_;
    double firstX_;
    
    friend class icecube::serialization::access;
    template <class Archive> void serialize(Archive & ar, unsigned version);
};


I3_CLASS_VERSION(I3CLSimRandomValueInterpolatedDistribution, i3clsimrandomvalueinterpolateddistribution_version_);

I3_POINTER_TYPEDEFS(I3CLSimRandomValueInterpolatedDistribution);

#endif //I3CLSIMRANDOMVALUEINTERPOLATEDDISTRIBUTION_H_INCLUDED
