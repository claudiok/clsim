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
 * @file I3CLSimRandomValue.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMRANDOMVALUE_H_INCLUDED
#define I3CLSIMRANDOMVALUE_H_INCLUDED

#include "icetray/serialization.h"
#include "icetray/I3TrayHeaders.h"

#include "phys-services/I3RandomService.h"

#include <string>

/**
 * @brief A value chosen from a random distribution
 */
static const unsigned i3clsimrandomvalue_version_ = 0;

struct I3CLSimRandomValue 
{
public:
    
    I3CLSimRandomValue();
    virtual ~I3CLSimRandomValue();

    /**
     * Return a random number sampled from the distribution.
     * This runs as host code and is mainly for cross-checking
     * the OpenCL implementation.
     *
     * The parameters vector size needs to be the same
     * as the number returned by NumberOfParameters().
     */
    virtual double SampleFromDistribution(const I3RandomServicePtr &random,
                                          const std::vector<double> &parameters
                                         ) const = 0;

    /**
     * This should return the number of parameters this distribution
     * requires. For a gaussian this would be something like the
     * mean and sigma.
     */
    virtual std::size_t NumberOfParameters() const = 0;
    
    /**
     * If the OpenCL function will only use a single random number,
     * the random number can be passed directly as a value instead of
     * passing an random number generator. This may improve performance,
     * so it can be made known to the caller.
     */
    virtual bool OpenCLFunctionWillOnlyUseASingleRandomNumber() const = 0;
    
    /**
     * Shall return an OpenCL-compatible function.
     * The declaration of the form "float {functionName}({functionArgs})"
     * is provided by the caller in functionDecl.
     * The function call to generate uniformly distributed
     * random numbers between 0 and 1 is
     * provided in uniformRandomCall_{co|oc}.
     * (co: closed-open, 0 included, 1 not included
     *  oc: open-closed, 0 not included, 1 included)
     */
    virtual std::string GetOpenCLFunction(const std::string &functionName,
                                          const std::string &functionArgs,
                                          const std::string &functionArgsToCall,
                                          const std::string &uniformRandomCall_co,
                                          const std::string &uniformRandomCall_oc
                                          ) const = 0;

    /**
     * Shall compare the internal state another I3CLSimRandomValue object
     */
    virtual bool CompareTo(const I3CLSimRandomValue &other) const = 0;
    
private:
    friend class icecube::serialization::access;
    template <class Archive> void serialize(Archive & ar, unsigned version);
};

inline bool operator==(const I3CLSimRandomValue& a, const I3CLSimRandomValue& b)
{
    return a.CompareTo(b);
}


I3_CLASS_VERSION(I3CLSimRandomValue, i3clsimrandomvalue_version_);

I3_POINTER_TYPEDEFS(I3CLSimRandomValue);

#endif //I3CLSIMRANDOMVALUE_H_INCLUDED
