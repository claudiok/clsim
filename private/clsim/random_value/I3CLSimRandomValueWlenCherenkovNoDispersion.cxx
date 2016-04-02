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
 * @file I3CLSimRandomValueWlenCherenkovNoDispersion.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <icetray/serialization.h>
#include <icetray/I3Logging.h>
#include <clsim/random_value/I3CLSimRandomValueWlenCherenkovNoDispersion.h>

#include <boost/lexical_cast.hpp>

#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <limits>

#include "clsim/I3CLSimHelperToFloatString.h"

I3CLSimRandomValueWlenCherenkovNoDispersion::
I3CLSimRandomValueWlenCherenkovNoDispersion(double fromWlen,
                                            double toWlen)
:
fromWlen_(fromWlen),
toWlen_(toWlen)
{
    if (std::isnan(fromWlen_)) log_fatal("The \"fromWlen\" argument must not be NaN!");
    if (std::isnan(toWlen_)) log_fatal("The \"toWlen\" argument must not be NaN!");
              
    if (fromWlen_ > toWlen_) 
        log_fatal("The \"fromWlen\" argument must not be greater than \"toWlen\".");
}

I3CLSimRandomValueWlenCherenkovNoDispersion::~I3CLSimRandomValueWlenCherenkovNoDispersion() 
{ 
}

I3CLSimRandomValueWlenCherenkovNoDispersion::I3CLSimRandomValueWlenCherenkovNoDispersion() {;}

std::size_t I3CLSimRandomValueWlenCherenkovNoDispersion::NumberOfParameters() const {return 0;}

double I3CLSimRandomValueWlenCherenkovNoDispersion::SampleFromDistribution(const I3RandomServicePtr &random,
                                                                           const std::vector<double> &parameters) const
{
    if (!random) log_fatal("random service is NULL!");
    if (parameters.size() != 0) log_fatal("This distribution expects 0 parameters. Got %zu.", parameters.size());

    const double minVal = 1./toWlen_;
    const double range = (1./fromWlen_) - minVal;

    const double r = random->Uniform();
    return 1./(minVal + r * range);
}

std::string I3CLSimRandomValueWlenCherenkovNoDispersion::GetOpenCLFunction
(const std::string &functionName,
 const std::string &functionArgs,
 const std::string &functionArgsToCall,
 const std::string &uniformRandomCall_co,
 const std::string &uniformRandomCall_oc
 ) const
{
    const double minVal = 1./toWlen_;
    const double range = (1./fromWlen_) - minVal;

    const std::string functionDecl = std::string("inline float ") + functionName + "(" + functionArgs + ")";
    
    return functionDecl + ";\n\n" + functionDecl + "\n"
    "{\n"
    "    const float r = " + uniformRandomCall_oc + ";\n"
    "#ifdef USE_NATIVE_MATH\n"
    "    return native_recip(" + I3CLSimHelper::ToFloatString(minVal) + " + r * " + I3CLSimHelper::ToFloatString(range) + ");\n"
    "#else\n"
    "    return 1.f/(" + I3CLSimHelper::ToFloatString(minVal) + " + r * " + I3CLSimHelper::ToFloatString(range) + ");\n"
    "#endif\n"
    "}\n"
    ;
}

bool I3CLSimRandomValueWlenCherenkovNoDispersion::CompareTo(const I3CLSimRandomValue &other) const
{
    try
    {
        const I3CLSimRandomValueWlenCherenkovNoDispersion &other_ = dynamic_cast<const I3CLSimRandomValueWlenCherenkovNoDispersion &>(other);

        if (other_.fromWlen_ != fromWlen_) return false;
        if (other_.toWlen_ != toWlen_) return false;
        
        return true;
    }
    catch (const std::bad_cast& e)
    {
        // not of the same type, treat it as non-equal
        return false;
    }
    
}

template <class Archive>
void I3CLSimRandomValueWlenCherenkovNoDispersion::serialize(Archive &ar, unsigned version)
{
    if (version>i3clsimrandomvaluewlencherenkovnodispersion_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimRandomValueWlenCherenkovNoDispersion class.",
                  version,
                  i3clsimrandomvaluewlencherenkovnodispersion_version_);

    ar & make_nvp("I3CLSimRandomValue", base_object<I3CLSimRandomValue>(*this));
    ar & make_nvp("fromWlen", fromWlen_);
    ar & make_nvp("toWlen", toWlen_);
}     


I3_SERIALIZABLE(I3CLSimRandomValueWlenCherenkovNoDispersion);
