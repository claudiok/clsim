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
 * @file I3CLSimFunctionDeltaPeak.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <icetray/serialization.h>
#include <clsim/function/I3CLSimFunctionDeltaPeak.h>

#include <typeinfo>
#include <cmath>
#include <math.h>
#include <sstream>
#include <stdexcept>
#include <limits>

#include "clsim/I3CLSimHelperToFloatString.h"
using namespace I3CLSimHelper;

I3CLSimFunctionDeltaPeak::
I3CLSimFunctionDeltaPeak(double peakPosition)
:
peakPosition_(peakPosition)
{
}

I3CLSimFunctionDeltaPeak::I3CLSimFunctionDeltaPeak() {;}

I3CLSimFunctionDeltaPeak::~I3CLSimFunctionDeltaPeak() 
{;}

bool I3CLSimFunctionDeltaPeak::HasNativeImplementation() const 
{
    return true;
}

bool I3CLSimFunctionDeltaPeak::HasDerivative() const 
{
    return false;
}

double I3CLSimFunctionDeltaPeak::GetValue(double wlen) const
{
    if (wlen==peakPosition_) {
        return std::numeric_limits<double>::infinity();
    } else {
        return 0.;
    }
}

double I3CLSimFunctionDeltaPeak::GetMinWlen() const
{
    return -std::numeric_limits<double>::infinity();
}

double I3CLSimFunctionDeltaPeak::GetMaxWlen() const
{
    return std::numeric_limits<double>::infinity();
}

std::string I3CLSimFunctionDeltaPeak::GetOpenCLFunction(const std::string &functionName) const
{
    std::string funcDef = 
    std::string("inline float ") + functionName + std::string("(float wavelength)\n");
    
    std::ostringstream output(std::ostringstream::out);
    output.setf(std::ios::scientific,std::ios::floatfield);
    output.precision(std::numeric_limits<float>::digits10+4); // maximum precision for a float
    
    output << peakPosition_ << "f";
    std::string peakPositionString = output.str();
    
    std::string funcBody = std::string() + 
    "if (wlen==" + peakPositionString + ") {\n"
    "    return INFINITY;\n"
    "} else {\n"
    "    return 0.;\n"
    "}\n"
    ;
    
    return funcDef + ";\n\n" + funcDef + funcBody;
}

bool I3CLSimFunctionDeltaPeak::CompareTo(const I3CLSimFunction &other) const
{
    try
    {
        const I3CLSimFunctionDeltaPeak &other_ = dynamic_cast<const I3CLSimFunctionDeltaPeak &>(other);
        if ((other_.peakPosition_ != peakPosition_))
            return false;
        
        return true;
    }
    catch (const std::bad_cast& e)
    {
        // not of the same type, treat it as non-equal
        return false;
    }
    
}



template <class Archive>
void I3CLSimFunctionDeltaPeak::serialize(Archive &ar, unsigned version)
{
    if (version>i3clsimfunctiondeltapeak_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimFunctionDeltaPeak class.",version,i3clsimfunctiondeltapeak_version_);

    ar & make_nvp("I3CLSimFunction", base_object<I3CLSimFunction>(*this));
    ar & make_nvp("peakPosition", peakPosition_);
}     


I3_SERIALIZABLE(I3CLSimFunctionDeltaPeak);
