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
 * @file I3CLSimFunctionScatLenIceCube.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <icetray/serialization.h>
#include <icetray/I3Units.h>
#include <clsim/function/I3CLSimFunctionScatLenIceCube.h>

#include <typeinfo>
#include <cmath>

#include "clsim/I3CLSimHelperToFloatString.h"
using namespace I3CLSimHelper;

I3CLSimFunctionScatLenIceCube::
I3CLSimFunctionScatLenIceCube(double alpha,
                                        double b400
                                        )
:
alpha_(alpha),
b400_(b400)
{ 
}

I3CLSimFunctionScatLenIceCube::I3CLSimFunctionScatLenIceCube() {;}

I3CLSimFunctionScatLenIceCube::~I3CLSimFunctionScatLenIceCube() 
{;}


double I3CLSimFunctionScatLenIceCube::GetValue(double wlen) const
{
    const double x = wlen/I3Units::nanometer;
    return I3Units::m/( b400_ * std::pow(x/400., -alpha_) );
}


std::string I3CLSimFunctionScatLenIceCube::GetOpenCLFunction(const std::string &functionName) const
{
    std::string funcDef = 
    std::string("inline float ") + functionName + std::string("(float wlen)\n");

    const std::string refWlenAsString = ToFloatString(1./(400.*I3Units::nanometer));
    
    
    std::string funcBody = std::string() + 
    "{\n"
    "    const float alpha = " + ToFloatString(alpha_) + ";\n"
    "    const float b400 = " + ToFloatString(b400_) + ";\n"
    "    \n"
    "#ifdef USE_NATIVE_MATH\n"
    "    return " + ToFloatString(I3Units::m) + "*native_recip( b400 * native_powr(wlen*" + refWlenAsString + ", -alpha) );\n"
    "#else\n"
    "    return " + ToFloatString(I3Units::m) + "/( b400 * powr(wlen*" + refWlenAsString + ", -alpha) );\n"
    "#endif\n"
    "}\n"
    ;
    
    return funcDef + ";\n\n" + funcDef + funcBody;
}

bool I3CLSimFunctionScatLenIceCube::CompareTo(const I3CLSimFunction &other) const
{
    try
    {
        const I3CLSimFunctionScatLenIceCube &other_ = dynamic_cast<const I3CLSimFunctionScatLenIceCube &>(other);
        return ((other_.alpha_ == alpha_) &&
                (other_.b400_ == b400_));
    }
    catch (const std::bad_cast& e)
    {
        // not of the same type, treat it as non-equal
        return false;
    }
    
}


template <class Archive>
void I3CLSimFunctionScatLenIceCube::serialize(Archive &ar, unsigned version)
{
    if (version>i3clsimfunctionscatlenicecube_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimFunctionScatLenIceCube class.",version,i3clsimfunctionscatlenicecube_version_);

    ar & make_nvp("I3CLSimFunction", base_object<I3CLSimFunction>(*this));
    ar & make_nvp("alpha", alpha_);
    ar & make_nvp("b400", b400_);
}


I3_SERIALIZABLE(I3CLSimFunctionScatLenIceCube);
