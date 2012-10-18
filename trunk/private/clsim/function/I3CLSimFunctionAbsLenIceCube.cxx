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
 * @file I3CLSimFunctionAbsLenIceCube.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <icetray/serialization.h>
#include <icetray/I3Units.h>
#include <clsim/function/I3CLSimFunctionAbsLenIceCube.h>

#include <typeinfo>
#include <cmath>

#include "clsim/I3CLSimHelperToFloatString.h"
using namespace I3CLSimHelper;

I3CLSimFunctionAbsLenIceCube::
I3CLSimFunctionAbsLenIceCube(double kappa,
                                       double A,
                                       double B,
                                       double D,
                                       double E,
                                       double aDust400,
                                       double deltaTau
                                       )
:
kappa_(kappa),
A_(A),
B_(B),
D_(D),
E_(E),
aDust400_(aDust400),
deltaTau_(deltaTau)
{ 
}

I3CLSimFunctionAbsLenIceCube::I3CLSimFunctionAbsLenIceCube() {;}

I3CLSimFunctionAbsLenIceCube::~I3CLSimFunctionAbsLenIceCube() 
{;}


double I3CLSimFunctionAbsLenIceCube::GetValue(double wlen) const
{
    const double x = wlen/I3Units::nanometer;
    return I3Units::m/( (D_*aDust400_+E_) * std::pow(x, -kappa_)  +  A_*std::exp(-B_/x) * (1.f + 0.01f*deltaTau_) );
}


std::string I3CLSimFunctionAbsLenIceCube::GetOpenCLFunction(const std::string &functionName) const
{
    std::string funcDef = 
    std::string("inline float ") + functionName + std::string("(float wlen)\n");

    std::string funcBody = std::string() + 
    "{\n"
    "    const float kappa = " + ToFloatString(kappa_) + ";\n"
    "    const float A = " + ToFloatString(A_) + ";\n"
    "    const float B = " + ToFloatString(B_) + ";\n"
    "    const float D = " + ToFloatString(D_) + ";\n"
    "    const float E = " + ToFloatString(E_) + ";\n"
    "    const float aDust400 = " + ToFloatString(aDust400_) + ";\n"
    "    const float deltaTau = " + ToFloatString(deltaTau_) + ";\n"
    "    \n"
    "    const float x = wlen/" + ToFloatString(I3Units::nanometer) + ";\n"
    "    \n"
    "#ifdef USE_NATIVE_MATH\n"
    "    return " + ToFloatString(I3Units::m) + "*native_recip( (D*aDust400+E) * native_powr(x, -kappa)  +  A*native_exp(-B/x) * (1.f + 0.01f*deltaTau) );\n"
    "#else\n"
    "    return " + ToFloatString(I3Units::m) + "/( (D*aDust400+E) * powr(x, -kappa)  +  A*exp(-B/x) * (1.f + 0.01f*deltaTau) );\n"
    "#endif\n"
    "}\n"
    ;
    
    return funcDef + ";\n\n" + funcDef + funcBody;
}

bool I3CLSimFunctionAbsLenIceCube::CompareTo(const I3CLSimFunction &other) const
{
    try
    {
        const I3CLSimFunctionAbsLenIceCube &other_ = dynamic_cast<const I3CLSimFunctionAbsLenIceCube &>(other);
        return ((other_.kappa_ == kappa_) &&
                (other_.A_ == A_) &&
                (other_.B_ == B_) &&
                (other_.D_ == D_) &&
                (other_.E_ == E_) &&
                (other_.aDust400_ == aDust400_) &&
                (other_.deltaTau_ == deltaTau_));
    }
    catch (const std::bad_cast& e)
    {
        // not of the same type, treat it as non-equal
        return false;
    }
    
}


template <class Archive>
void I3CLSimFunctionAbsLenIceCube::serialize(Archive &ar, unsigned version)
{
    if (version>i3clsimfunctionabslenicecube_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimFunctionAbsLenIceCube class.",version,i3clsimfunctionabslenicecube_version_);

    ar & make_nvp("I3CLSimFunction", base_object<I3CLSimFunction>(*this));
    ar & make_nvp("kappa", kappa_);
    ar & make_nvp("A", A_);
    ar & make_nvp("B", B_);
    ar & make_nvp("D", D_);
    ar & make_nvp("E", E_);
    ar & make_nvp("aDust400", aDust400_);
    ar & make_nvp("deltaTau", deltaTau_);
}


I3_SERIALIZABLE(I3CLSimFunctionAbsLenIceCube);
