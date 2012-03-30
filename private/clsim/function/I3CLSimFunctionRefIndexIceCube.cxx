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
 * @file I3CLSimFunctionRefIndexIceCube.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <icetray/serialization.h>
#include <icetray/I3Units.h>
#include <clsim/function/I3CLSimFunctionRefIndexIceCube.h>

#include <typeinfo>
#include <cmath>

#include "clsim/I3CLSimHelperToFloatString.h"
using namespace I3CLSimHelper;

const std::string I3CLSimFunctionRefIndexIceCube::default_mode = "phase";
const double I3CLSimFunctionRefIndexIceCube::default_n0= 1.55749;
const double I3CLSimFunctionRefIndexIceCube::default_n1=-1.57988;
const double I3CLSimFunctionRefIndexIceCube::default_n2= 3.99993;
const double I3CLSimFunctionRefIndexIceCube::default_n3=-4.68271;
const double I3CLSimFunctionRefIndexIceCube::default_n4= 2.09354;
const double I3CLSimFunctionRefIndexIceCube::default_g0= 1.227106;
const double I3CLSimFunctionRefIndexIceCube::default_g1=-0.954648;
const double I3CLSimFunctionRefIndexIceCube::default_g2= 1.42568;
const double I3CLSimFunctionRefIndexIceCube::default_g3=-0.711832;
const double I3CLSimFunctionRefIndexIceCube::default_g4= 0.00000;


I3CLSimFunctionRefIndexIceCube::
I3CLSimFunctionRefIndexIceCube(std::string mode,
                                         double n0,
                                         double n1,
                                         double n2,
                                         double n3,
                                         double n4,
                                         double g0,
                                         double g1,
                                         double g2,
                                         double g3,
                                         double g4
                                         )
:
mode_(mode),
n0_(n0),
n1_(n1),
n2_(n2),
n3_(n3),
n4_(n4),
g0_(g0),
g1_(g1),
g2_(g2),
g3_(g3),
g4_(g4)
{ 
    if ((mode_!="phase") && (mode_!="group"))
        log_fatal("In I3CLSimFunctionRefIndexIceCube(): \"mode\" argument has to be either \"phase\" or \"group\".");
}

I3CLSimFunctionRefIndexIceCube::~I3CLSimFunctionRefIndexIceCube() 
{;}


double I3CLSimFunctionRefIndexIceCube::GetValue(double wlen) const
{
    const double x = wlen/I3Units::micrometer;
    const double np = n0_ + x*(n1_ + x*(n2_ + x*(n3_ + x*n4_)));
    
    if (mode_=="phase")
    {
        return np;
    }
    else if (mode_=="group")
    {
        const double np_corr = g0_ + x*(g1_ + x*(g2_ + x*(g3_ + x*g4_)));
        return np * np_corr;
    }
    else
    {
        log_fatal("Internal error. Invalid mode.");
    }
}

double I3CLSimFunctionRefIndexIceCube::GetDerivative(double wlen) const
{
    const double x = wlen/I3Units::micrometer;

    const double dnp = (n1_ + x*(2.*n2_ + x*(3.*n3_ + x*4.*n4_)))/I3Units::micrometer;

    if (mode_=="phase")
    {
        return dnp;
    }
    else if (mode_=="group")
    {
        const double np = n0_ + x*(n1_ + x*(n2_ + x*(n3_ + x*n4_)));
        const double np_corr = g0_ + x*(g1_ + x*(g2_ + x*(g3_ + x*g4_)));
        const double dnp_corr = (g1_ + x*(2.*g2_ + x*(3.*g3_ + x*4.*g4_)))/I3Units::micrometer;

        return dnp*np_corr + np*dnp_corr;
    }
    else
    {
        log_fatal("Internal error. Invalid mode.");
    }
}

std::string I3CLSimFunctionRefIndexIceCube::GetOpenCLFunction(const std::string &functionName) const
{
    std::string funcDef = 
    std::string("inline float ") + functionName + std::string("(float wlen)\n");

    std::string funcBody = std::string() + 
    "{\n"
    "    const float n0 = " + ToFloatString(n0_) + ";\n"
    "    const float n1 = " + ToFloatString(n1_) + ";\n"
    "    const float n2 = " + ToFloatString(n2_) + ";\n"
    "    const float n3 = " + ToFloatString(n3_) + ";\n"
    "    const float n4 = " + ToFloatString(n4_) + ";\n";
    
    if (mode_=="group")
    {
        funcBody = funcBody +
        "    \n"
        "    const float g0 = " + ToFloatString(g0_) + ";\n"
        "    const float g1 = " + ToFloatString(g1_) + ";\n"
        "    const float g2 = " + ToFloatString(g2_) + ";\n"
        "    const float g3 = " + ToFloatString(g3_) + ";\n"
        "    const float g4 = " + ToFloatString(g4_) + ";\n";
    }
    
    funcBody = funcBody +
    "    \n"
    "    const float x = wlen/" + ToFloatString(I3Units::micrometer) + ";\n"
    "    const float np = n0 + x*(n1 + x*(n2 + x*(n3 + x*n4)));\n";

    if (mode_=="phase")
    {
        funcBody = funcBody +
        "    \n"
        "    return np;\n";
    }
    else if (mode_=="group")
    {
        funcBody = funcBody +
        "    const float np_corr = g0 + x*(g1 + x*(g2 + x*(g3 + x*g4)));\n"
        "    \n"
        "    return np*np_corr;\n";
    }
    else
    {
        log_fatal("Internal error. Invalid mode.");
    }

    funcBody = funcBody +
    "}\n"
    ;
    
    return funcDef + ";\n\n" + funcDef + funcBody;
}

std::string I3CLSimFunctionRefIndexIceCube::GetOpenCLFunctionDerivative(const std::string &functionName) const
{
    std::string funcDef = 
    std::string("inline float ") + functionName + std::string("(float wlen)\n");
    
    std::string funcBody = std::string() + 
    "{\n";

    if (mode_=="group")
    {
        funcBody = funcBody +
        "    const float n0 = " + ToFloatString(n0_) + ";\n";
    }
    
    funcBody = funcBody +
    "    const float n1 = " + ToFloatString(n1_) + ";\n"
    "    const float n2 = " + ToFloatString(n2_) + ";\n"
    "    const float n3 = " + ToFloatString(n3_) + ";\n"
    "    const float n4 = " + ToFloatString(n4_) + ";\n";
    
    if (mode_=="group")
    {
        funcBody = funcBody +
        "    \n"
        "    const float g0 = " + ToFloatString(g0_) + ";\n"
        "    const float g1 = " + ToFloatString(g1_) + ";\n"
        "    const float g2 = " + ToFloatString(g2_) + ";\n"
        "    const float g3 = " + ToFloatString(g3_) + ";\n"
        "    const float g4 = " + ToFloatString(g4_) + ";\n";
    }

    funcBody = funcBody +
    "    \n"
    "    const float x = wlen/" + ToFloatString(I3Units::micrometer) + ";\n"
    "    \n"
    "    const float dnp = (n1 + x*(2.f*n2 + x*(3.f*n3 + x*4.f*n4)))/" + ToFloatString(I3Units::micrometer) + ";\n";
    
    if (mode_=="phase")
    {
        funcBody = funcBody +
        "    \n"
        "    return dnp;\n";
    }
    else if (mode_=="group")
    {
        funcBody = funcBody +
        "    const float np = n0 + x*(n1 + x*(n2 + x*(n3 + x*n4)));\n"
        "    const float np_corr = g0 + x*(g1 + x*(g2 + x*(g3 + x*g4)));\n"
        "    const float dnp_corr = (g1 + x*(2.f*g2 + x*(3.f*g3 + x*4.f*g4)))/" + ToFloatString(I3Units::micrometer) + ";\n"
        "    \n"
        "    return dnp*np_corr + np*dnp_corr;\n";
    }
    else
    {
        log_fatal("Internal error. Invalid mode.");
    }
    
    funcBody = funcBody +
    "}\n"
    ;
    
    return funcDef + ";\n\n" + funcDef + funcBody;
}

bool I3CLSimFunctionRefIndexIceCube::CompareTo(const I3CLSimFunction &other) const
{
    try
    {
        const I3CLSimFunctionRefIndexIceCube &other_ = dynamic_cast<const I3CLSimFunctionRefIndexIceCube &>(other);
        return ((other_.n0_ == n0_) &&
                (other_.n1_ == n1_) &&
                (other_.n2_ == n2_) &&
                (other_.n3_ == n3_) &&
                (other_.n4_ == n4_) &&
                (other_.g0_ == g0_) &&
                (other_.g1_ == g1_) &&
                (other_.g2_ == g2_) &&
                (other_.g3_ == g3_) &&
                (other_.g4_ == g4_) &&
                (other_.mode_ == mode_));
    }
    catch (const std::bad_cast& e)
    {
        // not of the same type, treat it as non-equal
        return false;
    }
    
}


template <class Archive>
void I3CLSimFunctionRefIndexIceCube::serialize(Archive &ar, unsigned version)
{
    if (version>i3clsimfunctionrefindexicecube_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimFunctionRefIndexIceCube class.",version,i3clsimfunctionrefindexicecube_version_);

    ar & make_nvp("I3CLSimFunction", base_object<I3CLSimFunction>(*this));
    ar & make_nvp("mode", mode_);
    ar & make_nvp("n0", n0_);
    ar & make_nvp("n1", n1_);
    ar & make_nvp("n2", n2_);
    ar & make_nvp("n3", n3_);
    ar & make_nvp("n4", n4_);
    ar & make_nvp("g0", g0_);
    ar & make_nvp("g1", g1_);
    ar & make_nvp("g2", g2_);
    ar & make_nvp("g3", g3_);
    ar & make_nvp("g4", g4_);
}


I3_SERIALIZABLE(I3CLSimFunctionRefIndexIceCube);
