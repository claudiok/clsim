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
 * @file I3CLSimScalarFieldAnisotropyAbsLenScaling.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <icetray/serialization.h>
#include <clsim/function/I3CLSimScalarFieldAnisotropyAbsLenScaling.h>
#include <icetray/I3Units.h>

#include <typeinfo>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include "clsim/I3CLSimHelperToFloatString.h"
using namespace I3CLSimHelper;

const double I3CLSimScalarFieldAnisotropyAbsLenScaling::default_anisotropyDirAzimuth = 216.*I3Units::deg;
const double I3CLSimScalarFieldAnisotropyAbsLenScaling::default_magnitudeAlongDir = 0.04;
const double I3CLSimScalarFieldAnisotropyAbsLenScaling::default_magnitudePerpToDir = -0.08;

I3CLSimScalarFieldAnisotropyAbsLenScaling::
I3CLSimScalarFieldAnisotropyAbsLenScaling(
    double anisotropyDirAzimuth,
    double magnitudeAlongDir,
    double magnitudePerpToDir)
:
anisotropyDirAzimuth_(anisotropyDirAzimuth),
magnitudeAlongDir_(magnitudeAlongDir),
magnitudePerpToDir_(magnitudePerpToDir)
{
}

I3CLSimScalarFieldAnisotropyAbsLenScaling::~I3CLSimScalarFieldAnisotropyAbsLenScaling() 
{;}

bool I3CLSimScalarFieldAnisotropyAbsLenScaling::HasNativeImplementation() const 
{
    return true;
}

double I3CLSimScalarFieldAnisotropyAbsLenScaling::GetValue(double x, double y, double z) const
{
    const double azx = std::cos(anisotropyDirAzimuth_); // x-component of the "tilt" direction
    const double azy = std::sin(anisotropyDirAzimuth_); // y-component of the "tilt" direction
    const double k1  = std::exp(magnitudeAlongDir_);    // coefficient of anisotropy parallel to "tilt" direction
    const double k2  = std::exp(magnitudePerpToDir_);   // coefficient of anisotropy perpendicular to "tilt" direction
    const double kz  = 1./(k1*k2);                      // a normalizing factor for the z direction
    const double l1=k1*k1;
    const double l2=k2*k2;
    const double l3=kz*kz;
    const double B2=1./l1+1./l2+1./l3;


    const double n1= azx*x+azy*y;
    const double n2=-azy*x+azx*y;
    const double n3= z;

    const double s1=n1*n1;
    const double s2=n2*n2;
    const double s3=n3*n3;

    const double nB=s1/l1+s2/l2+s3/l3;
    const double An=s1*l1+s2*l2+s3*l3;

    const double nr=(B2-nB)*An/2.;

    return 1./nr; // the absorption length is multiplied with this factor in the kernel
}

std::string I3CLSimScalarFieldAnisotropyAbsLenScaling::GetOpenCLFunction(const std::string &functionName) const
{
    // the OpenCL interface takes a float4, but ignores the fourth component

    // pre-calculate constants used in the function
    const double azx = std::cos(anisotropyDirAzimuth_); // x-component of the "tilt" direction
    const double azy = std::sin(anisotropyDirAzimuth_); // y-component of the "tilt" direction
    const double k1  = std::exp(magnitudeAlongDir_);    // coefficient of anisotropy parallel to "tilt" direction
    const double k2  = std::exp(magnitudePerpToDir_);   // coefficient of anisotropy perpendicular to "tilt" direction
    const double kz  = 1./(k1*k2);                      // a normalizing factor for the z direction

    const double l1  = k1*k1;
    const double l2  = k2*k2;
    const double l3  = kz*kz;
    const double B2  = 1./l1+1./l2+1./l3;

    std::string funcDef = 
        std::string("inline float ") + functionName + std::string("(float4 vec)");
    
    std::ostringstream output(std::ostringstream::out);
    
    std::string funcBody = std::string() + 
    "{\n"
    "    const float4 l  = (float4)(" + ToFloatString(l1)    + ", " + ToFloatString(l2)    + ", " + ToFloatString(l3)    + ", 0.f);\n"
    "    const float4 rl = (float4)(" + ToFloatString(1./l1) + ", " + ToFloatString(1./l2) + ", " + ToFloatString(1./l3) + ", 0.f);\n"
    "    \n"
    "    const float4 n = (float4)\n"
    "        (\n"
    "         (" + ToFloatString( azx) + "*vec.x)+(" + ToFloatString(azy) + "*vec.y),\n"
    "         (" + ToFloatString(-azy) + "*vec.x)+(" + ToFloatString(azx) + "*vec.y),\n"
    "         vec.z,\n"
    "         0.f\n"
    "        );\n"
    "    const float4 s=n*n;\n"
    "    \n"
    "    const float nB = dot(s,rl);\n"
    "    const float An = dot(s,l);\n"
    "    \n"
    "    return 2.f/((" + ToFloatString(B2) + "-nB)*An);\n"
    "}\n"
    ;

    return funcDef + ";\n\n" + funcDef + "\n" + funcBody;
}

bool I3CLSimScalarFieldAnisotropyAbsLenScaling::CompareTo(const I3CLSimScalarField &other) const
{
    try
    {
        const I3CLSimScalarFieldAnisotropyAbsLenScaling &other_ = dynamic_cast<const I3CLSimScalarFieldAnisotropyAbsLenScaling &>(other);
        if ((other_.anisotropyDirAzimuth_ != anisotropyDirAzimuth_) ||
            (other_.magnitudeAlongDir_ != magnitudeAlongDir_) ||
            (other_.magnitudePerpToDir_ != magnitudePerpToDir_))
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
void I3CLSimScalarFieldAnisotropyAbsLenScaling::serialize(Archive &ar, unsigned version)
{
    if (version>i3clsimscalarfieldanisotropyabslenscaling_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimScalarFieldAnisotropyAbsLenScaling class.",version,i3clsimscalarfieldanisotropyabslenscaling_version_);

    ar & make_nvp("I3CLSimScalarField", base_object<I3CLSimScalarField>(*this));
    ar & make_nvp("anisotropyDirAzimuth", anisotropyDirAzimuth_);
    ar & make_nvp("magnitudeAlongDir", magnitudeAlongDir_);
    ar & make_nvp("magnitudePerpToDir", magnitudePerpToDir_);
}     


I3_SERIALIZABLE(I3CLSimScalarFieldAnisotropyAbsLenScaling);
