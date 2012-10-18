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
 * @file I3CLSimFunctionScatLenPartic.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <icetray/serialization.h>
#include <icetray/I3Units.h>
#include <clsim/function/I3CLSimFunctionScatLenPartic.h>

#include <typeinfo>
#include <cmath>

#include "clsim/I3CLSimHelperToFloatString.h"
using namespace I3CLSimHelper;

const double I3CLSimFunctionScatLenPartic::default_volumeConcentrationSmallParticles=0.0075f*I3Units::perMillion; 
const double I3CLSimFunctionScatLenPartic::default_volumeConcentrationLargeParticles=0.0075f*I3Units::perMillion; 

I3CLSimFunctionScatLenPartic::I3CLSimFunctionScatLenPartic(double volumeConcentrationSmallParticles,
                                                                               double volumeConcentrationLargeParticles)
:
volumeConcentrationSmallParticles_(volumeConcentrationSmallParticles/I3Units::perMillion), // in ppm
volumeConcentrationLargeParticles_(volumeConcentrationLargeParticles/I3Units::perMillion)  // in ppm
{ 
    
}

I3CLSimFunctionScatLenPartic::~I3CLSimFunctionScatLenPartic() 
{;}


double I3CLSimFunctionScatLenPartic::GetValue(double wlen) const
{
    const double refWlen = 550.f*I3Units::nanometer;
    const double x = refWlen/wlen;

    const double scatCoeff = 
    0.0017 * std::pow(x,4.3) + 
    1.34  * volumeConcentrationSmallParticles_ * std::pow(x,1.7) + 
    0.312 * volumeConcentrationLargeParticles_ * std::pow(x,0.3);

    return 1./scatCoeff;
}

std::string I3CLSimFunctionScatLenPartic::GetOpenCLFunction(const std::string &functionName) const
{
    std::string funcDef = 
    std::string("inline float ") + functionName + std::string("(float wavelength)\n");

    std::string funcBody = std::string() + 
    "{\n"
    "    const float volumeConcentrationSmallParticles=" + ToFloatString(volumeConcentrationSmallParticles_) + ";   // in ppm\n"
    "    const float volumeConcentrationLargeParticles=" + ToFloatString(volumeConcentrationLargeParticles_) + ";   // in ppm\n"
    "    const float refWlen = " + ToFloatString(550.*I3Units::nanometer) + ";\n"
    "    \n"
    "#ifdef USE_NATIVE_MATH\n"
    "    const float x = native_divide(refWlen,wavelength);\n"
    "    const float scatCoeff = 0.0017f * native_powr(x,4.3f) + \n"
    "                            1.34f   * volumeConcentrationSmallParticles * native_powr(x,1.7f) + \n"
    "                            0.312f  * volumeConcentrationLargeParticles * native_powr(x,0.3f);\n"
    "    return native_recip(scatCoeff);\n"
    "#else\n"
    "    const float x = refWlen/wavelength;\n"
    "    const float scatCoeff = 0.0017f * powr(x,4.3f) + \n"
    "                            1.34f   * volumeConcentrationSmallParticles * powr(x,1.7f) + \n"
    "                            0.312f  * volumeConcentrationLargeParticles * powr(x,0.3f);\n"
    "    return 1.f/scatCoeff;\n"
    "#endif\n"
    "}\n"
    ;

    return funcDef + ";\n\n" + funcDef + funcBody;
}

bool I3CLSimFunctionScatLenPartic::CompareTo(const I3CLSimFunction &other) const
{
    try
    {
        const I3CLSimFunctionScatLenPartic &other_ = dynamic_cast<const I3CLSimFunctionScatLenPartic &>(other);
        return ((other_.volumeConcentrationSmallParticles_ == volumeConcentrationSmallParticles_) &&
                (other_.volumeConcentrationLargeParticles_ == volumeConcentrationLargeParticles_));
    }
    catch (const std::bad_cast& e)
    {
        // not of the same type, treat it as non-equal
        return false;
    }

}




template <class Archive>
void I3CLSimFunctionScatLenPartic::serialize(Archive &ar, unsigned version)
{
    if (version>i3clsimfunctionscatlenpartic_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimFunctionScatLenPartic class.",version,i3clsimfunctionscatlenpartic_version_);

    ar & make_nvp("I3CLSimFunction", base_object<I3CLSimFunction>(*this));
    ar & make_nvp("volumeConcentrationSmallParticles", volumeConcentrationSmallParticles_);
    ar & make_nvp("volumeConcentrationLargeParticles", volumeConcentrationLargeParticles_);
}     


I3_SERIALIZABLE(I3CLSimFunctionScatLenPartic);
