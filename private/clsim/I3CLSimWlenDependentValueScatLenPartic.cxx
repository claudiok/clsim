#include <icetray/serialization.h>
#include <icetray/I3Units.h>
#include <clsim/I3CLSimWlenDependentValueScatLenPartic.h>

#include <typeinfo>
#include <cmath>

#include "clsim/to_float_string.h"
using namespace I3CLSimHelper;

const double I3CLSimWlenDependentValueScatLenPartic::default_volumeConcentrationSmallParticles=0.0075f*I3Units::perMillion; 
const double I3CLSimWlenDependentValueScatLenPartic::default_volumeConcentrationLargeParticles=0.0075f*I3Units::perMillion; 

I3CLSimWlenDependentValueScatLenPartic::I3CLSimWlenDependentValueScatLenPartic(double volumeConcentrationSmallParticles,
                                                                               double volumeConcentrationLargeParticles)
:
volumeConcentrationSmallParticles_(volumeConcentrationSmallParticles/I3Units::perMillion), // in ppm
volumeConcentrationLargeParticles_(volumeConcentrationLargeParticles/I3Units::perMillion)  // in ppm
{ 
    
}

I3CLSimWlenDependentValueScatLenPartic::~I3CLSimWlenDependentValueScatLenPartic() 
{;}


double I3CLSimWlenDependentValueScatLenPartic::GetValue(double wlen) const
{
    const double refWlen = 550.f*I3Units::nanometer;
	const double x = refWlen/wlen;
	
    const double scatCoeff = 
    0.0017 * std::pow(x,4.3) + 
    1.34  * volumeConcentrationSmallParticles_ * std::pow(x,1.7) + 
    0.312 * volumeConcentrationLargeParticles_ * std::pow(x,0.3);
	
    return 1./scatCoeff;
}

std::string I3CLSimWlenDependentValueScatLenPartic::GetOpenCLFunction(const std::string &functionName) const
{
    std::string funcDef = 
    std::string("inline float ") + functionName + std::string("(float wavelength)\n");
    
    std::string funcBody = std::string() + 
    "{\n"
    "    const float volumeConcentrationSmallParticles=" + to_float_string(volumeConcentrationSmallParticles_) + ";   // in ppm\n"
    "    const float volumeConcentrationLargeParticles=" + to_float_string(volumeConcentrationLargeParticles_) + ";   // in ppm\n"
    "    const float refWlen = " + to_float_string(550.*I3Units::nanometer) + ";\n"
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
    
    return funcDef + funcBody;
}

bool I3CLSimWlenDependentValueScatLenPartic::CompareTo(const I3CLSimWlenDependentValue &other) const
{
    try
    {
        const I3CLSimWlenDependentValueScatLenPartic &other_ = dynamic_cast<const I3CLSimWlenDependentValueScatLenPartic &>(other);
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
void I3CLSimWlenDependentValueScatLenPartic::serialize(Archive &ar, unsigned version)
{
    if (version>i3clsimwlendependentvaluescatlenpartic_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimWlenDependentValueScatLenPartic class.",version,i3clsimwlendependentvaluescatlenpartic_version_);

    ar & make_nvp("I3CLSimWlenDependentValue", base_object<I3CLSimWlenDependentValue>(*this));
    ar & make_nvp("volumeConcentrationSmallParticles", volumeConcentrationSmallParticles_);
    ar & make_nvp("volumeConcentrationLargeParticles", volumeConcentrationLargeParticles_);
}     


I3_SERIALIZABLE(I3CLSimWlenDependentValueScatLenPartic);
