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
    if (isnan(fromWlen_)) log_fatal("The \"fromWlen\" argument must not be NaN!");
    if (isnan(toWlen_)) log_fatal("The \"toWlen\" argument must not be NaN!");
              
    if (fromWlen_ > toWlen_) 
        log_fatal("The \"fromWlen\" argument must not be greater than \"toWlen\".");
}

I3CLSimRandomValueWlenCherenkovNoDispersion::~I3CLSimRandomValueWlenCherenkovNoDispersion() 
{ 
}

I3CLSimRandomValueWlenCherenkovNoDispersion::I3CLSimRandomValueWlenCherenkovNoDispersion() {;}



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
