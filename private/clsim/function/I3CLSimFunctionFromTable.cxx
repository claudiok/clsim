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
 * @file I3CLSimFunctionFromTable.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#define __CL_ENABLE_EXCEPTIONS
#include "clsim/cl.hpp"

#include <icetray/serialization.h>
#include <clsim/function/I3CLSimFunctionFromTable.h>

#include <typeinfo>
#include <cmath>
#include <math.h>
#include <stdexcept>

#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>

#include "clsim/I3CLSimHelperToFloatString.h"
using namespace I3CLSimHelper;

// when storing data in 16bits, this class can optionally
// use IEEE half precision. This is not used by default,
// but can be activated using this ifdef:
 //#define USE_OPENCL_HALF_PRECISION

#ifdef USE_OPENCL_HALF_PRECISION
#include "opencl/ieeehalfprecision.h"
#endif

const bool I3CLSimFunctionFromTable::default_storeDataAsHalfPrecision=false;

I3CLSimFunctionFromTable::
I3CLSimFunctionFromTable(const std::vector<double> &wlens,
                         const std::vector<double> &values,
                         bool storeDataAsHalfPrecision)
:
wlenStep_(NAN),
wlens_(wlens),
values_(values),
equalSpacingMode_(false),
storeDataAsHalfPrecision_(storeDataAsHalfPrecision)
{
    if (wlens_.size() < 2) throw std::range_error("wlens must contain at least 2 elements!");
    if (wlens_.size() != values_.size()) throw std::range_error("wlens and values must have the same size!");

    startWlen_ = wlens_[0];
}

I3CLSimFunctionFromTable::
I3CLSimFunctionFromTable(double startWlen,
                         double wlenStep,
                         const std::vector<double> &values,
                         bool storeDataAsHalfPrecision)
:
startWlen_(startWlen),
wlenStep_(wlenStep),
values_(values),
equalSpacingMode_(true),
storeDataAsHalfPrecision_(storeDataAsHalfPrecision)
{
    if (wlenStep_ <= 0.) throw std::range_error("wlenStep must not be <= 0!");
    if (values_.size() < 2) throw std::length_error("values must contain at least 2 elements!");

    // fill the wavelength array, even though we don't really need it in this mode
    wlens_.clear();
    for (std::size_t i=0;i<values_.size();++i) {
        wlens_.push_back(startWlen_ + static_cast<double>(i)*wlenStep_);
    }
}

I3CLSimFunctionFromTable::I3CLSimFunctionFromTable() {;}

I3CLSimFunctionFromTable::~I3CLSimFunctionFromTable() 
{;}

namespace {
    inline double mix(double min, double max, double t)
    {
        return min+(max-min)*t;
    }
}

double I3CLSimFunctionFromTable::GetValue(double wlen) const
{
    if (equalSpacingMode_) {
        double fbin;
        double fraction = modf((wlen-startWlen_)/wlenStep_, &fbin);
        int ibin=static_cast<int>(fbin);

        if ((ibin<0) || ((ibin==0) && (fraction<0)))  {
            ibin=0;
            fraction=0.;
        } else if (static_cast<std::size_t>(ibin)>=values_.size()-1) {
            ibin=values_.size()-2;
            fraction=1.;
        }
        const unsigned int bin = ibin;

        return mix(values_[bin], 
                   values_[bin+1],
                   fraction);
    } else {
        if (wlen <= wlens_[0]) {
            return values_[0];
        }

        std::size_t bin;
        double fraction;
        for (std::size_t i=1;i<wlens_.size();++i) {
            if (wlen <= wlens_[i]) {
                // got it!
                bin=i-1;
                fraction = (wlen-wlens_[i-1])/(wlens_[i]-wlens_[i-1]);

                return mix(values_[bin], 
                           values_[bin+1],
                           fraction);
            }
        }

        // nothing in range
        return values_[wlens_.size()-1];
    }
}

double I3CLSimFunctionFromTable::GetMinWlen() const
{
    if (equalSpacingMode_) {
        return startWlen_;
    } else {
        return wlens_[0];
    }
}

double I3CLSimFunctionFromTable::GetMaxWlen() const
{
    if (equalSpacingMode_) {
        return startWlen_ + wlenStep_ * static_cast<double>(values_.size()-1);
    } else {
        return wlens_[wlens_.size()-1];
    }
}

std::string I3CLSimFunctionFromTable::GetOpenCLFunction(const std::string &functionName) const
{
    if (!equalSpacingMode_)
        log_fatal("GetOpenCLFunction() has not yet been implemented for non-equal bin spacings.");

    // some names
    const std::string dataName = functionName + "_data";
    const std::string interpHelperName = functionName + "_getInterpolationBinAndFraction";

#ifndef USE_OPENCL_HALF_PRECISION
    double smallestEntry=0., largestEntry=0.;
#endif

    std::string dataDef;
    if (storeDataAsHalfPrecision_)
    {
#ifdef USE_OPENCL_HALF_PRECISION
        // define the actual data
        dataDef = dataDef +
        "// this is IEEE half precision values stored as 16bit integers:\n" +
        "__constant unsigned short " + dataName + "[" + boost::lexical_cast<std::string>(values_.size()) + "] = {\n";
        BOOST_FOREACH(const double &val, values_)
        {
            // convert the double value to IEEE half precision and store it as an unsigned integer (16bit)
            uint16_t dummy=0;
            doubles2halfp(&dummy, &val, 1);
            dataDef += boost::lexical_cast<std::string>(dummy) + ", ";
        }
        dataDef += "\n";
        dataDef += "};\n\n";
#else
        // find the smallest and largest entries.
        bool firstIteration=true;
        BOOST_FOREACH(const double &val, values_)
        {
            if (firstIteration) {
                firstIteration=false;
                smallestEntry=val;
                largestEntry=val;
                continue;
            }

            if (val < smallestEntry) smallestEntry=val;
            if (val > largestEntry) largestEntry=val;
        }

        // define the actual data
        dataDef = dataDef +
        "#define " + dataName + "_SMALLEST_ENTRY " + ToFloatString(smallestEntry) + " \n" +
        "#define " + dataName + "_LARGEST_ENTRY " + ToFloatString(largestEntry) + " \n" +
        "__constant unsigned short " + dataName + "[" + boost::lexical_cast<std::string>(values_.size()) + "] = {\n";
        BOOST_FOREACH(const double &val, values_)
        {
            uint16_t dummy = static_cast<uint16_t>(65535.*(val-smallestEntry)/(largestEntry-smallestEntry));
            dataDef += boost::lexical_cast<std::string>(dummy) + ", ";
        }
        dataDef += "\n";
        dataDef += "};\n\n";
#endif
    }
    else
    {
        // define the actual data
        dataDef = dataDef +
        "__constant float " + dataName + "[" + boost::lexical_cast<std::string>(values_.size()) + "] = {\n";
        BOOST_FOREACH(const double &val, values_)
        {
            dataDef += ToFloatString(val) + ", ";
        }
        dataDef += "\n";
        dataDef += "};\n\n";
    }

    std::string interpHelperDef =
    std::string("inline void ") + interpHelperName + "(float wavelength, int *bin, float *fraction)";

    std::string interpHelperBody =
    "{\n"
    "    float fbin;\n"
    "    *fraction = modf((wavelength - " + ToFloatString(startWlen_) + ")/" + ToFloatString(wlenStep_) + ", &fbin);\n"
    "    \n"
    "    int ibin=(int)fbin;\n"
    "    \n"
    "    if ((ibin<0) || ((ibin==0) && (*fraction<0))) {\n"
    "        ibin=0;\n"
    "        *fraction=0.f;\n"
    "    } else if (ibin>=" + boost::lexical_cast<std::string>(values_.size()) + "-1) {\n"
    "        ibin=" + boost::lexical_cast<std::string>(values_.size()) + "-2;\n"
    "        *fraction=1.f;\n"
    "    }\n"
    "    \n"
    "    *bin = ibin;\n"
    "}\n\n"
    ;

    std::string funcDef = 
    std::string("inline float ") + functionName + std::string("(float wavelength)\n");

    std::string funcBody;
    if (storeDataAsHalfPrecision_)
    {
#ifdef USE_OPENCL_HALF_PRECISION
        funcBody = funcBody + 
        "{\n"
        "    int bin; float fraction;\n"
        "    " + interpHelperName + "(wavelength, &bin, &fraction);\n"
        "    \n"
        "    return mix(vload_half(bin,   (const __constant half *)" + dataName + "),\n"
        "               vload_half(bin+1, (const __constant half *)" + dataName + "),\n"
        "               fraction);\n"
        "}\n"
        ;
#else
        funcBody = funcBody + 
        "{\n"
        "    int bin; float fraction;\n"
        "    " + interpHelperName + "(wavelength, &bin, &fraction);\n"
        "    \n"
            "    return mix(convert_float(" + dataName + "[bin])  *((" + dataName + "_LARGEST_ENTRY-" + dataName + "_SMALLEST_ENTRY)/65535.f) + " + dataName + "_SMALLEST_ENTRY,\n"
            "               convert_float(" + dataName + "[bin+1])*((" + dataName + "_LARGEST_ENTRY-" + dataName + "_SMALLEST_ENTRY)/65535.f) + " + dataName + "_SMALLEST_ENTRY,\n"
            "               fraction);\n"
        "}\n"
        ;
#endif
    }
    else
    {
        funcBody = funcBody + 
        "{\n"
        "    int bin; float fraction;\n"
        "    " + interpHelperName + "(wavelength, &bin, &fraction);\n"
        "    \n"
        "    return mix(" + dataName + "[bin],\n"
        "               " + dataName + "[bin+1],\n"
        "               fraction);\n"
        "}\n"
        ;
    }

    return funcDef + ";\n" + interpHelperDef + ";\n" + dataDef + interpHelperDef + interpHelperBody + funcDef + funcBody;
}

bool I3CLSimFunctionFromTable::CompareTo(const I3CLSimFunction &other) const
{
    try
    {
        const I3CLSimFunctionFromTable &other_ = dynamic_cast<const I3CLSimFunctionFromTable &>(other);
        if ((other_.equalSpacingMode_ != equalSpacingMode_) ||
            (other_.values_.size() != values_.size()))
            return false;

        for (std::size_t i=0;i<values_.size();++i)
        {
            if (other_.values_[i] != values_[i]) return false;
        }

        if (equalSpacingMode_) {
            if (other_.startWlen_ != startWlen_) return false;
            if (other_.wlenStep_ != wlenStep_) return false;
        } else {
            for (std::size_t i=0;i<wlens_.size();++i)
            {
                if (other_.wlens_[i] != wlens_[i]) return false;
            }
        }

        return true;
    }
    catch (const std::bad_cast& e)
    {
        // not of the same type, treat it as non-equal
        return false;
    }

}




template <class Archive>
void I3CLSimFunctionFromTable::serialize(Archive &ar, unsigned version)
{
    if (version>i3clsimfunctionfromtable_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimFunctionFromTable class.",version,i3clsimfunctionfromtable_version_);

    ar & make_nvp("I3CLSimFunction", base_object<I3CLSimFunction>(*this));
    ar & make_nvp("startWlen", startWlen_);
    ar & make_nvp("wlenStep", wlenStep_);
    ar & make_nvp("wlens", wlens_);
    ar & make_nvp("values", values_);
    ar & make_nvp("equalSpacingMode", equalSpacingMode_);
}


I3_SERIALIZABLE(I3CLSimFunctionFromTable);
