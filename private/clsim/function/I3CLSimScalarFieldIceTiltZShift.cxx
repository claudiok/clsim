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
 * @file I3CLSimScalarFieldIceTiltZShift.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include <icetray/serialization.h>
#include <clsim/function/I3CLSimScalarFieldIceTiltZShift.h>
#include <icetray/I3Units.h>

#include <typeinfo>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include "clsim/I3CLSimHelperToFloatString.h"
using namespace I3CLSimHelper;

const double I3CLSimScalarFieldIceTiltZShift::default_directionOfTiltAzimuth = 225.*I3Units::deg;

I3CLSimScalarFieldIceTiltZShift::
I3CLSimScalarFieldIceTiltZShift(
    const std::vector<double> &distancesFromOriginAlongTilt,
    const std::vector<double> &zCoordinates,
    const I3Matrix &zCorrections,
    double directionOfTiltAzimuth)
:
distancesFromOriginAlongTilt_(distancesFromOriginAlongTilt),
zCoordinates_(zCoordinates),
zCorrections_(zCorrections),
directionOfTiltAzimuth_(directionOfTiltAzimuth)
{
    if (distancesFromOriginAlongTilt_.size() != zCorrections_.size1())
        throw std::range_error("zCorrection array size in first dimension differs from size of distancesFromOriginAlongTilt!");
    if (zCoordinates_.size() != zCorrections_.size2())
        throw std::range_error("zCorrection array size in second dimension differs from size of zCoordinates!");

    if (distancesFromOriginAlongTilt_.size()<2)
        throw std::runtime_error("distancesFromOriginAlongTilt (dimension 1) needs at least 2 entries.");
    if (zCoordinates_.size()<2)
        throw std::runtime_error("zCoordinates (dimension 2) needs at least 2 entries.");

    // ensure that the z-coordinates are equally spaced
    double meanSpacing=0.;
    for (std::size_t i=0;i<zCoordinates_.size()-1;++i)
    {
        const double thisSpacing = zCoordinates_[i+1]-zCoordinates_[i];

        if (thisSpacing <= 0.)
            throw std::runtime_error("zCoordinates (dimension 2) are not in ascending order.");

        meanSpacing += thisSpacing;
    }
    meanSpacing /= static_cast<double>(zCoordinates_.size()-1);

    for (std::size_t i=0;i<zCoordinates_.size()-1;++i)
    {
        const double thisSpacing = zCoordinates_[i+1]-zCoordinates_[i];
        const double spacingDiff = std::abs(thisSpacing-meanSpacing);

        if (spacingDiff > 1e-5)
            throw std::runtime_error("zCoordinates (dimension 2) are not in equally spaced.");
    }

    firstZCoordinate_ = zCoordinates_[0];
    zCoordinateSpacing_ = meanSpacing;
    log_debug("%zu z-coordinates starting at %fm, spacing is %fm",
        zCoordinates_.size(),
        firstZCoordinate_/I3Units::m,
        zCoordinateSpacing_/I3Units::m);

    // check the horizontal distance
    for (std::size_t i=0;i<distancesFromOriginAlongTilt_.size()-1;++i)
    {
        const double thisSpacing = distancesFromOriginAlongTilt_[i+1]-distancesFromOriginAlongTilt_[i];

        if (thisSpacing <= 0.)
            throw std::runtime_error("distancesFromOriginAlongTilt (dimension 1) is not in ascending order.");
    }

}

I3CLSimScalarFieldIceTiltZShift::I3CLSimScalarFieldIceTiltZShift() {;}

I3CLSimScalarFieldIceTiltZShift::~I3CLSimScalarFieldIceTiltZShift() 
{;}

bool I3CLSimScalarFieldIceTiltZShift::HasNativeImplementation() const 
{
    return true;
}

double I3CLSimScalarFieldIceTiltZShift::GetValue(double x, double y, double z) const
{
    // this is more or less taken from PPC's zshift() (in pro.cu)
    const double lnx = std::cos(directionOfTiltAzimuth_);
    const double lny = std::sin(directionOfTiltAzimuth_);

    const double z_rescaled = (z-firstZCoordinate_)/zCoordinateSpacing_;
    const std::size_t k = std::min(std::max(std::floor(z_rescaled), 0.), static_cast<double>(zCoordinates_.size()-2));

    const double fraction_z_above = z_rescaled-static_cast<double>(k);
    const double fraction_z_below = static_cast<double>(k+1)-z_rescaled;

    const double nr = lnx*x+lny*y;

    for(int j=1; j<static_cast<int>(distancesFromOriginAlongTilt_.size()); j++) 
    {
        if((nr<distancesFromOriginAlongTilt_[j]) || (j==static_cast<int>(distancesFromOriginAlongTilt_.size())-1))
        {
            const double thisDistanceBinWidth = (distancesFromOriginAlongTilt_[j] - distancesFromOriginAlongTilt_[j-1]);

            const double frac_at_lower = (distancesFromOriginAlongTilt_[j] - nr  )/thisDistanceBinWidth;
            const double frac_at_upper = (nr - distancesFromOriginAlongTilt_[j-1])/thisDistanceBinWidth;

            const double val_at_lower = (zCorrections_(j-1, k+1)*fraction_z_above + zCorrections_(j-1, k)*fraction_z_below);
            const double val_at_upper = (zCorrections_(j,   k+1)*fraction_z_above + zCorrections_(j,   k)*fraction_z_below);

            return (val_at_upper * frac_at_upper + val_at_lower * frac_at_lower);
        }
    }
    return 0.;
}

std::string I3CLSimScalarFieldIceTiltZShift::GetOpenCLFunction(const std::string &functionName) const
{
    // some names
    const std::string dataName = functionName + "_data";

    std::string dataDef;
    
    dataDef += std::string("\n") + 
    "#define " + dataName + "_numDistances  " + boost::lexical_cast<std::string>(zCorrections_.size1()) + "\n" +
    "#define " + dataName + "_numZCoords    " + boost::lexical_cast<std::string>(zCorrections_.size2()) + "\n" +
    "#define " + dataName + "_firstZCoord   " + ToFloatString(firstZCoordinate_) + "\n" +
    "#define " + dataName + "_zCoordSpacing " + ToFloatString(zCoordinateSpacing_) + "\n\n";

    dataDef +=
    "__constant float " + dataName + "_distancesFromOriginAlongTilt[" + dataName + "_numDistances] = {\n";
    for (std::size_t i=0;i<distancesFromOriginAlongTilt_.size();++i)
    {
        dataDef += ToFloatString(distancesFromOriginAlongTilt_[i]) + ", \n";
    }
    dataDef += "};\n\n";

    dataDef +=
    "__constant float " + dataName + "_zCorrections[" + dataName + "_numDistances*" + dataName + "_numZCoords] = {\n";
    for (std::size_t i=0;i<zCorrections_.size1();++i)
    {
        for (std::size_t j=0;j<zCorrections_.size2();++j)
        {
            dataDef += ToFloatString(zCorrections_(i,j)) + ", \n";

        }
        dataDef += " // distances[" + boost::lexical_cast<std::string>(i) + "]\n";
    }
    dataDef += "};\n\n";


    std::string funcDef = 
        std::string("inline float ") + functionName + std::string("(float4 vec)");
    
    const double lnx = std::cos(directionOfTiltAzimuth_);
    const double lny = std::sin(directionOfTiltAzimuth_);

    std::string funcBody = std::string() + 
    "{\n"
    "    const float z_rescaled = (vec.z-" + dataName + "_firstZCoord)/" + dataName + "_zCoordSpacing;\n"
    "    const int k = min(max(convert_int_rtn(z_rescaled), 0), " + dataName + "_numZCoords-2);\n"
    "    \n"
    "    const float fraction_z_above = z_rescaled-convert_float(k);\n"
    "    const float fraction_z_below = 1.-fraction_z_above;\n"
    "    \n"
    "    const float nr = " + ToFloatString(lnx) + "*vec.x + " + ToFloatString(lny) + "*vec.y;\n"
    "    \n"
    "    for(int j=1; j<" + dataName + "_numDistances; j++)\n"
    "    {\n"
    "        const float thisDist = " + dataName + "_distancesFromOriginAlongTilt[j];\n"
    "        if((nr<thisDist) || (j==" + dataName + "_numDistances-1))\n"
    "        {\n"
    "            const float previousDist = " + dataName + "_distancesFromOriginAlongTilt[j-1];\n"
    "            const float thisDistanceBinWidth = thisDist - previousDist;\n"
    "            \n"
    "            const float frac_at_lower = (thisDist - nr    )/thisDistanceBinWidth;\n"
    "            const float frac_at_upper = 1.-frac_at_lower;\n"
    "            \n"
    "            const float val_at_lower = (" + dataName + "_zCorrections[(j-1)*" + dataName + "_numZCoords + k+1]*fraction_z_above + " + dataName + "_zCorrections[(j-1)*" + dataName + "_numZCoords + k]*fraction_z_below);\n"
    "            const float val_at_upper = (" + dataName + "_zCorrections[j    *" + dataName + "_numZCoords + k+1]*fraction_z_above + " + dataName + "_zCorrections[j    *" + dataName + "_numZCoords + k]*fraction_z_below);\n"
    "            \n"
    "            return (val_at_upper * frac_at_upper + val_at_lower * frac_at_lower);\n"
    "        }\n"
    "    }\n"
    "}\n"
    ;

    return dataDef + "\n\n" + funcDef + ";\n\n" + funcDef + "\n" + funcBody;
}

bool I3CLSimScalarFieldIceTiltZShift::CompareTo(const I3CLSimScalarField &other) const
{
    try
    {
        const I3CLSimScalarFieldIceTiltZShift &other_ = dynamic_cast<const I3CLSimScalarFieldIceTiltZShift &>(other);

        if (other_.distancesFromOriginAlongTilt_.size() != distancesFromOriginAlongTilt_.size())
            return false;
        for (std::size_t i=0;i<distancesFromOriginAlongTilt_.size();++i)
        {
            if (other_.distancesFromOriginAlongTilt_[i] != distancesFromOriginAlongTilt_[i])
                return false;
        }

        if (other_.zCoordinates_.size() != zCoordinates_.size())
            return false;
        for (std::size_t i=0;i<zCoordinates_.size();++i)
        {
            if (other_.zCoordinates_[i] != zCoordinates_[i])
                return false;
        }


        if ((other_.zCorrections_.size1() != zCorrections_.size1()) ||
            (other_.zCorrections_.size2() != zCorrections_.size2()))
            return false;
        for (std::size_t i=0;i<zCorrections_.size1();++i)
            for (std::size_t j=0;j<zCorrections_.size2();++j)
            {
                if (other_.zCorrections_(i,j) != zCorrections_(i,j))
                    return false;
            }

        if (other_.directionOfTiltAzimuth_ != directionOfTiltAzimuth_)
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
void I3CLSimScalarFieldIceTiltZShift::serialize(Archive &ar, unsigned version)
{
    if (version>i3clsimscalarfieldicetiltzshift_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimScalarFieldIceTiltZShift class.",version,i3clsimscalarfieldicetiltzshift_version_);

    ar & make_nvp("I3CLSimScalarField", base_object<I3CLSimScalarField>(*this));
    ar & make_nvp("distancesFromOriginAlongTilt", distancesFromOriginAlongTilt_);
    ar & make_nvp("zCoordinates", zCoordinates_);
    ar & make_nvp("zCorrections", zCorrections_);
    ar & make_nvp("directionOfTiltAzimuth", directionOfTiltAzimuth_);

    ar & make_nvp("firstZCoordinate", firstZCoordinate_);
    ar & make_nvp("zCoordinateSpacing", zCoordinateSpacing_);

}     


I3_SERIALIZABLE(I3CLSimScalarFieldIceTiltZShift);
