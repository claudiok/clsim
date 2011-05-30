#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include "clsim/I3CLSimSimpleGeometryFromI3Geometry.h"

#include <stdexcept>
#include <limits>

#include <boost/foreach.hpp>

// TODO: these defaults are IceCube-specific!
const int32_t I3CLSimSimpleGeometryFromI3Geometry::default_ignoreStringIDsSmallerThan = 1;
const int32_t I3CLSimSimpleGeometryFromI3Geometry::default_ignoreStringIDsLargerThan = std::numeric_limits<int32_t>::max();
const uint32_t I3CLSimSimpleGeometryFromI3Geometry::default_ignoreDomIDsSmallerThan = 1;
const uint32_t I3CLSimSimpleGeometryFromI3Geometry::default_ignoreDomIDsLargerThan = 60;

I3CLSimSimpleGeometryFromI3Geometry::
I3CLSimSimpleGeometryFromI3Geometry(double OMRadius,
                                    const I3GeometryConstPtr &geometry,
                                    const std::set<int> &ignoreStrings,
                                    const std::set<unsigned int> &ignoreDomIDs,
                                    int32_t ignoreStringIDsSmallerThan,
                                    int32_t ignoreStringIDsLargerThan,
                                    uint32_t ignoreDomIDsSmallerThan,
                                    uint32_t ignoreDomIDsLargerThan)
:
OMRadius_(OMRadius),
ignoreStrings_(ignoreStrings),
ignoreDomIDs_(ignoreDomIDs),
ignoreStringIDsSmallerThan_(ignoreStringIDsSmallerThan),
ignoreStringIDsLargerThan_(ignoreStringIDsLargerThan),
ignoreDomIDsSmallerThan_(ignoreDomIDsSmallerThan),
ignoreDomIDsLargerThan_(ignoreDomIDsLargerThan)
{
    if (!geometry) throw std::runtime_error("Received NULL geometry pointer!");
    
    log_debug("Ignoring StringNum<%" PRIi32 ", StringNum>%" PRIi32 ", OMNum<%" PRIu32 ", OMNum>%" PRIu32 ".",
              ignoreStringIDsSmallerThan, ignoreStringIDsLargerThan,
              ignoreDomIDsSmallerThan, ignoreDomIDsLargerThan);
    
    numOMs_=0;
    BOOST_FOREACH(const I3OMGeoMap::value_type &i, geometry->omgeo)
    {
        const OMKey &key = i.first;
        const I3OMGeo &geo = i.second;
        
        int32_t string=key.GetString();
        uint32_t dom=key.GetOM();

        if ((string < ignoreStringIDsSmallerThan_) ||
            (string > ignoreStringIDsLargerThan_) ||
            (dom < ignoreDomIDsSmallerThan_) ||
            (dom > ignoreDomIDsLargerThan_))
            continue;

        if (ignoreStrings_.count(string)!=0) continue;
        if (ignoreDomIDs_.count(dom)!=0) continue;
        
        stringIDs_.push_back(string);
        domIDs_.push_back(dom);
        posX_.push_back(geo.position.GetX());
        posY_.push_back(geo.position.GetY());
        posZ_.push_back(geo.position.GetZ());
        ++numOMs_;
    }
    
}

I3CLSimSimpleGeometryFromI3Geometry::
~I3CLSimSimpleGeometryFromI3Geometry()
{
    
}
