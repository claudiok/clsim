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
 * @file I3CLSimSimpleGeometryFromI3Geometry.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>

#include <icetray/I3Units.h>
#include "clsim/I3CLSimSimpleGeometryFromI3Geometry.h"

#include "dataclasses/geometry/I3Geometry.h"
#include "dataclasses/geometry/I3OMGeo.h"

#include "dataclasses/geometry/I3ModuleGeo.h"

#include <stdexcept>
#include <limits>

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

const std::set<int> I3CLSimSimpleGeometryFromI3Geometry::default_ignoreStrings;
const std::set<unsigned int> I3CLSimSimpleGeometryFromI3Geometry::default_ignoreDomIDs;
const std::set<std::string> I3CLSimSimpleGeometryFromI3Geometry::default_ignoreSubdetectors;
const int32_t I3CLSimSimpleGeometryFromI3Geometry::default_ignoreStringIDsSmallerThan = 1;
const int32_t I3CLSimSimpleGeometryFromI3Geometry::default_ignoreStringIDsLargerThan = std::numeric_limits<int32_t>::max();
const uint32_t I3CLSimSimpleGeometryFromI3Geometry::default_ignoreDomIDsSmallerThan = 1;
const uint32_t I3CLSimSimpleGeometryFromI3Geometry::default_ignoreDomIDsLargerThan = 60;
const bool I3CLSimSimpleGeometryFromI3Geometry::default_splitIntoPartsAccordingToPosition=false;
const bool I3CLSimSimpleGeometryFromI3Geometry::default_useHardcodedDeepCoreSubdetector=false;


I3CLSimSimpleGeometryFromI3Geometry::
I3CLSimSimpleGeometryFromI3Geometry(double OMRadius,
                                    double oversizeFactor,
                                    const I3FramePtr &frame,
                                    const std::set<int> &ignoreStrings,
                                    const std::set<unsigned int> &ignoreDomIDs,
                                    const std::set<std::string> &ignoreSubdetectors,
                                    int32_t ignoreStringIDsSmallerThan,
                                    int32_t ignoreStringIDsLargerThan,
                                    uint32_t ignoreDomIDsSmallerThan,
                                    uint32_t ignoreDomIDsLargerThan,
                                    bool splitIntoPartsAccordingToPosition,
                                    bool useHardcodedDeepCoreSubdetector)
:
OMRadius_(OMRadius),
oversizeFactor_(oversizeFactor),
splitIntoPartsAccordingToPosition_(splitIntoPartsAccordingToPosition),
ignoreStrings_(ignoreStrings),
ignoreDomIDs_(ignoreDomIDs),
ignoreSubdetectors_(ignoreSubdetectors),
ignoreStringIDsSmallerThan_(ignoreStringIDsSmallerThan),
ignoreStringIDsLargerThan_(ignoreStringIDsLargerThan),
ignoreDomIDsSmallerThan_(ignoreDomIDsSmallerThan),
ignoreDomIDsLargerThan_(ignoreDomIDsLargerThan),
useHardcodedDeepCoreSubdetector_(useHardcodedDeepCoreSubdetector)
{
    if (!frame) throw std::runtime_error("Received NULL frame pointer!");
    
    log_debug("Ignoring StringNum<%" PRIi32 ", StringNum>%" PRIi32 ", OMNum<%" PRIu32 ", OMNum>%" PRIu32 ".",
              ignoreStringIDsSmallerThan, ignoreStringIDsLargerThan,
              ignoreDomIDsSmallerThan, ignoreDomIDsLargerThan);
    
    numOMs_=0;
    
    I3ModuleGeoMapConstPtr moduleGeoMap = frame->Get<I3ModuleGeoMapConstPtr>("I3ModuleGeoMap");
    
    if (!moduleGeoMap) {
        if (frame->Has("I3Geometry")) {
            log_fatal("No I3ModuleGeoMap found in frame. There *is* an I3Geometry object. Please run the \"I3GeometryDecomposer\" module before this.");
        } else {
            log_fatal("No I3ModuleGeoMap found in frame. There does not seem to be a geometry.");
        }
    }

    I3MapModuleKeyStringConstPtr subdetectors = frame->Get<I3MapModuleKeyStringConstPtr>("Subdetectors");
    if (!subdetectors) log_error("No subdetector configuration in frame. Missing a \"Subdetectors\" object. Assuming all modules are on the same detector.");
    
    BOOST_FOREACH(const I3ModuleGeoMap::value_type &i, *moduleGeoMap)
    {
        const ModuleKey &key = i.first;
        const I3ModuleGeo &geo = i.second;
        
        std::string subdetectorName = "Unknown"; // use this if there is no Subdetectors object
        if (subdetectors) {
            I3MapModuleKeyString::const_iterator subdetector_it =
            subdetectors->find(key);
            if (subdetector_it == subdetectors->end()) {
                log_fatal("ModuleKey(%i/%u) not found in \"Subdetectors\".",
                          key.GetString(), key.GetOM());
            }
            subdetectorName = subdetector_it->second;
        }
        
        int32_t string=key.GetString();
        uint32_t dom=key.GetOM();

        if (useHardcodedDeepCoreSubdetector_) {
            // special hack for DeepCore
            if ((subdetectorName=="IceCube") || (subdetectorName=="DeepCore"))
            {
                if ((string>=79) && (string<=86)) // these are the DeepCore strings
                {
                    if (geo.GetPos().GetZ()>-30.*I3Units::m) // z=30m is about halfway between the upper and lower parts of DeepCore
                        subdetectorName="DeepCoreUpper";
                    else
                        subdetectorName="DeepCoreLower";
                }
                else if (string > 86)
                {
                    subdetectorName="PINGU";
                }
            }
        }
        
        if ((string < ignoreStringIDsSmallerThan_) ||
            (string > ignoreStringIDsLargerThan_) ||
            (dom < ignoreDomIDsSmallerThan_) ||
            (dom > ignoreDomIDsLargerThan_))
            continue;

        if (ignoreStrings_.count(string)!=0) continue;
        if (ignoreDomIDs_.count(dom)!=0) continue;
        if (ignoreSubdetectors_.count(subdetectorName)!=0) continue;

        
        // sanity check
        if (std::abs(geo.GetRadius()-OMRadius_) > 0.001*I3Units::mm)
            log_fatal("This version of clsim does only support DOMs with one single size. Configured size=%fmm, size in geometry=%fmm",
                      OMRadius_/I3Units::mm, geo.GetRadius()/I3Units::mm);
        
        stringIDs_.push_back(string);
        domIDs_.push_back(dom);
        posX_.push_back(geo.GetPos().GetX());
        posY_.push_back(geo.GetPos().GetY());
        posZ_.push_back(geo.GetPos().GetZ());
        subdetectors_.push_back(subdetectorName);

        ++numOMs_;
    }
    
}

I3CLSimSimpleGeometryFromI3Geometry::
~I3CLSimSimpleGeometryFromI3Geometry()
{
    
}
