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
 * @file I3ShadowedPhotonRemoverModule.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif 
#include <inttypes.h>

#include "clsim/shadow/I3ShadowedPhotonRemoverModule.h"

#include <boost/foreach.hpp>

#include "clsim/I3Photon.h"

//#include "dataclasses/geometry/I3Geometry.h"

#include "dataclasses/I3Constants.h"

// The module
I3_MODULE(I3ShadowedPhotonRemoverModule);

/**
 * This module is not ready for use yet.
 */
I3ShadowedPhotonRemoverModule::I3ShadowedPhotonRemoverModule(const I3Context& context) 
: I3ConditionalModule(context)
{
    inputPhotonSeriesMapName_="PropagatedPhotons";
    AddParameter("InputPhotonSeriesMapName",
                 "Name of the input I3PhotonSeriesMap frame object.",
                 inputPhotonSeriesMapName_);

    outputPhotonSeriesMapName_="PropagatedPhotonsWithShadow";
    AddParameter("OutputPhotonSeriesMapName",
                 "Name of the output I3PhotonSeriesMap frame object.",
                 outputPhotonSeriesMapName_);

    // add an outbox
    AddOutBox("OutBox");

}

I3ShadowedPhotonRemoverModule::~I3ShadowedPhotonRemoverModule()
{
    log_trace("%s", __PRETTY_FUNCTION__);
}


void I3ShadowedPhotonRemoverModule::Configure()
{
    log_trace("%s", __PRETTY_FUNCTION__);

    GetParameter("InputPhotonSeriesMapName", inputPhotonSeriesMapName_);
    GetParameter("OutputPhotonSeriesMapName", outputPhotonSeriesMapName_);

    // set up the worker class
    shadowedPhotonRemover_ = I3ShadowedPhotonRemoverPtr(new I3ShadowedPhotonRemover());

}


void I3ShadowedPhotonRemoverModule::DAQ(I3FramePtr frame)
{
    log_trace("%s", __PRETTY_FUNCTION__);
    
    //const I3Geometry& geometry = frame->Get<I3Geometry>();

    I3PhotonSeriesMapConstPtr inputPhotonSeriesMap = frame->Get<I3PhotonSeriesMapConstPtr>(inputPhotonSeriesMapName_);
    if (!inputPhotonSeriesMap) log_fatal("Frame does not contain an I3PhotonSeriesMap named \"%s\".",
                                         inputPhotonSeriesMapName_.c_str());
    
    // allocate the output hitSeriesMap
    I3PhotonSeriesMapPtr outputPhotonSeriesMap(new I3PhotonSeriesMap());
    
    BOOST_FOREACH(const I3PhotonSeriesMap::value_type &it, *inputPhotonSeriesMap)
    {
        const ModuleKey &key = it.first;
        const I3PhotonSeries &photons = it.second;

        //// Find the current OM in the geometry map
        //I3OMGeoMap::const_iterator geo_it = geometry.omgeo.find(key);
        //if (geo_it == geometry.omgeo.end())
        //        log_fatal("OM (%i/%u) not found in the current geometry map!", key.GetString(), key.GetOM());
        //const I3OMGeo &om = geo_it->second;

        
        // a pointer to the output vector. The vector will be allocated 
        // by the map, this is merely a pointer to it in case we have multiple
        // photons per OM.
        I3PhotonSeries *out_photons = NULL;

        BOOST_FOREACH(const I3Photon &photon, photons)
        {
            
            I3Photon out_photon = photon;
            const bool isShadowed = 
            shadowedPhotonRemover_->IsPhotonShadowed(out_photon);

            if (isShadowed) continue;

            // allocate the output vector if not already done
            if (!out_photons) out_photons = &(outputPhotonSeriesMap->insert(std::make_pair(key, I3PhotonSeries())).first->second);

            // add a new copy of the input photon to the output list
            out_photons->push_back(out_photon);
        }
        
        
    }
    
    // store the output I3MCHitSeriesMap
    frame->Put(outputPhotonSeriesMapName_, outputPhotonSeriesMap);
    
    // that's it!
    PushFrame(frame);
}
