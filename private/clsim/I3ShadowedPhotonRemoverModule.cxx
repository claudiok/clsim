/**
 * copyright (C) 2011
 * Claudio Kopper <claudio.kopper@nikhef.nl>
 * $Id$
 *
 * @file I3ShadowedPhotonRemoverModule.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 *
 *
 *  This file is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>
 *  
 */

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif 
#include <inttypes.h>

#include "clsim/I3ShadowedPhotonRemoverModule.h"

#include <boost/foreach.hpp>

#include "clsim/I3Photon.h"

#include "dataclasses/geometry/I3Geometry.h"

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


#ifdef IS_Q_FRAME_ENABLED
void I3ShadowedPhotonRemoverModule::DAQ(I3FramePtr frame)
#else
void I3ShadowedPhotonRemoverModule::Physics(I3FramePtr frame)
#endif
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
        const OMKey &key = it.first;
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
