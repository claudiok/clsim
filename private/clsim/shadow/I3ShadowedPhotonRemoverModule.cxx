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

#include <vector>

#include "clsim/shadow/I3ShadowedPhotonRemoverModule.h"

#include "clsim/shadow/I3ExtraGeometryItemCylinder.h"

#include <boost/foreach.hpp>

#include "simclasses/I3CompressedPhoton.h"

#include "dataclasses/I3Double.h"

//#include "dataclasses/geometry/I3Geometry.h"

#include "dataclasses/I3Constants.h"

#include "simclasses/I3CylinderMap.h"

// The module
I3_MODULE(I3ShadowedPhotonRemoverModule);


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

    AddParameter("Cable" , 
		 "Cable that is represented as a cylinder" , 
		 cylinder_name_ );

    AddParameter("Cable_Map",
		 "Map containing all the cables found in the geometry",
		 cylinder_map_name_);

    distance = 10.0;

    AddParameter("Distance" ,
		 "Distance from where photon hits DOM to extended distance to last scatter" , 
		 distance) ;

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
    GetParameter("Cable" , cylinder_name_);
    GetParameter("Cable_Map" , cylinder_map_name_);
    GetParameter("Distance",distance);

    // set up the worker class
}


void I3ShadowedPhotonRemoverModule::DAQ(I3FramePtr frame)
{
    log_trace("%s", __PRETTY_FUNCTION__);
    
    //const I3Geometry& geometry = frame->Get<I3Geometry>();



    I3CompressedPhotonSeriesMapConstPtr inputCompressedPhotonSeriesMap = frame->Get<I3CompressedPhotonSeriesMapConstPtr>(inputPhotonSeriesMapName_);
    if (!inputCompressedPhotonSeriesMap) log_fatal("Frame does not contain an I3CompressedPhotonSeriesMap named \"%s\".",
                                         inputPhotonSeriesMapName_.c_str());
    
    // allocate the output hitSeriesMap
    I3CompressedPhotonSeriesMapPtr outputCompressedPhotonSeriesMap(new I3CompressedPhotonSeriesMap());
    
    BOOST_FOREACH(const I3CompressedPhotonSeriesMap::value_type &it, *inputCompressedPhotonSeriesMap)
    {
        const ModuleKey &key = it.first;
        const I3CompressedPhotonSeries &photons = it.second;

        //// Find the current OM in the geometry map
        //I3OMGeoMap::const_iterator geo_it = geometry.omgeo.find(key);
        //if (geo_it == geometry.omgeo.end())
        //        log_fatal("OM (%i/%u) not found in the current geometry map!", key.GetString(), key.GetOM());
        //const I3OMGeo &om = geo_it->second;

        
        // a pointer to the output vector. The vector will be allocated 
        // by the map, this is merely a pointer to it in case we have multiple
        // photons per OM.
        I3CompressedPhotonSeries *out_photons = NULL;

        BOOST_FOREACH(const I3CompressedPhoton &photon, photons)
        {
            
            I3CompressedPhoton out_photon = photon;
            const bool isShadowed = 
            shadowedPhotonRemover_->IsPhotonShadowed(out_photon);

            if (isShadowed) continue;

            // allocate the output vector if not already done
            if (!out_photons) out_photons = &(outputCompressedPhotonSeriesMap->insert(std::make_pair(key, I3CompressedPhotonSeries())).first->second);

            // add a new copy of the input photon to the output list
            out_photons->push_back(out_photon);
        }
        
        
    }

    // store the output I3MCHitSeriesMap
    frame->Put(outputPhotonSeriesMapName_, outputCompressedPhotonSeriesMap);
    
    // that's it!
    PushFrame(frame);
}

void I3ShadowedPhotonRemoverModule::Geometry(I3FramePtr frame)
{
  if (!cylinder_map_name_.empty() && !cylinder_name_.empty() )
    log_fatal("Both Cylinder and Cylinder Map Name were specified.");

  if (cylinder_map_name_.empty() && cylinder_name_.empty() )
    log_fatal("Both Cylinder and Cyilnder Map were not specified.");

  log_trace("%s", __PRETTY_FUNCTION__);
  
  if(!cylinder_name_.empty()){
    const I3ExtraGeometryItemCylinder& cylinder = frame ->Get<I3ExtraGeometryItemCylinder>(cylinder_name_);
    shadowedPhotonRemover_ = I3ShadowedPhotonRemoverPtr(new I3ShadowedPhotonRemover(cylinder , distance ));
  }
  else if (!cylinder_map_name_.empty()){
    const I3CylinderMap& cylinder_map = frame ->Get<I3CylinderMap>(cylinder_map_name_);
    I3CylinderMap::const_iterator it2 = cylinder_map.begin();
      while (it2 != cylinder_map.end()){
      shadowedPhotonRemover_ = I3ShadowedPhotonRemoverPtr(new I3ShadowedPhotonRemover(it2->second , distance));
    }
  }

  PushFrame(frame);
};
  
