/**
 * copyright (C) 2011
 * Claudio Kopper <claudio.kopper@nikhef.nl>
 * $Id$
 *
 * @file I3HoleIceSimulatorModule.cxx
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

#include "clsim/I3HoleIceSimulatorModule.h"

#include <boost/foreach.hpp>

#include "clsim/I3Photon.h"

#include "dataclasses/geometry/I3Geometry.h"
#include "dataclasses/physics/I3MCHit.h"
#include "dataclasses/physics/I3MCTree.h"

#include "dataclasses/I3Constants.h"

// The module
I3_MODULE(I3HoleIceSimulatorModule);

I3HoleIceSimulatorModule::I3HoleIceSimulatorModule(const I3Context& context) 
: I3ConditionalModule(context)
{
    AddParameter("RandomService",
                 "A random number generating service (derived from I3RandomService).",
                 randomService_);

    inputPhotonSeriesMapName_="PropagatedPhotons";
    AddParameter("InputPhotonSeriesMapName",
                 "Name of the input I3PhotonSeriesMap frame object.",
                 inputPhotonSeriesMapName_);

    outputPhotonSeriesMapName_="PropagatedPhotonsAfterHoleIce";
    AddParameter("OutputPhotonSeriesMapName",
                 "Name of the output I3PhotonSeriesMap frame object.",
                 outputPhotonSeriesMapName_);

    DOMRadiusWithoutOversize_=0.16510*I3Units::m; // 13 inch diameter
    AddParameter("DOMRadiusWithoutOversize",
                 "Specifiy the DOM radius. Do not include oversize factors here.",
                 DOMRadiusWithoutOversize_);

    holeRadius_=0.6*I3Units::m;
    AddParameter("HoleRadius",
                 "Specifiy the hole radius.",
                 holeRadius_);

    AddParameter("MediumProperties",
                 "An instance of I3CLSimMediumProperties describing the bulk ice properties.",
                 mediumProperties_);

    AddParameter("HoleIceAbsorptionLength",
                 "An instance of I3CLSimWlenDependentValue describing the hole ice absorption lengths.",
                 holeIceAbsorptionLength_);

    AddParameter("HoleIceScatteringLength",
                 "An instance of I3CLSimMediumProperties describing the hole ice scattering lengths.",
                 holeIceScatteringLength_);

    // add an outbox
    AddOutBox("OutBox");

}

I3HoleIceSimulatorModule::~I3HoleIceSimulatorModule()
{
    log_trace("%s", __PRETTY_FUNCTION__);
}


void I3HoleIceSimulatorModule::Configure()
{
    log_trace("%s", __PRETTY_FUNCTION__);

    GetParameter("RandomService", randomService_);

    GetParameter("InputPhotonSeriesMapName", inputPhotonSeriesMapName_);
    GetParameter("OutputPhotonSeriesMapName", outputPhotonSeriesMapName_);

    GetParameter("DOMRadiusWithoutOversize", DOMRadiusWithoutOversize_);
    GetParameter("HoleRadius", holeRadius_);
    GetParameter("MediumProperties", mediumProperties_);
    GetParameter("HoleIceAbsorptionLength", holeIceAbsorptionLength_);
    GetParameter("HoleIceScatteringLength", holeIceScatteringLength_);

    if (holeRadius_ < DOMRadiusWithoutOversize_)
        log_fatal("DOM does not fit in hole.");
    
    if (!mediumProperties_) log_fatal("You have to specify the \"MediumProperties\" parameter!");
    if (!holeIceAbsorptionLength_) log_fatal("You have to specify the \"HoleIceAbsorptionLength\" parameter!");
    if (!holeIceScatteringLength_) log_fatal("You have to specify the \"HoleIceScatteringLength\" parameter!");
    
    if (!randomService_) {
        log_info("No random service provided as a parameter, trying to get one from the context..");
        randomService_ = context_.Get<I3RandomServicePtr>();
    }
    
    if (!randomService_)
        log_fatal("No random service provided using the \"RandomService\" parameter (and there is none installed on the context).");


    // set up the hole ice simulator
    holeIceSimulator_ = I3HoleIceSimulatorPtr(new I3HoleIceSimulator(randomService_,
                                                                     DOMRadiusWithoutOversize_,
                                                                     holeRadius_,
                                                                     mediumProperties_,
                                                                     holeIceAbsorptionLength_,
                                                                     holeIceScatteringLength_));

}


#ifdef IS_Q_FRAME_ENABLED
void I3HoleIceSimulatorModule::DAQ(I3FramePtr frame)
#else
void I3HoleIceSimulatorModule::Physics(I3FramePtr frame)
#endif
{
    log_trace("%s", __PRETTY_FUNCTION__);
    
    // First we need to get our geometry
    const I3Geometry& geometry = frame->Get<I3Geometry>();

    I3PhotonSeriesMapConstPtr inputPhotonSeriesMap = frame->Get<I3PhotonSeriesMapConstPtr>(inputPhotonSeriesMapName_);
    if (!inputPhotonSeriesMap) log_fatal("Frame does not contain an I3PhotonSeriesMap named \"%s\".",
                                         inputPhotonSeriesMapName_.c_str());
    
    // allocate the output hitSeriesMap
    I3PhotonSeriesMapPtr outputPhotonSeriesMap(new I3PhotonSeriesMap());
    
    BOOST_FOREACH(const I3PhotonSeriesMap::value_type &it, *inputPhotonSeriesMap)
    {
        const OMKey &key = it.first;
        const I3PhotonSeries &photons = it.second;

        // Find the current OM in the geometry map
        I3OMGeoMap::const_iterator geo_it = geometry.omgeo.find(key);
        if (geo_it == geometry.omgeo.end())
                log_fatal("OM (%i/%u) not found in the current geometry map!", key.GetString(), key.GetOM());
        const I3OMGeo &om = geo_it->second;

        
        // a pointer to the output vector. The vector will be allocated 
        // by the map, this is merely a pointer to it in case we have multiple
        // photons per OM.
        I3PhotonSeries *out_photons = NULL;

        BOOST_FOREACH(const I3Photon &photon, photons)
        {
            
            I3Photon out_photon = photon;
            const bool arrivedAtDom = 
            holeIceSimulator_->TrackPhoton(out_photon, om.position);

            if (!arrivedAtDom) continue;

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
