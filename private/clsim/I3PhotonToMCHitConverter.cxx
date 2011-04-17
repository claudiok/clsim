/**
 * copyright (C) 2011
 * Claudio Kopper <claudio.kopper@nikhef.nl>
 * $Id$
 *
 * @file I3PhotonToMCHitConverter.cxx
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

#define __STDC_FORMAT_MACROS 
#include <inttypes.h>

#include "clsim/I3PhotonToMCHitConverter.h"

#include <boost/foreach.hpp>

#include "clsim/I3Photon.h"

#include "dataclasses/geometry/I3Geometry.h"
#include "dataclasses/physics/I3MCHit.h"
#include "dataclasses/physics/I3MCTree.h"

#include "dataclasses/I3Constants.h"

// The module
I3_MODULE(I3PhotonToMCHitConverter);

I3PhotonToMCHitConverter::I3PhotonToMCHitConverter(const I3Context& context) 
: I3ConditionalModule(context)
{
    AddParameter("RandomService",
                 "A random number generating service (derived from I3RandomService).",
                 randomService_);

    inputPhotonSeriesMapName_="PropagatedPhotons";
    AddParameter("InputPhotonSeriesMapName",
                 "Name of the input I3PhotonSeriesMap frame object.",
                 inputPhotonSeriesMapName_);

    outputMCHitSeriesMapName_="MCHitSeriesMap";
    AddParameter("OutputMCHitSeriesMapName",
                 "Name of the output I3MCHitSeries frame object. ",
                 outputMCHitSeriesMapName_);

    MCTreeName_="I3MCTree";
    AddParameter("MCTreeName",
                 "Name of the I3MCTree frame object. All photon particle IDs are checked against this tree.",
                 MCTreeName_);

    AddParameter("WavelengthAcceptance",
                 "Wavelength acceptance of the (D)OM as a I3WlenDependedValue object.",
                 wavelengthAcceptance_);

    AddParameter("AngularAcceptance",
                 "Angular acceptance of the (D)OM as a I3WlenDependedValue object.",
                 angularAcceptance_);


    // add an outbox
    AddOutBox("OutBox");

}

I3PhotonToMCHitConverter::~I3PhotonToMCHitConverter()
{
    log_trace("%s", __PRETTY_FUNCTION__);
}


void I3PhotonToMCHitConverter::Configure()
{
    log_trace("%s", __PRETTY_FUNCTION__);

    GetParameter("RandomService", randomService_);

    GetParameter("InputPhotonSeriesMapName", inputPhotonSeriesMapName_);
    GetParameter("OutputMCHitSeriesMapName", outputMCHitSeriesMapName_);

    GetParameter("MCTreeName", MCTreeName_);

    GetParameter("WavelengthAcceptance", wavelengthAcceptance_);
    GetParameter("AngularAcceptance", angularAcceptance_);

    if (!wavelengthAcceptance_)
        log_fatal("The \"WavelengthAcceptance\" parameter must not be empty.");
    if (!angularAcceptance_)
        log_fatal("The \"AngularAcceptance\" parameter must not be empty.");
    
    if (!wavelengthAcceptance_->HasNativeImplementation())
        log_fatal("The wavelength acceptance function must have a native (i.e. non-OpenCL) implementation!");
    if (!angularAcceptance_->HasNativeImplementation())
        log_fatal("The angular acceptance function must have a native (i.e. non-OpenCL) implementation!");
    
    
    if (!randomService_) {
        log_info("No random service provided as a parameter, trying to get one from the context..");
        randomService_ = context_.Get<I3RandomServicePtr>();
    }
    
    if (!randomService_)
        log_fatal("No random service provided using the \"RandomService\" parameter (and there is none installed on the context).");
    
}


void I3PhotonToMCHitConverter::Physics(I3FramePtr frame)
{
    log_trace("%s", __PRETTY_FUNCTION__);
    
    // First we need to get our geometry
	const I3Geometry& geometry = frame->Get<I3Geometry>();

    I3PhotonSeriesMapConstPtr inputPhotonSeriesMap = frame->Get<I3PhotonSeriesMapConstPtr>(inputPhotonSeriesMapName_);
    if (!inputPhotonSeriesMap) log_fatal("Frame does not contain an I3PhotonSeriesMap named \"%s\".",
                                         inputPhotonSeriesMapName_.c_str());
    
    // urrently, the only reason we need the MCTree is that I3MCHit does
    // only allow setting the major/minor particle IDs using an existing
    // I3Particle instance with that ID combination.
    I3MCTreeConstPtr MCTree = frame->Get<I3MCTreeConstPtr>(MCTreeName_);
    if (!MCTree) log_fatal("Frame does not contain an I3MCTree named \"%s\".",
                           MCTreeName_.c_str());

    // allocate the output hitSeriesMap
    I3MCHitSeriesMapPtr outputMCHitSeriesMap(new I3MCHitSeriesMap());
    
    // build an index into the I3MCTree
    std::map<std::pair<uint64_t, int>, const I3Particle &> mcTreeIndex;
    for (I3MCTree::iterator it = MCTree->begin();
         it != MCTree->end(); ++it)
    {
        const I3Particle &particle = *it;
        mcTreeIndex.insert(std::make_pair(std::make_pair(particle.GetMajorID(), particle.GetMinorID()), particle));
    }
    
    
    
    // DOM is looking downwards (this is currently specific to IceCube)
    const double DOMDir_x = 0.;
    const double DOMDir_y = 0.;
    const double DOMDir_z = -1.;
    
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
        // hits per OM.
        I3MCHitSeries *hits = NULL;

        BOOST_FOREACH(const I3Photon &photon, photons)
        {
            double hitProbability = photon.GetWeight();
            if (hitProbability < 0.) log_fatal("Photon with negative weight found.");
            if (hitProbability == 0.) continue;

            double photonCosAngle = -(photon.GetDir().GetX() * DOMDir_x +
                                      photon.GetDir().GetY() * DOMDir_y +
                                      photon.GetDir().GetZ() * DOMDir_z);
            photonCosAngle = std::max(-1., std::min(1., photonCosAngle));
            
#ifndef I3_OPTIMIZE
            const double photonAngle = std::acos(photonCosAngle);
            const double distFromDOMCenter = std::sqrt(std::pow(photon.GetPos().GetX()-om.position.GetX(),2) + 
                                                       std::pow(photon.GetPos().GetY()-om.position.GetY(),2) + 
                                                       std::pow(photon.GetPos().GetZ()-om.position.GetZ(),2));
            
            log_trace("Photon (lambda=%fnm, angle=%fdeg, dist=%fm) has weight %g",
                     photon.GetWavelength()/I3Units::nanometer,
                     photonAngle/I3Units::deg,
                     distFromDOMCenter/I3Units::m,
                     hitProbability);
#endif
            
            hitProbability *= wavelengthAcceptance_->GetValue(photon.GetWavelength());
            log_trace("After wlen acceptance: prob=%g (wlen acceptance is %f)",
                     hitProbability, wavelengthAcceptance_->GetValue(photon.GetWavelength()));

            hitProbability *= angularAcceptance_->GetValue(photonCosAngle);
            log_trace("After wlen&angular acceptance: prob=%g (angular acceptance is %f)",
                      hitProbability, angularAcceptance_->GetValue(photonCosAngle));

            // does it survive?
            if (hitProbability <= randomService_->Uniform()) continue;

            // find the particle
            std::map<std::pair<uint64_t, int>, const I3Particle &>::const_iterator it = 
            mcTreeIndex.find(std::make_pair(photon.GetParticleMajorID(), photon.GetParticleMinorID()));
            if (it==mcTreeIndex.end())
                log_fatal("Particle with id maj=%" PRIu64 ", min=%i does not exist in MC tree, but we have a photon that claims it was created by that particle..",
                          photon.GetParticleMajorID(), photon.GetParticleMinorID());
            const I3Particle &particle = it->second;
            
            // allocate the output vector if not already done
            if (!hits) hits = &(outputMCHitSeriesMap->insert(std::make_pair(key, I3MCHitSeries())).first->second);

            // add a new hit
            hits->push_back(I3MCHit());
            I3MCHit &hit = hits->back();
            
            // fill in all information
            hit.SetTime(photon.GetTime());
            hit.SetHitID(photon.GetID());
            hit.SetWeight(1.);
            hit.SetParticleID(particle);
            hit.SetCherenkovDistance(NAN);
            hit.SetHitSource(I3MCHit::SPE); // SPE for now, may be changed by afterpulse simulation

        }
        
        
    }
    
    // store the output I3MCHitSeriesMap
    frame->Put(outputMCHitSeriesMapName_, outputMCHitSeriesMap);
    
    // that's it!
    PushFrame(frame);
}
