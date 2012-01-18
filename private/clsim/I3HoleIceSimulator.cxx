/**
 * copyright (C) 2011
 * Claudio Kopper <claudio.kopper@nikhef.nl>
 * $Id$
 *
 * @file I3HoleIceSimulator.cxx
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

#include <limits>

#include "clsim/I3HoleIceSimulator.h"

#include "dataclasses/I3Constants.h"

I3HoleIceSimulator::I3HoleIceSimulator(I3RandomServicePtr random,
                                       double DOMRadius, double holeRadius,
                                       I3CLSimMediumPropertiesConstPtr mediumProperties,
                                       I3CLSimWlenDependentValueConstPtr holeIceAbsorptionLength,
                                       I3CLSimWlenDependentValueConstPtr holeIceScatteringLength) 
:
random_(random),
DOMRadius_(DOMRadius), holeRadius_(holeRadius),
mediumProperties_(mediumProperties),
holeIceAbsorptionLength_(holeIceAbsorptionLength),
holeIceScatteringLength_(holeIceScatteringLength)
{
    log_trace("%s", __PRETTY_FUNCTION__);

}

I3HoleIceSimulator::~I3HoleIceSimulator()
{
    log_trace("%s", __PRETTY_FUNCTION__);
}


#define EPSILON 1e-5


namespace {
    inline double lengthToNearestIntersectionWithCylinder(double px, double py, double pz,
                                                          double dx, double dy, double dz,
                                                          double radius, bool &currentlyInside,
                                                          double &lengthToNextToNearestIntersection)
    {
        const double px2 = px*px;
        const double py2 = py*py;
        const double radius2 = radius*radius;
        currentlyInside = px2+py2 <= radius2;

        const double dx2 = dx*dx;
        const double dy2 = dy*dy;
        const double dz2 = dz*dz;
        
        if (std::abs(1.-(dx2+dy2+dz2)) > 1e-10) log_fatal("direction vector not normalized!");
        if (std::abs(1.-dz) < 1e-10) {
            log_trace("no intersection for vertical tracks");
            lengthToNextToNearestIntersection = std::numeric_limits<double>::infinity();
            return std::numeric_limits<double>::infinity(); // no intersection for vertical tracks
        }
        
        const double pxdx = px*dx;
        const double pydy = py*dy;

        const double D = 2.*pxdx*pydy + dx2*(radius2-py2) + dy2*(radius2-px2);
        if (D<0.) {
            log_trace("can't be inside if there's no intersection (for a non-vertical track)");
            currentlyInside=false; // can't be inside if there's no intersection (for a non-vertical track)
            lengthToNextToNearestIntersection = std::numeric_limits<double>::infinity();
            return std::numeric_limits<double>::infinity(); // no intersection
        }
        
        const double sqrtD = std::sqrt(D);
        
        const double A = -pxdx-pydy;
        const double B = 1.-dz2;
        
        const double lambda1 = (A+sqrtD)/B;
        const double lambda2 = (A-sqrtD)/B;
        
        if (lambda1<=0.) {
            lengthToNextToNearestIntersection = std::numeric_limits<double>::infinity();

            if (lambda2<=0.) {
                log_trace("can't be inside if there's no intersection (for a non-vertical track)");
                currentlyInside=false; // can't be inside if there's no intersection (for a non-vertical track)
                return std::numeric_limits<double>::infinity(); // no intersection
            }
            // lambda2 > 0.:
            
            return lambda2;
        } else {
            // lambda1 > 0.
            if (lambda2 <= 0.) {
                lengthToNextToNearestIntersection = std::numeric_limits<double>::infinity();
                return lambda1;
            }
            
            if (lambda1 < lambda2)
            {
                lengthToNextToNearestIntersection = lambda2;
                return lambda1;
            }
            else
            {
                lengthToNextToNearestIntersection = lambda1;
                return lambda2;
            }
            
            return std::min(lambda1, lambda2); // return the first of two intersections
        }
        
    }
    
}

bool I3HoleIceSimulator::TrackPhoton(I3Photon &photon,
                                     const I3Position &om_pos)
{
    // determine current ice layer
    
    // determine the photon layer per DOM only
    // assume that the photon will only be in the layer closest to the measuring DOM
    int64_t photonLayerRaw = static_cast<int64_t>((om_pos.GetZ()-mediumProperties_->GetLayersZStart())/mediumProperties_->GetLayersHeight());
    if (photonLayerRaw < 0) photonLayerRaw=0;
    if (photonLayerRaw >= static_cast<int64_t>(mediumProperties_->GetLayersNum())) photonLayerRaw=mediumProperties_->GetLayersNum()-1;
    
    const uint32_t photonLayer = static_cast<uint32_t>(photonLayerRaw);

    // look up ice properties for current photon wavelength and layer
    const double wavelength = photon.GetWavelength();

    const double bulkAbsorptionLength = mediumProperties_->GetAbsorptionLength(photonLayer)->GetValue(wavelength);
    const double bulkScatteringLength = mediumProperties_->GetScatteringLength(photonLayer)->GetValue(wavelength);

    const double holeAbsorptionLength = holeIceAbsorptionLength_->GetValue(wavelength);
    const double holeScatteringLength = holeIceScatteringLength_->GetValue(wavelength);

    
    double groupVelocity=NAN;
    
    if (mediumProperties_->GetGroupRefractiveIndexOverride(photonLayer)) {
        groupVelocity = I3Constants::c / mediumProperties_->GetGroupRefractiveIndexOverride(photonLayer)->GetValue(wavelength);
    } else {
        const double phaseRefIndex = mediumProperties_->GetPhaseRefractiveIndex(photonLayer)->GetValue(wavelength);
        const double n_inv = 1./phaseRefIndex;

        if (!mediumProperties_->GetPhaseRefractiveIndex(photonLayer)->HasDerivative())
            log_fatal("Need derivative of phase refractive index, which is not available!");
        
        const double dispersion = mediumProperties_->GetPhaseRefractiveIndex(photonLayer)->GetDerivative(wavelength);
        
        groupVelocity = I3Constants::c * (1.0 + dispersion*wavelength*n_inv) * n_inv;
    }

    if (isnan(photon.GetGroupVelocity()) || (photon.GetGroupVelocity() <= 0.))
    {
        photon.SetGroupVelocity(groupVelocity);
    }
    else
    {   // paranoid sanity check
        if (std::abs(photon.GetGroupVelocity() - groupVelocity) > EPSILON)
            log_fatal("Group velocity of photon is different from group velocity defined by ice properties. v_phot=%f, v_ice=%f",
                      photon.GetGroupVelocity(), groupVelocity);
    }
    
    if (isnan(photon.GetCherenkovDist())) photon.SetCherenkovDist(0.);
    if (isnan(photon.GetStartTime())) photon.SetStartTime(photon.GetTime());
    
    double dx=photon.GetDir().GetX();
    double dy=photon.GetDir().GetY();
    double dz=photon.GetDir().GetZ();
    double px=photon.GetPos().GetX()-om_pos.GetX();
    double py=photon.GetPos().GetY()-om_pos.GetY();
    double pz=photon.GetPos().GetZ()-om_pos.GetZ();
    double totalTrackedLength=0.;

    // determine the track lengths in the very beginning
    double abs_lens_left=-std::log(random_->Uniform());
    double sca_lens_left=-std::log(random_->Uniform());

    enum whatToDo_t {doChangeGeometry, doScatter, doAbsorb, doPhotonAtDOM};
    bool photonAtDOM=false;
    bool justChangedGeometry=false; // prevent from changing twice
    uint32_t numScatters=0;
    
    for (;;)
    {
        log_trace("old_pos=(%f,%f,%f)m", px/I3Units::m, py/I3Units::m, pz/I3Units::m);
        
        // find next intersection with hole ice cylinder
        bool currentlyInCylinder;
        double lengthToNextToNearestGeometricBoundary;
        double lengthToNextGeometricBoundary = 
        lengthToNearestIntersectionWithCylinder(px,py,pz,dx,dy,dz,holeRadius_,currentlyInCylinder,lengthToNextToNearestGeometricBoundary);
        
        // prevent frequent geometry changes at the same boundary due to limited precision
        if ((justChangedGeometry) && (lengthToNextGeometricBoundary < EPSILON))
        {
            log_trace("just changed geometry, will not change it again at d=%fm, advancing to next boundary at d=%fm",
                      lengthToNextGeometricBoundary/I3Units::m, lengthToNextToNearestGeometricBoundary/I3Units::m);
            lengthToNextGeometricBoundary = lengthToNextToNearestGeometricBoundary;
        }
        justChangedGeometry=false;
        
        log_trace("lengthToNextGeometricBoundary=%fm [%s]", lengthToNextGeometricBoundary/I3Units::m, currentlyInCylinder?"in cyl":"outside");
        
        // paranoid sanity check
        if (std::isfinite(lengthToNextGeometricBoundary))
        {
            const double px2 = px+lengthToNextGeometricBoundary*dx;
            const double py2 = py+lengthToNextGeometricBoundary*dy;
            const double rhoCalc = px2*px2+py2*py2;
            const double rhoExpect = holeRadius_*holeRadius_;
            
            if ( std::abs(rhoCalc-rhoExpect) > EPSILON )
                log_fatal("Internal logic error: should be on cylinder boundary! rhoCalc=%f, rhoExpect=%f",
                          rhoCalc, rhoExpect);
        }
        
        const double currentAbsorptionLength = (currentlyInCylinder?holeAbsorptionLength:bulkAbsorptionLength);
        const double currentScatteringLength = (currentlyInCylinder?holeScatteringLength:bulkScatteringLength);
        
        const double lengthToAbsorption = abs_lens_left*currentAbsorptionLength;
        const double lengthToScatter    = sca_lens_left*currentScatteringLength;

        log_trace("lengthToAbsorption=%fm", lengthToAbsorption/I3Units::m);
        log_trace("lengthToScatter=%fm", lengthToScatter/I3Units::m);

        // the photon will now either:
        //  - be absorbed
        //  - be scattered
        //  - or change from inside to outside the hole ice (or vice versa)

        // determine what to do
        whatToDo_t whatToDo;
        double propagationStepLength;
        if (lengthToScatter < lengthToAbsorption) {
            if (lengthToNextGeometricBoundary < lengthToScatter) {
                whatToDo = doChangeGeometry;
                propagationStepLength=lengthToNextGeometricBoundary;
            } else {
                whatToDo = doScatter;
                propagationStepLength=lengthToScatter;
            }
        } else {
            if (lengthToNextGeometricBoundary < lengthToAbsorption) {
                whatToDo = doChangeGeometry;
                propagationStepLength=lengthToNextGeometricBoundary;
            } else {
                whatToDo = doAbsorb;
                propagationStepLength=lengthToAbsorption;
            }
        }

        // we have a step length now. see if the photon hit anything (the DOM or the cable)
        // (assume DOM is at (0,0,0) )
        {
            const double dr2 = px*px + py*py + pz*pz;
            if (dr2 < DOMRadius_*DOMRadius_)
            {
                // photon is inside the DOM.
                propagationStepLength=0.; // limit step length
                whatToDo = doPhotonAtDOM;
            }
            else
            {
                const double urdot = px*dx + py*dy + pz*dz;
                double discr   = urdot*urdot - dr2 + DOMRadius_*DOMRadius_;   // (discr)^2

            
                if (discr >= 0.0f) 
                {
                    discr = std::sqrt(discr);
                    
                    double smin = -urdot - discr;
                    if (smin < 0.0f) smin = -urdot + discr;
                    
                    // check if distance to intersection <= thisStepLength; if not then no detection 
                    if ((smin >= 0.0f) && (smin < propagationStepLength))
                    {
                        propagationStepLength=smin; // limit step length
                        whatToDo = doPhotonAtDOM;
                    }
                }
            }
        
        }
        
        
        // update photon position
        px += propagationStepLength*dx;
        py += propagationStepLength*dy;
        pz += propagationStepLength*dz;
        totalTrackedLength += propagationStepLength;

        log_trace("new_pos=(%f,%f,%f)m", px/I3Units::m, py/I3Units::m, pz/I3Units::m);

        // update absorption and scattering lengths
        abs_lens_left -= propagationStepLength/currentAbsorptionLength;
        sca_lens_left -= propagationStepLength/currentScatteringLength;

        if (whatToDo==doPhotonAtDOM) {
            log_trace("-> photon @ DOM");

            // sanity check
            photonAtDOM=true;

            break; // that's it, the photon reached the DOM
        }
        
        if (whatToDo==doAbsorb) {
            log_trace("-> photon absorbed");

            // sanity check
            if (std::abs(abs_lens_left) > EPSILON) log_fatal("doAbsorb, however: abs_lens_left==%f", abs_lens_left);
            
            break; // that's it, no need to update the photon, it got absorbed
        }

        if (whatToDo==doChangeGeometry) {
            log_trace("-> photon changed geometry");
            justChangedGeometry=true;
            
            continue; // continue with the next iteration
        }

        // perform photon scattering
        if (whatToDo==doScatter) {
            log_trace("-> photon scattered");

            // sanity check
            if (std::abs(sca_lens_left) > EPSILON) log_fatal("doScatter, however: sca_lens_left==%f", sca_lens_left);
            
            // scatter the photon (find a new photon direction)
            // flat scattering angle distribution
            double cosScatAngle;
            if (currentlyInCylinder) {
                cosScatAngle = (random_->Uniform()*2.) - 1.; // TODO: FIXME!!!!!
            } else {
                cosScatAngle = 1.; // TODO: FIXME!!!!!
            }
            const double sinScatAngle = std::sqrt(1. - cosScatAngle*cosScatAngle);
            
            // randomize direction of scattering (rotation around old direction axis)
            const double b=2.*M_PI*random_->Uniform();
            const double cosb=std::cos(b);
            const double sinb=std::sin(b);
            
            // Rotate new direction into absolute frame of reference 
            const double sinth = std::sqrt(std::max(0., 1.-dz*dz));
            
            if(sinth>0.){  // Current direction not vertical, so rotate 
                const double old_dx = dx;
                const double old_dy = dy;
                const double old_dz = dz;
                
                dx=old_dx*cosScatAngle - ( (old_dy*cosb+old_dz*old_dx*sinb)*sinScatAngle/sinth );
                dy=old_dy*cosScatAngle + ( (old_dx*cosb-old_dz*old_dy*sinb)*sinScatAngle/sinth );
                dz=old_dz*cosScatAngle + sinScatAngle*sinb*sinth;
            }else{         // Current direction is vertical, so this is trivial
                dx=sinScatAngle * cosb;
                dy=sinScatAngle * sinb;
                dz=cosScatAngle * ((dz>0.)?1.:-1.);
            }

            // normalize the length to compensate for numerical errors
            const double recip_length = 1./( std::sqrt(dx*dx + dy*dy + dz*dz) );
            dx *= recip_length;
            dy *= recip_length;
            dz *= recip_length;
            
            // find a new scattering point
            sca_lens_left = -std::log(random_->Uniform());

            numScatters++;
            
            continue;
        }

        log_error("should not reach this point");
    }
    
    // update the photon
    photon.SetDir(I3Direction(dx,dy,dz));
    photon.SetPos(I3Position(px+om_pos.GetX(),py+om_pos.GetY(),pz+om_pos.GetZ()));
    photon.SetCherenkovDist( photon.GetCherenkovDist() + totalTrackedLength );
    photon.SetTime( photon.GetTime() + totalTrackedLength/photon.GetGroupVelocity() );
    photon.SetNumScattered( photon.GetNumScattered() + numScatters );

    
    return photonAtDOM; // photon did or did not hit the DOM
}
