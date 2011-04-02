/**
 * copyright (C) 2011
 * Claudio Kopper <claudio.kopper@nikhef.nl>
 * $Id$
 *
 * @file I3Photon.h
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

#ifndef I3PHOTON_H_INCLUDED
#define I3PHOTON_H_INCLUDED

#include "icetray/I3FrameObject.h"

#include "dataclasses/I3Vector.h"
#include "dataclasses/I3Map.h"
#include "dataclasses/OMKey.h"
#include "dataclasses/I3Direction.h"
#include "dataclasses/I3Position.h"
#include "dataclasses/physics/I3Particle.h"

/**
 * @brief This class contains a photon with
 * simulated arrival time, the direction
 * the photon arrived from, the position
 * on the OM the photon hit (stored as an 
 * direction from the OM center to the hit
 * position) and the photon's wavelength.
 */
static const unsigned i3photon_version_ = 0;

class I3Photon : public I3FrameObject
{
    
public:
    
    I3Photon() : 
    time_(NAN), 
    weight_(NAN), 
    ID_(-1),
    particleID_(-1), 
    particleMajorID_(0), 
    cherenkovDist_(NAN),
    cherenkovTime_(NAN),
    wavelength_(NAN),
    direction_(),
    position_(),
    numScattered_(0)
    { }
    
    I3Photon(uint64_t mid, int id) : 
    time_(NAN), 
    weight_(NAN), 
    ID_(-1),
    particleID_(id), 
    particleMajorID_(mid), 
    cherenkovDist_(NAN),
    cherenkovTime_(NAN),
    wavelength_(NAN),
    direction_(),
    position_(),
    numScattered_(0)
    { }
    
    virtual ~I3Photon();
    
    double GetTime() const {return time_;}
    
    void SetTime(double time) {time_ = time;}
    
    int32_t GetID() const {return ID_;}
    
    void SetID(int32_t ID) {ID_ = ID;}
    
    double GetWeight() const {return weight_; }
    
    void SetWeight(double weight) {weight_ = weight; }
    
    int32_t GetParticleMinorID() const { return particleID_; }
    
    uint64_t GetParticleMajorID() const { return particleMajorID_; }

    void SetParticleMinorID(int32_t minorID) { particleID_ = minorID; }
    
    void SetParticleMajorID(uint64_t majorID) { particleMajorID_ = majorID; }

    void SetParticleID(const I3Particle&);
    
    /**
     * @return the path distance to the track which emitted this photon.
     * This is the full path with all possible scatters.
     */
    double GetCherenkovDist() const {return cherenkovDist_; }
    
    /**
     * @param CherenkovDist set the path distance to track which 
     * emitted this photon.
     * This is the full path with all possible scatters.
     */
    void SetCherenkovDist(double CherenkovDist) {cherenkovDist_ = CherenkovDist; }
    
    /**
     * @return The full time it took the photon to travel from its emission
     * point to the OM.
     */
    double GetCherenkovTime() const {return cherenkovTime_; }
    
    /**
     * @param CherenkovTime The full time it took the photon to travel from its emission
     * point to the OM.
     */
    void SetCherenkovTime(double CherenkovTime) {cherenkovTime_ = CherenkovTime; }
    
    double GetWavelength() const {return wavelength_;}
    
    void SetWavelength(double wlen) {wavelength_ = wlen;}
    
    /**
     * @return The position of this photon on the OM. This
     * position is relative to the global coordinate system, where the
     * z-axis is facing upwards. It is NOT relative to the OM axis!!
     */
    const I3Position& GetPos() const { return position_; }
    
    /**
     * @param p The position of this on the OM. This
     * position is relative to the global coordinate system, where the
     * z-axis is facing upwards. It is NOT relative to the OM axis!!
     */
    void SetPos(const I3Position& p) { position_.SetPosition(p); }
    
    /**
     * @return The direction vector of the photon on the OM. 
     * This direction is relative to the global coordinate system, where the
     * z-axis is facing upwards. It is NOT relative to the OM axis!!
     */
    const I3Direction& GetDir() const { return direction_; }
    
    /**
     * @param d The direction vector of the photon on the OM. 
     * This direction is relative to the global coordinate system, where the
     * z-axis is facing upwards. It is NOT relative to the OM axis!!
     */
    void SetDir(const I3Direction& d) { direction_.SetDirection(d); }
    
    /** 
     * @return this returns the number of times this photon was scattered.
     */
    inline uint32_t GetNumScattered() const {return numScattered_;}
    
    /** 
     * this sets the number of times this photon was scattered.
     */
    inline void SetNumScattered(uint32_t numScattered) {numScattered_=numScattered;}
    
    bool operator==(const I3Photon& rhs) {
        return time_ == rhs.time_
        && ID_ == rhs.ID_
        && weight_ == rhs.weight_
        && particleID_ == rhs.particleID_
        && particleMajorID_ == rhs.particleMajorID_
        && cherenkovDist_ == rhs.cherenkovDist_
        && cherenkovTime_ == rhs.cherenkovTime_
        && wavelength_ == rhs.wavelength_
        && direction_.GetAzimuth() == rhs.direction_.GetAzimuth()
        && direction_.GetZenith() == rhs.direction_.GetZenith()
        && position_.GetX() == rhs.position_.GetX()
        && position_.GetY() == rhs.position_.GetY()
        && position_.GetZ() == rhs.position_.GetZ()
        && numScattered_ == rhs.numScattered_;
    }
    
private:
    double time_;
    double weight_;
    int32_t ID_;
    int32_t particleID_;
    uint64_t particleMajorID_;
    double cherenkovDist_;
    double cherenkovTime_;
    double wavelength_;
    
    I3Direction direction_;
    I3Position position_;

    uint32_t numScattered_;

    friend class boost::serialization::access;
    
    template <class Archive> void serialize(Archive & ar, unsigned version);
    
};

BOOST_CLASS_VERSION(I3Photon, i3photon_version_);

typedef I3Vector<I3Photon> I3PhotonSeries;
typedef I3Map<OMKey, I3PhotonSeries> I3PhotonSeriesMap; 

I3_POINTER_TYPEDEFS(I3Photon);
I3_POINTER_TYPEDEFS(I3PhotonSeries);
I3_POINTER_TYPEDEFS(I3PhotonSeriesMap);

#endif //I3PHOTON_H_INCLUDED
