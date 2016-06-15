/**
 * Copyright (c) 2013
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
 * $Id: I3CompressedPhoton.h$
 *
 * @file I3CompressedPhoton.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3COMPRESSEDPHOTON_H_INCLUDED
#define I3COMPRESSEDPHOTON_H_INCLUDED

#include <stdexcept>

#include "icetray/I3FrameObject.h"

#include "dataclasses/I3Vector.h"
#include "dataclasses/I3Map.h"
#include "dataclasses/ModuleKey.h"
#include "dataclasses/I3Direction.h"
#include "dataclasses/I3Position.h"
#include "dataclasses/physics/I3Particle.h"

#include "clsim/I3Photon.h"

/**
 * @brief A photon with only the most important
 * information.
 */
static const unsigned i3compressedphoton_version_ = 0;

class I3CompressedPhoton
{
    
public:
    I3CompressedPhoton() : 
    time_(NAN), 
    weight_(NAN), 
    wavelength_(NAN),
    zenith_(NAN),
    azimuth_(NAN),
    x_(NAN),
    y_(NAN),
    z_(NAN),
    particleMajorID_(0),
    particleID_(-1)
    { }
    
    I3CompressedPhoton(uint64_t mid, int id) : 
    time_(NAN), 
    weight_(NAN), 
    wavelength_(NAN),
    zenith_(NAN),
    azimuth_(NAN),
    x_(NAN),
    y_(NAN),
    z_(NAN),
    particleMajorID_(mid),
    particleID_(id)
    { }
    
    I3CompressedPhoton(const I3Photon &p) : 
    time_(p.GetTime()),
    weight_(p.GetWeight()), 
    wavelength_(p.GetWavelength()),
    zenith_(p.GetDir().GetZenith()),
    azimuth_(p.GetDir().GetAzimuth()),
    x_(p.GetPos().GetX()),
    y_(p.GetPos().GetY()),
    z_(p.GetPos().GetZ()),
    particleMajorID_(p.GetParticleMajorID()),
    particleID_(p.GetParticleMinorID())
    { }

    operator I3Photon() const
    {
        I3Photon p(particleMajorID_, particleID_);
        // I3Photon p;

        p.SetTime(time_);
        p.SetWeight(weight_);
        p.SetWavelength(wavelength_);
        p.SetDir(I3Direction(zenith_, azimuth_));
        p.SetPos(I3Position(x_, y_, z_));

        return p;
    }

    ~I3CompressedPhoton();
    
    double GetTime() const {return time_;}
    
    void SetTime(double time) {time_ = time;}
    
    float GetWeight() const {return weight_; }
    
    void SetWeight(float weight) {weight_ = weight; }
    
    int32_t GetParticleMinorID() const { return particleID_; }
    uint64_t GetParticleMajorID() const { return particleMajorID_; }

    void SetParticleMinorID(int32_t minorID) { particleID_ = minorID; }
    void SetParticleMajorID(uint64_t majorID) { particleMajorID_ = majorID; }

    void SetParticleID(const I3Particle&);
    
    float GetWavelength() const {return wavelength_;}
    
    void SetWavelength(float wlen) {wavelength_ = wlen;}

    /**
     * @return The position of this photon on the OM. This
     * position is relative to the global coordinate system, where the
     * z-axis is facing upwards. It is NOT relative to the OM axis!!
     */
    const I3Position GetPos() const { return I3Position(x_, y_, z_); }
    
    /**
     * @param p The position of this on the OM. This
     * position is relative to the global coordinate system, where the
     * z-axis is facing upwards. It is NOT relative to the OM axis!!
     */
    void SetPos(const I3Position& p) { x_=p.GetX(); y_=p.GetY(); z_=p.GetZ(); }
    
    /**
     * @return The direction vector of the photon on the OM. 
     * This direction is relative to the global coordinate system, where the
     * z-axis is facing upwards. It is NOT relative to the OM axis!!
     */
    const I3Direction GetDir() const { return I3Direction(zenith_, azimuth_); }
    
    /**
     * @param d The direction vector of the photon on the OM. 
     * This direction is relative to the global coordinate system, where the
     * z-axis is facing upwards. It is NOT relative to the OM axis!!
     */
    void SetDir(const I3Direction& d) { zenith_=d.GetZenith(); azimuth_=d.GetAzimuth(); }

    
    bool operator==(const I3CompressedPhoton& rhs) const {
        if (!( time_ == rhs.time_
        && weight_ == rhs.weight_
        && wavelength_ == rhs.wavelength_
        && azimuth_ == rhs.azimuth_
        && zenith_ == rhs.zenith_
        && x_ == rhs.x_
        && y_ == rhs.y_
        && z_ == rhs.z_
        && particleID_ == rhs.particleID_
        && particleMajorID_ == rhs.particleMajorID_))
            return false;
        
        return true;
    }
    
private:
    double time_;
    float weight_;
    float wavelength_;
    
    float zenith_;
    float azimuth_;
    double x_;
    double y_;
    double z_;

    uint64_t particleMajorID_;
    int32_t particleID_;

    friend class icecube::serialization::access;
    template <class Archive> void serialize(Archive & ar, unsigned version);
    
};

I3_CLASS_VERSION(I3CompressedPhoton, i3compressedphoton_version_);

typedef I3Vector<I3CompressedPhoton> I3CompressedPhotonSeries;

typedef I3Map<ModuleKey, I3CompressedPhotonSeries> I3CompressedPhotonSeriesMap; 

I3_POINTER_TYPEDEFS(I3CompressedPhoton);
I3_POINTER_TYPEDEFS(I3CompressedPhotonSeries);
I3_POINTER_TYPEDEFS(I3CompressedPhotonSeriesMap);

#endif //I3COMPRESSEDPHOTON_H_INCLUDED
