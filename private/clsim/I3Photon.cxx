/**
 * copyright (C) 2011
 * Claudio Kopper <claudio.kopper@nikhef.nl>
 * $Id$
 *
 * @file I3Photon.cxx
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

#include <icetray/serialization.h>
#include <clsim/I3Photon.h>

I3Photon::~I3Photon() { }

void I3Photon::SetParticleID(const I3Particle& p) { 
    particleID_ = p.GetMinorID(); 
    particleMajorID_ = p.GetMajorID();
}

template <class Archive>
void I3Photon::serialize (Archive &ar, unsigned version)
{
    if (version > i3photon_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3Photon class.",version,i3photon_version_);
    
    ar & make_nvp("time", time_);
    ar & make_nvp("weight", weight_);
    ar & make_nvp("ID",ID_);
    ar & make_nvp("particleID", particleID_);
    ar & make_nvp("particleMajorID", particleMajorID_);
    ar & make_nvp("cherenkovDist", cherenkovDist_);
    ar & make_nvp("wavelength", wavelength_);
    ar & make_nvp("groupVelocity", groupVelocity_);
    ar & make_nvp("dir", direction_);
    ar & make_nvp("pos", position_);
    ar & make_nvp("startTime", startTime_);
    ar & make_nvp("startDir", startDirection_);
    ar & make_nvp("startPos", startPosition_);
    ar & make_nvp("numScattered", numScattered_);
}     

I3_SERIALIZABLE(I3Photon);
I3_SERIALIZABLE(I3PhotonSeriesMap);
