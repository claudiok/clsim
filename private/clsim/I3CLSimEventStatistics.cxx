/**
 * copyright (C) 2011
 * Claudio Kopper <claudio.kopper@nikhef.nl>
 * $Id$
 *
 * @file I3CLSimEventStatistics.cxx
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
#include <clsim/I3CLSimEventStatistics.h>

I3CLSimEventStatistics::~I3CLSimEventStatistics() { }

I3CLSimEventStatistics::I3CLSimEventStatistics()
{
    Reset();
}

template <class Archive>
void I3CLSimEventStatistics::serialize (Archive &ar, unsigned version)
{
	if (version > i3clsimeventstatistics_version_)
		log_fatal("Attempting to read version %u from file but running version %u of I3CLSimEventStatistics class.",version,i3clsimeventstatistics_version_);

    ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    
    ar & make_nvp("numberOfPhotonsGeneratedPerParticle",numberOfPhotonsGeneratedPerParticle_);
    ar & make_nvp("sumOfWeightsPhotonsGeneratedPerParticle",sumOfWeightsPhotonsGeneratedPerParticle_);
    ar & make_nvp("totalNumberOfPhotonsGenerated",totalNumberOfPhotonsGenerated_);
    ar & make_nvp("totalSumOfWeightsPhotonsGenerated",totalSumOfWeightsPhotonsGenerated_);

    ar & make_nvp("numberOfPhotonsAtDOMsPerParticle",numberOfPhotonsAtDOMsPerParticle_);
    ar & make_nvp("sumOfWeightsPhotonsAtDOMsPerParticle",sumOfWeightsPhotonsAtDOMsPerParticle_);
    ar & make_nvp("totalNumberOfPhotonsAtDOMs",totalNumberOfPhotonsAtDOMs_);
    ar & make_nvp("totalSumOfWeightsPhotonsAtDOMs",totalSumOfWeightsPhotonsAtDOMs_);

    
}     

I3_SERIALIZABLE(I3CLSimEventStatistics);
