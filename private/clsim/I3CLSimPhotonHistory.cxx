/**
 * copyright (C) 2012
 * Claudio Kopper <claudio.kopper@icecube.wisc.edu>
 * $Id$
 *
 * @file I3CLSimPhotonHistory.cxx
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
#include <clsim/I3CLSimPhotonHistory.h>

I3CLSimPhotonHistory::~I3CLSimPhotonHistory() { }


template <class Archive>
void I3CLSimPhotonHistory::serialize(Archive &ar, unsigned version)
{
    if (version > i3clsimphotonhistory_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimPhotonHistory class.",version,i3clsimphotonhistory_version_);
    
    ar & make_nvp("posX", posX_);
    ar & make_nvp("posY", posY_);
    ar & make_nvp("posZ", posZ_);
}     


I3_SERIALIZABLE(I3CLSimPhotonHistory);
I3_SERIALIZABLE(I3CLSimPhotonHistorySeries);
