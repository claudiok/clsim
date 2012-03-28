/**
 * copyright (C) 2011
 * Claudio Kopper <claudio.kopper@nikhef.nl>
 * $Id$
 *
 * @file I3MuonSlicer.h
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

#ifndef I3MUONSLICER_H_INCLUDED
#define I3MUONSLICER_H_INCLUDED

#include "icetray/I3Module.h"
#include "icetray/I3ConditionalModule.h"

#include <string>


/**
 * @brief takes an MC tree with a muon that has been
 * propagated using MMC (i.e. a muon with cascades as daughter
 * particles). It chops the muon track into slices 
 * with estimated energies. All the original muon
 * will be retained and all slices will be added to
 * the muon as daughter particles. So you will end up with
 * a long muon track with cascades and shorter muons
 * as daughter particles. Make sure your light propagation
 * knows how to handle this (you run the risk of having light
 * generated twice, once for the parent (=long) track and
 * again for the daughters muon tracks).
 *
 * I3CLSimModule knows how to deal with these tracks.
 *
 * The muon starting/stopping points and energies are taken
 * from an I3MMCTrackList object in case it exists
 * and a name is configured. 
 *
 */
class I3MuonSlicer : public I3ConditionalModule
{
public:
    /**
     * Builds an instance of this class
     */
    I3MuonSlicer(const I3Context& ctx);
    
    /**
     * Destroys an instance of this class
     */
    virtual ~I3MuonSlicer();
    
    /**
     * This module takes a configuration parameter and so it must be configured.
     */
    virtual void Configure();
    
    /**
     * The module needs to process Physics frames
     */
#ifdef IS_Q_FRAME_ENABLED
    virtual void DAQ(I3FramePtr frame);
#else
    virtual void Physics(I3FramePtr frame);
#endif

    
private:
    // parameters
    
    /// Parameter: Name of the I3MCTree frame object. 
    std::string inputMCTreeName_;

    /// Parameter: Name of the I3MMCTrackList frame object. 
    std::string MMCTrackListName_;

    /// Parameter: Name of the output I3MCTree frame object. If identical to the
    /// input or empty, the input object will be replaced.
    std::string outputMCTreeName_;
    
    
private:
    // default, assignment, and copy constructor declared private
    I3MuonSlicer();
    I3MuonSlicer(const I3MuonSlicer&);
    I3MuonSlicer& operator=(const I3MuonSlicer&);

    
    SET_LOGGER("I3MuonSlicer");
};

#endif //I3MUONSLICER_H_INCLUDED
