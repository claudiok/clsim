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
 * @file I3MuonSlicer.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
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
    virtual void DAQ(I3FramePtr frame);

    
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
