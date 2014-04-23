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
 * @file I3TauSanitizer.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>

#include "icetray/I3Module.h"
#include "icetray/I3ConditionalModule.h"

#include <string>

#include <boost/foreach.hpp>

#include "dataclasses/physics/I3Particle.h"
#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3MCTreeUtils.h"
#include "dataclasses/I3TreeUtils.h"
#include "dataclasses/I3Constants.h"

#include "simclasses/I3MMCTrack.h"

/**
 * Sets taus with length==0 or NaN to shape "Dark".
 */
class I3TauSanitizer : public I3ConditionalModule
{
public:
    /**
     * Builds an instance of this class
     */
    I3TauSanitizer(const I3Context& ctx);
    
    /**
     * Destroys an instance of this class
     */
    virtual ~I3TauSanitizer();
    
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

    // /// Parameter: Name of the I3MMCTrackList frame object. 
    // std::string MMCTrackListName_;

    /// Parameter: Name of the output I3MCTree frame object. If identical to the
    /// input or empty, the input object will be replaced.
    std::string outputMCTreeName_;
    
    
private:
    // default, assignment, and copy constructor declared private
    I3TauSanitizer();
    I3TauSanitizer(const I3TauSanitizer&);
    I3TauSanitizer& operator=(const I3TauSanitizer&);
    
    SET_LOGGER("I3TauSanitizer");
};



// The module
I3_MODULE(I3TauSanitizer);

I3TauSanitizer::I3TauSanitizer(const I3Context& context) 
: I3ConditionalModule(context)
{
    
    inputMCTreeName_="I3MCTree";
    AddParameter("InputMCTreeName",
                 "Name of the I3MCTree frame object.",
                 inputMCTreeName_);

    // inputMCTreeName_="MMCTrackList";
    // AddParameter("MMCTrackListName",
    //              "Name of the I3MMCTrackList frame object.",
    //              MMCTrackListName_);

    inputMCTreeName_="I3MCTree_sanitized";
    AddParameter("OutputMCTreeName",
                 "Name of the output I3MCTree frame object. If identical to the\n"
                 "input or empty, the input object will be replaced.",
                 outputMCTreeName_);


    // add an outbox
    AddOutBox("OutBox");

}

I3TauSanitizer::~I3TauSanitizer()
{
    log_trace("%s", __PRETTY_FUNCTION__);
}


void I3TauSanitizer::Configure()
{
    log_trace("%s", __PRETTY_FUNCTION__);

    GetParameter("InputMCTreeName", inputMCTreeName_);
    // GetParameter("MMCTrackListName", MMCTrackListName_);
    GetParameter("OutputMCTreeName", outputMCTreeName_);

    if (inputMCTreeName_=="")
        log_fatal("The \"InputMCTreeName\" parameter must not be empty.");
}


void I3TauSanitizer::DAQ(I3FramePtr frame)
{
    log_trace("%s", __PRETTY_FUNCTION__);
    
    I3MCTreeConstPtr inputMCTree = frame->Get<I3MCTreeConstPtr>(inputMCTreeName_);
    if (!inputMCTree) {
        log_debug("Frame does not contain an I3MCTree named \"%s\".",
                  inputMCTreeName_.c_str());
        PushFrame(frame);
        return;
    }

    // I3MMCTrackListConstPtr MMCTrackList = frame->Get<I3MMCTrackListConstPtr>(MMCTrackListName_);
    // if (!MMCTrackList) log_fatal("Frame does not contain an I3MMCTrackList named \"%s\".",
    //                              MMCTrackListName_.c_str());

    // allocate the output I3MCTree
    I3MCTreePtr outputMCTree(new I3MCTree(*inputMCTree));

    // set all taus without lengths to "Dark"
    for (I3MCTree::iterator it = outputMCTree->begin();
         it != outputMCTree->end();
         ++it)
    {
        if ((it->GetType() != I3Particle::TauPlus) && (it->GetType() != I3Particle::TauMinus))
            continue;

        if (isnan(it->GetLength()))
        {
            log_warn("Particle (%" PRIu64 ",%i) has NaN length. setting to shape \"Dark\".",
                it->GetMajorID(),
                it->GetMinorID());
            it->SetShape(I3Particle::Dark);
        }
        else if (it->GetLength() <= 0.)
        {
            log_warn("Particle (%" PRIu64 ",%i) has length %fm. setting to shape \"Dark\".",
                it->GetMajorID(),
                it->GetMinorID(),
                it->GetLength()/I3Units::m);
            it->SetShape(I3Particle::Dark);
        }
    }

    // store the output I3MCTree
    if ((outputMCTreeName_=="") || (outputMCTreeName_==inputMCTreeName_)) {
        frame->Delete(inputMCTreeName_);
        frame->Put(inputMCTreeName_, outputMCTree);
    } else {
        frame->Put(outputMCTreeName_, outputMCTree);
    }
    
    // that's it!
    PushFrame(frame);
}
