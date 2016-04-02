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
 * @file I3ShadowedPhotonRemoverModule.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3SHADOWEDPHOTONREMOVERMODULE_H_INCLUDED
#define I3SHADOWEDPHOTONREMOVERMODULE_H_INCLUDED

#include "icetray/I3Module.h"
#include "icetray/I3ConditionalModule.h"

#include "dataclasses/geometry/I3Geometry.h"

#include "clsim/shadow/I3ShadowedPhotonRemover.h"

#include <string>


/**
 * @brief This module removes photons that have paths intersecting 
 *   with any shadowing part of the detecor (such as cables).
 *   This code is NOT functional at the moment.
 */
class I3ShadowedPhotonRemoverModule : public I3ConditionalModule
{
public:
    /**
     * Builds an instance of this class
     */
    I3ShadowedPhotonRemoverModule(const I3Context& ctx);
    
    /**
     * Destroys an instance of this class
     */
    ~I3ShadowedPhotonRemoverModule();
    
    /**
     * This module takes a configuration parameter and so it must be configured.
     */
    virtual void Configure();
    
    /**
     * The module needs to process Physics frames
     */
    void DAQ(I3FramePtr frame);

    
private:
    // parameters
    
    /// Parameter: Name of the input I3PhotonSeriesMap frame object. 
    std::string inputPhotonSeriesMapName_;

    /// Parameter: Name of the output I3PhotonSeriesMap frame object. 
    std::string outputPhotonSeriesMapName_;

    I3ShadowedPhotonRemoverPtr shadowedPhotonRemover_;
    
private:
    // default, assignment, and copy constructor declared private
    I3ShadowedPhotonRemoverModule();
    I3ShadowedPhotonRemoverModule(const I3ShadowedPhotonRemoverModule&);
    I3ShadowedPhotonRemoverModule& operator=(const I3ShadowedPhotonRemoverModule&);

    SET_LOGGER("I3ShadowedPhotonRemoverModule");
};

#endif //I3SHADOWEDPHOTONREMOVERMODULE_H_INCLUDED
