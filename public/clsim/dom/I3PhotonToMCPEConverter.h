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
 * @file I3PhotonToMCPEConverter.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3PHOTONTOMCPECONVERTER_H_INCLUDED
#define I3PHOTONTOMCPECONVERTER_H_INCLUDED

#include "icetray/I3Module.h"
#include "icetray/I3ConditionalModule.h"

#include "dataclasses/geometry/I3Geometry.h"
#include "dataclasses/calibration/I3Calibration.h"
#include "dataclasses/status/I3DetectorStatus.h"

#include "phys-services/I3RandomService.h"

#include "clsim/function/I3CLSimFunction.h"

#include <string>


/**
 * @brief This module reads I3PhotonSeriesMaps generated
 * by CLSim, applies (D)OM acceptances (wavelength&angular)
 * to the photons and stores the results in an I3MCPESeriesMap.
 *
 */
class I3PhotonToMCPEConverter : public I3ConditionalModule
{
public:
    /**
     * Builds an instance of this class
     */
    I3PhotonToMCPEConverter(const I3Context& ctx);
    
    /**
     * Destroys an instance of this class
     */
    virtual ~I3PhotonToMCPEConverter();
    
    /**
     * This module takes a configuration parameter and so it must be configured.
     */
    virtual void Configure();
    
    /**
     * The module needs to process Physics frames
     */
    void DAQ(I3FramePtr frame);

    /**
     * The module needs to process DetectorStatus frames
     */
    virtual void DetectorStatus(I3FramePtr frame);
    
    /**
     * The module needs to process Calibration frames
     */
    virtual void Calibration(I3FramePtr frame);

    virtual void Finish();

    
private:
    // parameters
    
    /// Parameter: A random number generating service (derived from I3RandomService).
    I3RandomServicePtr randomService_;

    /// Parameter: Name of the input I3PhotonSeriesMap frame object. 
    std::string inputPhotonSeriesMapName_;

    /// Parameter: Name of the output I3MCPESeries frame object. 
    std::string outputMCPESeriesMapName_;

    /// Parameter: Name of the I3MCTree frame object. All photon particle IDs are checked against this tree.
    std::string MCTreeName_;

    /// Parameter: Wavelength acceptance of the (D)OM as a I3WlenDependedValue object.
    I3CLSimFunctionConstPtr wavelengthAcceptance_;

    /// Parameter: Angular acceptance of the (D)OM as a I3WlenDependedValue object.
    I3CLSimFunctionConstPtr angularAcceptance_;

    /// Parameter: Specifiy the \"oversize factor\" (i.e. DOM radius scaling factor) you used during the
    ///            CLSim run. The photon arrival times will be corrected.
    double DOMOversizeFactor_;

    /// Parameter: Specifiy the \"pancake factor\" of a DOM. This is the factor a DOM has been *shrunk* again
    ///            (in the direction of the photon) after oversizing. You should set this to whatever
    ///            value you used during running I3CLSimModule. And most of the time this is the same as the
    ///            oversize factor.
    double DOMPancakeFactor_;

    /// Parameter: Specifiy the DOM radius. Do not include oversize factors here.
    double DOMRadiusWithoutOversize_;

    /// Parameter: Default relative efficiency. This value is used if no entry is available from I3Calibration.
    double defaultRelativeDOMEfficiency_;

    /// Parameter: Always use the default relative efficiency, ignore other values from I3Calibration.
    bool replaceRelativeDOMEfficiencyWithDefault_;

    /// Parameter: Do not generate hits for OMKeys not found in the I3DetectorStatus.I3DOMStatusMap
    bool ignoreDOMsWithoutDetectorStatusEntry_;

    /// Parameter: Make photon position/radius check a warning only (instead of a fatal condition)
    bool onlyWarnAboutInvalidPhotonPositions_;

    
private:
    // default, assignment, and copy constructor declared private
    I3PhotonToMCPEConverter();
    I3PhotonToMCPEConverter(const I3PhotonToMCPEConverter&);
    I3PhotonToMCPEConverter& operator=(const I3PhotonToMCPEConverter&);

    I3CalibrationConstPtr calibration_;
    I3DetectorStatusConstPtr status_;
    
    // record some statistics
    uint64_t numGeneratedHits_;
    
    SET_LOGGER("I3PhotonToMCPEConverter");
};

#endif //I3PHOTONTOMCPECONVERTER_H_INCLUDED
