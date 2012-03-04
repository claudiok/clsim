/**
 * copyright (C) 2011
 * Claudio Kopper <claudio.kopper@nikhef.nl>
 * $Id$
 *
 * @file I3PhotonToMCHitConverter.h
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

#ifndef I3PHOTONTOMCHITCONVERTER_H_INCLUDED
#define I3PHOTONTOMCHITCONVERTER_H_INCLUDED

#include "icetray/I3Module.h"
#include "icetray/I3ConditionalModule.h"

#include "dataclasses/geometry/I3Geometry.h"
#include "dataclasses/calibration/I3Calibration.h"
#include "dataclasses/status/I3DetectorStatus.h"

#include "phys-services/I3RandomService.h"

#include "clsim/I3CLSimWlenDependentValue.h"
#include "clsim/I3CLSimPMTPhotonSimulator.h"

#include <string>


/**
 * @brief This module reads I3PhotonSeriesMaps generated
 * by CLSim, applies (D)OM acceptances (wavelength&angular)
 * to the photons and stores the results in an I3MCHitSeriesMap.
 * After-pulsing is not taken into account, please
 * use the I3AfterPulseSimulatorIceCube on the resulting
 * I3MCHitSeriesMap if you want to have it simulated.
 *
 */
class I3PhotonToMCHitConverter : public I3ConditionalModule
{
public:
    /**
     * Builds an instance of this class
     */
    I3PhotonToMCHitConverter(const I3Context& ctx);
    
    /**
     * Destroys an instance of this class
     */
    ~I3PhotonToMCHitConverter();
    
    /**
     * This module takes a configuration parameter and so it must be configured.
     */
    virtual void Configure();
    
    /**
     * The module needs to process Physics frames
     */
#ifdef IS_Q_FRAME_ENABLED
    void DAQ(I3FramePtr frame);
#else
    void Physics(I3FramePtr frame);
#endif

    /**
     * The module needs to process DetectorStatus frames
     */
    virtual void DetectorStatus(I3FramePtr frame);
    
    /**
     * The module needs to process Calibration frames
     */
    virtual void Calibration(I3FramePtr frame);

    
private:
    // parameters
    
    /// Parameter: A random number generating service (derived from I3RandomService).
    I3RandomServicePtr randomService_;

    /// Parameter: Name of the input I3PhotonSeriesMap frame object. 
    std::string inputPhotonSeriesMapName_;

    /// Parameter: Name of the output I3MCHitSeries frame object. 
    std::string outputMCHitSeriesMapName_;

    /// Parameter: Name of the I3MCTree frame object. All photon particle IDs are checked against this tree.
    std::string MCTreeName_;

    /// Parameter: Wavelength acceptance of the (D)OM as a I3WlenDependedValue object.
    I3CLSimWlenDependentValueConstPtr wavelengthAcceptance_;

    /// Parameter: Angular acceptance of the (D)OM as a I3WlenDependedValue object.
    I3CLSimWlenDependentValueConstPtr angularAcceptance_;

    /// Parameter: Specifiy the \"oversize factor\" (i.e. DOM radius scaling factor) you used during the
    ///            CLSim run. The photon arrival times will be corrected. In practice this means your
    ///            large spherical DOMs will become ellipsoids.
    double DOMOversizeFactor_;

    /// Parameter: Specifiy the DOM radius. Do not include oversize factors here.
    double DOMRadiusWithoutOversize_;

    /// Parameter: Optional after-pulse, late-pulse and jitter simulator object, an instance of I3CLSimPMTPhotonSimulator
    I3CLSimPMTPhotonSimulatorPtr pmtPhotonSimulator_;

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
    I3PhotonToMCHitConverter();
    I3PhotonToMCHitConverter(const I3PhotonToMCHitConverter&);
    I3PhotonToMCHitConverter& operator=(const I3PhotonToMCHitConverter&);

    I3CalibrationConstPtr calibration_;
    I3DetectorStatusConstPtr status_;
    
    SET_LOGGER("I3PhotonToMCHitConverter");
};

#endif //I3PHOTONTOMCHITCONVERTER_H_INCLUDED
