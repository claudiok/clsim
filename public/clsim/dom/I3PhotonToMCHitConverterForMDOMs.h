/**
 * Copyright (c) 2012
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
 * @file I3PhotonToMCHitConverterForMDOMs.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3PHOTONTOMCHITCONVERTERFORMDOMS_H_INCLUDED
#define I3PHOTONTOMCHITCONVERTERFORMDOMS_H_INCLUDED

#include <string>
#include <vector>

#include <icetray/I3Module.h>
#include <icetray/I3ConditionalModule.h>
#include <icetray/I3TrayHeaders.h>
#include <icetray/I3Logging.h>

#include "phys-services/I3RandomService.h"

#include "clsim/function/I3CLSimFunction.h"

/**
 * This module uses PMT and OM acceptance information from the
 * multiPMT-patched I3Geometry class to convert from an I3PhotonMap
 * to an I3MCHitSeriesMapMap (or I3MCHitSeriesMultiOMMap) which 
 * can be passed to a PMT and readout simulation.
 *
 * Functionality from this module should probably integrated in
 * the I3PhotonToMCHitConverter module..
 *
 */
class I3PhotonToMCHitConverterForMDOMs : public I3ConditionalModule
    {
    public:
        
        /**
         * @brief Constructor
         */
        I3PhotonToMCHitConverterForMDOMs(const I3Context& ctx);
        
        /**
         * @brief Destructor
         */
        virtual ~I3PhotonToMCHitConverterForMDOMs();
        
        /**
         * @brief Configure this module
         */
        void Configure();
        
        /**
         * @brief We subscribe to physics events
         *
         * @param frame The frame
         */
        void DAQ(I3FramePtr frame);
        
        /**
         * @brief To clean up
         */
        void Finish();
        
    private:
        /**
         * Parameter:A random number generating service (derived from I3RandomService).
         */
        I3RandomServicePtr randomService_;

        /**
         * Parameter: Name of the input I3PhotonSeriesMap frame object
         */
        std::string inputPhotonSeriesMapName_;

        /**
         * Parameter: Name of the I3MCHitSeriesMap which will be put into the frame.
         */
        std::string outputMCHitSeriesMapName_;
        
        /**
         * Parameter: Name of the I3MCTree frame object. All photon particle IDs are checked against this tree.
         */
        std::string MCTreeName_;

        /**
         * Parameter: A list of subdetectors to ignore
         */
        std::vector<std::string> ignoreSubdetectors_;
        
        /**
         * Parameter: Wavelength acceptance of a PMT as a I3CLSimFunction object. This should include
         * quantum efficiency and collection efficiency effects.
         */
        I3CLSimFunctionConstPtr pmtWavelengthAcceptance_;
        
        /**
         * Parameter: Angular acceptance of a PMT as a I3CLSimFunction object. This is a correction factor
         * to the wavelength acceptance depending on the angle of incidence of a photon. It is assumed to
         * be a function of cos(angle), with angle==0/cos==1 corresponding to a head-on photon.
         */
        I3CLSimFunctionConstPtr pmtAngularAcceptance_;

       /**
        * Parameter: The thickness of the DOM pressure housing glass.
        */
       double glassThickness_;
       
       double DOMOversizeFactor_;
       double DOMPancakeFactor_;

        /**
         * Parameter: The absorption length of the DOM pressure housing glass.
         */
        I3CLSimFunctionConstPtr glassAbsorptionLength_;

        /**
         * Parameter: The absorption length of the optical gel between the DOM sphere and the PMT.
         */
        I3CLSimFunctionConstPtr gelAbsorptionLength_;

        /**
         * @brief The logger can also be used for this module
         */
        SET_LOGGER("I3PhotonToMCHitConverterForMDOMs");
        
    };

#endif

