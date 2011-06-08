/**
 * copyright (C) 2011
 * Claudio Kopper <claudio.kopper@nikhef.nl>
 * $Id$
 *
 * @file I3PhotonToMCHitConverterForMultiPMT.h
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

#ifndef I3PHOTONTOMCHITCONVERTERFORMULTIPMT_H_INCLUDED
#define I3PHOTONTOMCHITCONVERTERFORMULTIPMT_H_INCLUDED

#include <icetray/I3Module.h>
#include <icetray/I3ConditionalModule.h>
#include <icetray/I3TrayHeaders.h>
#include <icetray/I3Logging.h>

#include "phys-services/I3RandomService.h"

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
class I3PhotonToMCHitConverterForMultiPMT : public I3ConditionalModule
	{
	public:
		
		/**
		 * @brief Constructor
		 */
		I3PhotonToMCHitConverterForMultiPMT(const I3Context& ctx);
		
		/**
		 * @brief Destructor
		 */
		~I3PhotonToMCHitConverterForMultiPMT();
		
		/**
		 * @brief Configure this module
		 */
		void Configure();
		
		/**
		 * @brief We subscribe to physics events
		 *
		 * @param frame The frame
		 */
		void Physics(I3FramePtr frame);
		
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
		 * @brief The logger can also be used for this module
		 */
		SET_LOGGER("I3PhotonToMCHitConverterForMultiPMT");
		
	};

#endif

