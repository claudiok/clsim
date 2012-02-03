/**
 * copyright (C) 2011
 * Claudio Kopper <claudio.kopper@nikhef.nl>
 * $Id$
 *
 * @file I3CLSimModule.h
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

#ifndef I3CLSIMMODULEHELPER_H_INCLUDED
#define I3CLSIMMODULEHELPER_H_INCLUDED


#include "phys-services/I3RandomService.h"

#include "clsim/I3CLSimRandomValue.h"
#include "clsim/I3CLSimWlenDependentValue.h"

#include "clsim/I3CLSimMediumProperties.h"
#include "clsim/I3CLSimSimpleGeometryFromI3Geometry.h"

#include "clsim/I3CLSimStepToPhotonConverterOpenCL.h"
#include "clsim/I3CLSimLightSourceToStepConverterGeant4.h"

#include "clsim/I3CLSimLightSourceParameterization.h"

#include "clsim/I3CLSimOpenCLDevice.h"

#include <vector>
#include <string>

namespace I3CLSimModuleHelper {
    I3CLSimStepToPhotonConverterOpenCLPtr
    initializeOpenCL(const I3CLSimOpenCLDevice &device,
                     I3RandomServicePtr rng,
                     I3CLSimSimpleGeometryFromI3GeometryPtr geometry,
                     I3CLSimMediumPropertiesPtr medium,
                     I3CLSimWlenDependentValueConstPtr wavelengthGenerationBias,
                     const std::vector<I3CLSimRandomValueConstPtr> &wavelengthGenerators);
    
    I3CLSimLightSourceToStepConverterGeant4Ptr
    initializeGeant4(I3RandomServicePtr rng,
                     I3CLSimMediumPropertiesPtr medium,
                     I3CLSimWlenDependentValueConstPtr wavelengthGenerationBias,
                     uint64_t bunchSizeGranularity,
                     uint64_t maxBunchSize,
                     const I3CLSimLightSourceParameterizationSeries &parameterizationList,
                     const std::string &physicsListName,
                     double maxBetaChangePerStep,
                     uint32_t maxNumPhotonsPerStep,
                     bool multiprocessor=false);

    I3CLSimRandomValueConstPtr
    makeWavelengthGenerator(I3CLSimWlenDependentValueConstPtr wavelengthGenerationBias,
                                     bool generateCherenkovPhotonsWithoutDispersion,
                                     I3CLSimMediumPropertiesPtr mediumProperties);


};


#endif //I3CLSIMMODULEHELPER_H_INCLUDED
