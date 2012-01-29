/**
 * copyright (C) 2011
 * Claudio Kopper <claudio.kopper@nikhef.nl>
 * $Id$
 *
 * @file I3CLSimModule.cxx
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

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>

#include "clsim/I3CLSimModuleHelper.h"

#include "clsim/I3CLSimWlenDependentValueConstant.h"
#include "clsim/I3CLSimWlenDependentValueFromTable.h"
#include "clsim/I3CLSimRandomValueInterpolatedDistribution.h"
#include "clsim/I3CLSimRandomValueWlenCherenkovNoDispersion.h"

#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/variant/get.hpp>

#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3MCTreeUtils.h"


namespace I3CLSimModuleHelper {
    
    namespace {
        double CherenkovYieldDistribution(double wlen, I3CLSimMediumPropertiesPtr mediumProperties, double beta=1.)
        {
            I3CLSimWlenDependentValueConstPtr nPhaseDist =
            mediumProperties->GetPhaseRefractiveIndex(0); // this assumes the refractive index does not change between layers
            
            if (!nPhaseDist->HasNativeImplementation()) 
                log_fatal("The refractive index distribution needs a native implementation to be usable!");
            
            const double nPhase = nPhaseDist->GetValue(wlen);
            
            return (2.*M_PI/(137.*(wlen*wlen)))*(1. - 1./ ( std::pow(beta*nPhase,2.) ) ); // dN/dxdwlen
        }

        // the normalization will not be correct here
        double CherenkovYieldDistributionNoDispersion(double wlen)
        {
            return 1./(wlen*wlen); // dN/dxdwlen
        }

    };
    
    
    I3CLSimRandomValueConstPtr
    makeWavelengthGenerator(I3CLSimWlenDependentValueConstPtr wavelengthGenerationBias,
                            bool generateCherenkovPhotonsWithoutDispersion,
                            I3CLSimMediumPropertiesPtr mediumProperties)
    {
        const double minWlen = mediumProperties->GetMinWavelength();
        const double maxWlen = mediumProperties->GetMaxWavelength();
        const double wlenRange = maxWlen-minWlen;
        if (wlenRange <= 0.) log_fatal("Internal error, wavelength range <= 0!");

        if (wavelengthGenerationBias->GetMinWlen() > minWlen)
            log_fatal("wavelength generation bias has to have a wavelength range larger or equal to the medium property range!");
        if (wavelengthGenerationBias->GetMaxWlen() < maxWlen)
            log_fatal("wavelength generation bias has to have a wavelength range larger or equal to the medium property range!");
        
        bool noBias=false;
        bool biasIsConstant=false;
        
        {
            I3CLSimWlenDependentValueConstantConstPtr wavelengthGenerationBiasConstant =
            dynamic_pointer_cast<const I3CLSimWlenDependentValueConstant>(wavelengthGenerationBias);
            
            if (wavelengthGenerationBiasConstant)
            {
                biasIsConstant=true;
                
                if ( std::abs(wavelengthGenerationBiasConstant->GetValue((minWlen+maxWlen)/2.)-1.) < 1e-10 )
                    noBias=true;
            }
        }
        
        I3CLSimWlenDependentValueFromTableConstPtr wavelengthGenerationBiasFromTable;
        if (!biasIsConstant)
        {
            wavelengthGenerationBiasFromTable =
            dynamic_pointer_cast<const I3CLSimWlenDependentValueFromTable>(wavelengthGenerationBias);
        }
        
        

        if ((!noBias) && (generateCherenkovPhotonsWithoutDispersion))
        {
            log_warn("**********");
            log_warn(" Using the \"GenerateCherenkovPhotonsWithoutDispersion\" option");
            log_warn(" with a biased photon spectrum generation does not yield a performance");
            log_warn(" increase. You might consider turning this option off to get a better");
            log_warn(" approximation of the Cherenkov spectrum.");
            log_warn("**********");
        }
        
        // Check if the bias values are from a tabulated distribution.
        // If yes, use the table binning, if no, make up a binning.
        
        if (wavelengthGenerationBiasFromTable)
        {
            std::size_t wlenPoints = wavelengthGenerationBiasFromTable->GetWavelengthNumValues();
            const double firstWlen = wavelengthGenerationBiasFromTable->GetFirstWavelength();
            const double wlenStep = wavelengthGenerationBiasFromTable->GetWavelengthStepping();
            
            std::vector<double> spectrum(wlenPoints, NAN);
            for (std::size_t i=0;i<wlenPoints;++i)
            {
                const double wavelength = firstWlen + static_cast<double>(i)*wlenStep;
                const double bias = wavelengthGenerationBiasFromTable->GetWavelengthValue(i);
                
                if (generateCherenkovPhotonsWithoutDispersion)
                {
                    spectrum[i] =
                    bias * CherenkovYieldDistributionNoDispersion(wavelength);
                }
                else
                {
                    spectrum[i] =
                    bias * CherenkovYieldDistribution(wavelength, mediumProperties);
                }
            }
            
            return I3CLSimRandomValueInterpolatedDistributionConstPtr
            (new I3CLSimRandomValueInterpolatedDistribution(firstWlen,
                                                            wlenStep,
                                                            spectrum));
        }
        else if ((noBias) && (generateCherenkovPhotonsWithoutDispersion))
        {
            return I3CLSimRandomValueWlenCherenkovNoDispersionConstPtr
            (new I3CLSimRandomValueWlenCherenkovNoDispersion(minWlen, maxWlen));
        }
        else
        {
            std::size_t wlenPoints = static_cast<std::size_t>(wlenRange/(10.*I3Units::nanometer))+2;
            const double firstWlen = minWlen;
            const double wlenStep = wlenRange/static_cast<double>(wlenPoints-1);
            
            std::vector<double> spectrum(wlenPoints, NAN);
            for (std::size_t i=0;i<wlenPoints;++i)
            {
                const double wavelength = firstWlen + static_cast<double>(i)*wlenStep;
                const double bias = wavelengthGenerationBias->GetValue(wavelength);
                
                if (generateCherenkovPhotonsWithoutDispersion)
                {
                    spectrum[i] =
                    bias * CherenkovYieldDistributionNoDispersion(wavelength);
                }
                else
                {
                    spectrum[i] =
                    bias * CherenkovYieldDistribution(wavelength, mediumProperties);
                }
            }

            return I3CLSimRandomValueInterpolatedDistributionConstPtr
            (new I3CLSimRandomValueInterpolatedDistribution(firstWlen,
                                                            wlenStep,
                                                            spectrum));
        }

    
    }

    
    I3CLSimStepToPhotonConverterOpenCLPtr initializeOpenCL(const I3CLSimOpenCLDevice &device,
                                                           I3RandomServicePtr rng,
                                                           I3CLSimSimpleGeometryFromI3GeometryPtr geometry,
                                                           I3CLSimMediumPropertiesPtr medium,
                                                           I3CLSimWlenDependentValueConstPtr wavelengthGenerationBias,
                                                           I3CLSimRandomValueConstPtr wavelengthGenerator)
    {
        I3CLSimStepToPhotonConverterOpenCLPtr conv(new I3CLSimStepToPhotonConverterOpenCL(rng, device.GetUseNativeMath()));

        conv->SetDevice(device);

        conv->SetWlenGenerator(wavelengthGenerator);
        conv->SetWlenBias(wavelengthGenerationBias);

        conv->SetMediumProperties(medium);
        conv->SetGeometry(geometry);
        
        conv->Compile();
        //log_trace("%s", conv.GetFullSource().c_str());
        
        const std::size_t maxWorkgroupSize = conv->GetMaxWorkgroupSize();
        conv->SetWorkgroupSize(maxWorkgroupSize);

        const std::size_t workgroupSize = conv->GetWorkgroupSize();
        
        // use approximately the given number of work items, convert to a multiple of the workgroup size
        std::size_t maxNumWorkitems = (static_cast<std::size_t>(device.GetApproximateNumberOfWorkItems())/workgroupSize)*workgroupSize;
        if (maxNumWorkitems==0) maxNumWorkitems=workgroupSize;
        
        conv->SetMaxNumWorkitems(maxNumWorkitems);

        log_info("maximum workgroup size is %zu", maxWorkgroupSize);
        log_info("configured workgroup size is %zu", workgroupSize);
        if (maxNumWorkitems != device.GetApproximateNumberOfWorkItems()) {
            log_warn("maximum number of work items is %zu (user configured was %" PRIu32 ")", maxNumWorkitems, device.GetApproximateNumberOfWorkItems());
        } else {
            log_info("maximum number of work items is %zu (user configured was %" PRIu32 ")", maxNumWorkitems, device.GetApproximateNumberOfWorkItems());
        }

        conv->Initialize();
        
        return conv;
    }

    I3CLSimLightSourceToStepConverterGeant4Ptr initializeGeant4(I3RandomServicePtr rng,
                                                             I3CLSimMediumPropertiesPtr medium,
                                                             I3CLSimWlenDependentValueConstPtr wavelengthGenerationBias,
                                                             uint64_t bunchSizeGranularity,
                                                             uint64_t maxBunchSize,
                                                             const I3CLSimLightSourceParameterizationSeries &parameterizationList,
                                                             const std::string &physicsListName,
                                                             double maxBetaChangePerStep,
                                                             uint32_t maxNumPhotonsPerStep,
                                                             bool multiprocessor)
    {
        I3CLSimLightSourceToStepConverterGeant4Ptr conv
        (
         new I3CLSimLightSourceToStepConverterGeant4
         (
          physicsListName,
          maxBetaChangePerStep,
          maxNumPhotonsPerStep
         )
        );
        
        conv->SetRandomService(rng);
        conv->SetWlenBias(wavelengthGenerationBias);
        conv->SetMediumProperties(medium);
        conv->SetMaxBunchSize(maxBunchSize);
        conv->SetBunchSizeGranularity(bunchSizeGranularity);
        
        conv->SetParticleParameterizationSeries(parameterizationList);
        
        conv->Initialize();
        
        return conv;
    }

}
