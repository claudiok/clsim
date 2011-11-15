#include "clsim/I3CLSimParticleToStepConverterUtils.h"

#include "clsim/I3CLSimWlenDependentValue.h"
#include "icetray/I3Units.h"
#include "dataclasses/I3Constants.h"

#include <gsl/gsl_integration.h>

#ifdef HAS_ACCELERATE_FRAMEWORK
#include <Accelerate/Accelerate.h>
#endif

namespace I3CLSimParticleToStepConverterUtils {
#define H_TIMES_C 1.
    
    struct f_params_t {
        const I3CLSimWlenDependentValue *phaseRefIndex;
        const I3CLSimWlenDependentValue *wavelengthGenerationBias;
    };
    
    static double f(double energy, void *params)
    {
        f_params_t *f_params = static_cast<f_params_t *>(params);
        
        const I3CLSimWlenDependentValue &phaseRefIndex = *(f_params->phaseRefIndex);
        const I3CLSimWlenDependentValue &wavelengthGenerationBias = *(f_params->wavelengthGenerationBias);
        
        const double beta=1.;
        const double wlen = H_TIMES_C/energy;
        
        const double refIndexValue = phaseRefIndex.GetValue(wlen);
        const double bias = wavelengthGenerationBias.GetValue(wlen);
        
        const double retval = bias*(2.*M_PI/(137.*H_TIMES_C))*(1. - 1./( pow(beta*refIndexValue, 2.) ));
        
        
        log_trace("r(%fnm)=%g  bias(%fnm)=%g  f(%fnm)=%g",
                  wlen/I3Units::nanometer,
                  refIndexValue,
                  wlen/I3Units::nanometer,
                  bias,
                  wlen/I3Units::nanometer,
                  retval);
        
        return retval;
    }
    
    double NumberOfPhotonsPerMeter(const I3CLSimWlenDependentValue &phaseRefIndex,
                                   const I3CLSimWlenDependentValue &wavelengthGenerationBias,
                                   double fromWlen, double toWlen)
    {
        log_trace("integrating in range [%f;%f]nm",
                  fromWlen/I3Units::nanometer,
                  toWlen/I3Units::nanometer);
        
        f_params_t f_params;
        f_params.phaseRefIndex = &phaseRefIndex;
        f_params.wavelengthGenerationBias = &wavelengthGenerationBias;
        
        gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
        
        gsl_function F;
        F.function = &f;
        F.params = static_cast<void *>(&f_params); // I'm not proud of this..
        
        double result, error;
        
        gsl_integration_qag(&F,
                            H_TIMES_C/toWlen,
                            H_TIMES_C/fromWlen,
                            0,
                            1e-5,
                            1000,
                            GSL_INTEG_GAUSS61,
                            w,
                            &result,
                            &error); 
        
        gsl_integration_workspace_free(w);
        
        return result/(1./I3Units::meter);
    }
#undef H_TIMES_C    
    

    GenerateStepPreCalculator::GenerateStepPreCalculator(I3RandomServicePtr randomService,
                                                         float angularDist_a,
                                                         float angularDist_b,
                                                         std::size_t numberOfValues)
    :
    randomService_(randomService),
    angularDist_a_(angularDist_a),
    one_over_angularDist_a_(numberOfValues, static_cast<float>(1./angularDist_a)),
    angularDist_b_(angularDist_b),
    angularDist_I_(static_cast<float>( 1.-std::exp(-static_cast<double>(angularDist_b)*std::pow(2., static_cast<double>(angularDist_a))) ) ),
    angular_sin_cache_(numberOfValues, NAN),
    angular_cos_cache_(numberOfValues, NAN),
    randomNumber_workspace_(numberOfValues, NAN),
    scratch_space1_(numberOfValues, NAN),
    scratch_space2_(numberOfValues, NAN),
    numberOfValues_(numberOfValues),
    index_(numberOfValues)
    {
        
    }
    
    GenerateStepPreCalculator::~GenerateStepPreCalculator()
    {
    }

    void GenerateStepPreCalculator::RegenerateValues()
    {
        // fill random number array
        for (std::size_t i=0;i<numberOfValues_;++i)
        {
            randomNumber_workspace_[i] = randomService_->Uniform();
        }
        
        // calculate values
#ifdef HAS_ACCELERATE_FRAMEWORK
        const int numVals = numberOfValues_;
        
        for (int i=0;i<numVals;++i)
        {
            scratch_space1_[i] = 1.f-randomNumber_workspace_[i]*angularDist_I_;
        }
        
        vvlogf(&(scratch_space2_[0]), &(scratch_space1_[0]), &numVals);
        
        for (int i=0;i<numVals;++i)
        {
            scratch_space1_[i] = -scratch_space2_[i]/angularDist_b_;
        }
        
        vvpowf(&(scratch_space2_[0]), &(one_over_angularDist_a_[0]), &(scratch_space1_[0]), &numVals);
        
        for (int i=0;i<numVals;++i)
        {
            angular_cos_cache_[i] = std::max(1.f-scratch_space2_[i], -1.0f);
            scratch_space1_[i] = 1.f-angular_cos_cache_[i]*angular_cos_cache_[i];
        }
        
        vvsqrtf(&(angular_sin_cache_[0]), &(scratch_space1_[0]), &numVals);
        
#else
        for (std::size_t i=0;i<numberOfValues_;++i)
        {
            angular_cos_cache_[i]=std::max(1.f-std::pow(-std::log(1.f-randomNumber_workspace_[i]*angularDist_I_)/angularDist_b_, one_over_angularDist_a_[i]), -1.0f);
            angular_sin_cache_[i]=std::sqrt(1.f-angular_cos_cache_[i]*angular_cos_cache_[i]);
        }
#endif
        
        index_=0;
    }


}
