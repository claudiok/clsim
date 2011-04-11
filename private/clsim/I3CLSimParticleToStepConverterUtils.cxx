#include "clsim/I3CLSimParticleToStepConverterUtils.h"

#include "clsim/I3CLSimWlenDependentValue.h"
#include "icetray/I3Units.h"
#include "dataclasses/I3Constants.h"

#include <gsl/gsl_integration.h>

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
    

}
