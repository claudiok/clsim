#include "clsim/I3CLSimLightSourceToStepConverterUtils.h"

#include "clsim/function/I3CLSimFunction.h"
#include "icetray/I3Units.h"
#include "dataclasses/I3Constants.h"

#include <gsl/gsl_integration.h>

namespace I3CLSimLightSourceToStepConverterUtils {
#define H_TIMES_C 1.
    
    struct f_params_t {
        const I3CLSimFunction *phaseRefIndex;
        const I3CLSimFunction *wavelengthGenerationBias;
    };
    
    static double f(double energy, void *params)
    {
        f_params_t *f_params = static_cast<f_params_t *>(params);
        
        const I3CLSimFunction &phaseRefIndex = *(f_params->phaseRefIndex);
        const I3CLSimFunction &wavelengthGenerationBias = *(f_params->wavelengthGenerationBias);
        
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
    
    double NumberOfPhotonsPerMeter(const I3CLSimFunction &phaseRefIndex,
                                   const I3CLSimFunction &wavelengthGenerationBias,
                                   double fromWlen, double toWlen)
    {
        log_trace("integrating in range [%f;%f]nm",
                  fromWlen/I3Units::nanometer,
                  toWlen/I3Units::nanometer);
        
        f_params_t f_params;
        f_params.phaseRefIndex = &phaseRefIndex;
        f_params.wavelengthGenerationBias = &wavelengthGenerationBias;
        
        gsl_integration_workspace *w = gsl_integration_workspace_alloc(10000);
        
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
    

    
    
    
    
    struct f2_params_t {
        const I3CLSimFunction *unbiasedSpectrum;
        const I3CLSimFunction *wavelengthGenerationBias; // can be NULL
    };
    
    static double f2(double wlen, void *params)
    {
        f2_params_t *f2_params = static_cast<f2_params_t *>(params);
        
        double bias = 1.; // no bias by default
        if (f2_params->wavelengthGenerationBias) {
            const I3CLSimFunction &wavelengthGenerationBias = *(f2_params->wavelengthGenerationBias);
            bias = wavelengthGenerationBias.GetValue(wlen);
        }

        const I3CLSimFunction &unbiasedSpectrum = *(f2_params->unbiasedSpectrum);
        double spectrumValue = unbiasedSpectrum.GetValue(wlen);
        
        const double retval = bias*spectrumValue;
        
        log_trace("spectrumValue(%fnm)=%g  bias(%fnm)=%g  f(%fnm)=%g",
                  wlen/I3Units::nanometer,
                  spectrumValue,
                  wlen/I3Units::nanometer,
                  bias,
                  wlen/I3Units::nanometer,
                  retval);
        
        return retval;
    }
    
    double PhotonNumberCorrectionFactorAfterBias(const I3CLSimFunction &unbiasedSpectrum,
                                                 const I3CLSimFunction &wavelengthGenerationBias,
                                                 double fromWlen, double toWlen)
    {
        log_trace("integrating in range [%f;%f]nm",
                  fromWlen/I3Units::nanometer,
                  toWlen/I3Units::nanometer);

        gsl_integration_workspace *w = gsl_integration_workspace_alloc(10000);

        f2_params_t f2_params;
        gsl_function F;
        F.function = &f2;
        F.params = static_cast<void *>(&f2_params); // I'm not proud of this..

        double result, error;

        f2_params.unbiasedSpectrum = &unbiasedSpectrum;
        
        // without bias first
        f2_params.wavelengthGenerationBias = NULL;
        
        gsl_integration_qag(&F,
                            fromWlen,
                            toWlen,
                            0,
                            1e-5,
                            10000,
                            GSL_INTEG_GAUSS61,
                            w,
                            &result,
                            &error); 

        const double integralUnbiased = result;
        
        // and again, this time with the bias
        f2_params.wavelengthGenerationBias = &wavelengthGenerationBias;

        gsl_integration_qag(&F,
                            fromWlen,
                            toWlen,
                            0,
                            1e-5,
                            10000,
                            GSL_INTEG_GAUSS61,
                            w,
                            &result,
                            &error); 
        const double integralWithBias = result;
        
        gsl_integration_workspace_free(w);
        
        return integralWithBias/integralUnbiased;
    }


}
