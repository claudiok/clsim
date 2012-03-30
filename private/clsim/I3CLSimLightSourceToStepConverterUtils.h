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
 * @file I3CLSimLightSourceToStepConverterUtils.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMLIGHTSOURCETOSTEPCONVERTERUTILS_H_INCLUDED
#define I3CLSIMLIGHTSOURCETOSTEPCONVERTERUTILS_H_INCLUDED

#include "dataclasses/physics/I3Particle.h"
#include "phys-services/I3RandomService.h"
#include "dataclasses/I3Constants.h"

#include "clsim/I3CLSimStep.h"
#include "clsim/function/I3CLSimFunction.h"

#include <cmath>

namespace I3CLSimLightSourceToStepConverterUtils 
{
    
    double NumberOfPhotonsPerMeter(const I3CLSimFunction &phaseRefIndex,
                                   const I3CLSimFunction &wavelengthGenerationBias,
                                   double fromWlen, double toWlen);
    
    double PhotonNumberCorrectionFactorAfterBias(const I3CLSimFunction &unbiasedSpectrum,
                                                 const I3CLSimFunction &wavelengthGenerationBias,
                                                 double fromWlen, double toWlen);
    
    inline uint64_t mwcRngInitState(I3RandomServicePtr randomService, uint32_t a)
    {
        uint64_t x=0;
        while( (x==0) | (((uint32_t)(x>>32))>=(a-1)) | (((uint32_t)x)>=0xfffffffful))
        {
            // generate random numbers for x and c (both are stored in "x")
            x = static_cast<uint32_t>(randomService->Integer(0xffffffff));
            x=x<<32;
            x += static_cast<uint32_t>(randomService->Integer(0xffffffff));
        }
        return x;
    }
    
    inline double mwcRngRandomNumber_co(uint64_t &state, uint32_t a)
    {
        state=(state&0xfffffffful)*a+(state>>32);
        return static_cast<double>(((uint32_t)(state&0xfffffffful)))/(double)0x100000000;
    }

    inline double mwcRngRandomNumber_oc(uint64_t &state, uint32_t a)
    {
        return 1.0-mwcRngRandomNumber_co(state,a);
    } 

    
    
    
    // stolen from PPC by D. Chirkin
    inline double gammaDistributedNumber(double shape, uint64_t &rngState, uint32_t rngA)
    {
        double x;
        if(shape<1.){  // Weibull algorithm
            double c=1./shape;
            double d=(1.-shape)*std::pow(shape, shape / (1.-shape) );
            double z, e;
            do
            {
                z=-std::log(mwcRngRandomNumber_oc(rngState, rngA));
                e=-std::log(mwcRngRandomNumber_oc(rngState, rngA));
                x=std::pow(z, c);
            } while(z+e<d+x); // or here
        }
        else  // Cheng's algorithm
        {
            double b=shape-std::log(4.0);
            double l=std::sqrt(2.*shape-1.0);
            const double cheng=1.0+std::log(4.5);
            
            //float u, v;
            float y, z, r;
            do
            {
                const double rx = mwcRngRandomNumber_oc(rngState, rngA);
                const double ry = mwcRngRandomNumber_oc(rngState, rngA);
                
                y=log( ry/(1.-ry) ) / l;
                x=shape*std::exp(y);
                z=rx*ry*ry;
                r=b+(shape+l)*y-x;
            } while(r<4.5*z-cheng && r<std::log(z));
        }
        
        return x;
    }

    // stolen from PPC by D. Chirkin
    inline double gammaDistributedNumber(double shape, I3RandomService &randomService)
    {
        double x;
        if(shape<1.){  // Weibull algorithm
            double c=1./shape;
            double d=(1.-shape)*std::pow(shape, shape / (1.-shape) );
            double z, e;
            do
            {
                z=-std::log(randomService.Uniform());
                e=-std::log(randomService.Uniform());
                x=std::pow(z, c);
            } while(z+e<d+x); // or here
        }
        else  // Cheng's algorithm
        {
            double b=shape-std::log(4.0);
            double l=std::sqrt(2.*shape-1.0);
            const double cheng=1.0+std::log(4.5);
            
            //float u, v;
            float y, z, r;
            do
            {
                const double rx = randomService.Uniform();
                const double ry = randomService.Uniform();
                
                y=log( ry/(1.-ry) ) / l;
                x=shape*std::exp(y);
                z=rx*ry*ry;
                r=b+(shape+l)*y-x;
            } while(r<4.5*z-cheng && r<std::log(z));
        }
        
        return x;
    }
    
    // smart pointer version
    inline double gammaDistributedNumber(double shape, I3RandomServicePtr randomService)
    {
        return gammaDistributedNumber(shape, *randomService);
    }
    
    // in-place rotation
    inline void scatterDirectionByAngle(double cosa, double sina,
                                        double &x, double &y, double &z,
                                        double randomValue)
    {
        // randomize direction of scattering (rotation around old direction axis)
        const double b=2.0*M_PI*randomValue;
        const double cosb=std::cos(b);
        const double sinb=std::sin(b);
        
        // Rotate new direction into absolute frame of reference 
        const double sinth = std::sqrt(std::max(0., 1.-z*z));
        
        if(sinth>0.){  // Current direction not vertical, so rotate 
            const double old_x=x;
            const double old_y=y;
            const double old_z=z;
            
            x=old_x*cosa-(old_y*cosb+old_z*old_x*sinb)*sina/sinth;
            y=old_y*cosa+(old_x*cosb-old_z*old_y*sinb)*sina/sinth;
            z=old_z*cosa+sina*sinb*sinth;
        }else{         // Current direction is vertical, so this is trivial
            x=sina*cosb;
            y=sina*sinb;
            if (z>=0.) {
                z=cosa;
            } else {
                z=-cosa;
            }
        }
        
        {
            const double recip_length = 1./std::sqrt( x*x + y*y + z*z );
            
            x *= recip_length;
            y *= recip_length;
            z *= recip_length;
        }
    }
    
    
}


#endif //I3CLSIMLIGHTSOURCETOSTEPCONVERTERUTILS_H_INCLUDED
