#ifndef I3CLSIMPARTICLETOSTEPCONVERTERUTILS_H_INCLUDED
#define I3CLSIMPARTICLETOSTEPCONVERTERUTILS_H_INCLUDED

#include "dataclasses/physics/I3Particle.h"
#include "phys-services/I3RandomService.h"
#include "dataclasses/I3Constants.h"

#include "clsim/I3CLSimStep.h"
#include "clsim/I3CLSimWlenDependentValue.h"

namespace I3CLSimParticleToStepConverterUtils 
{
    
    double NumberOfPhotonsPerMeter(const I3CLSimWlenDependentValue &phaseRefIndex,
                                   const I3CLSimWlenDependentValue &wavelengthGenerationBias,
                                   double fromWlen, double toWlen);
    
    // stolen from PPC by D. Chirkin
    inline double gammaDistributedNumber(float shape, I3RandomServicePtr randomService_)
    {
        double x;
        if(shape<1.){  // Weibull algorithm
            double c=1./shape;
            double d=(1.-shape)*pow(shape, shape / (1.-shape) );
            double z, e;
            do
            {
                z=-log(randomService_->Uniform());
                e=-log(randomService_->Uniform());
                x=pow(z, c);
            } while(z+e<d+x); // or here
        }
        else  // Cheng's algorithm
        {
            double b=shape-log(4.0);
            double l=sqrt(2.*shape-1.0);
            const double cheng=1.0+log(4.5);
            
            //float u, v;
            float y, z, r;
            do
            {
                const double rx = randomService_->Uniform();
                const double ry = randomService_->Uniform();
                
                y=log( ry/(1.-ry) ) / l;
                x=shape*exp(y);
                z=rx*ry*ry;
                r=b+(shape+l)*y-x;
            } while(r<4.5*z-cheng && r<log(z));
        }
        
        return x;
    }
    
    // in-place rotation
    inline void scatterDirectionByAngle(double cosa, double sina,
                                        double &x, double &y, double &z,
                                        I3RandomServicePtr randomService_)
    {
        // randomize direction of scattering (rotation around old direction axis)
        const double b=2.0*M_PI*randomService_->Uniform();
        const double cosb=cos(b);
        const double sinb=sin(b);
        
        // Rotate new direction into absolute frame of reference 
        const double sinth = sqrt(max(0., 1.-z*z));
        
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
            const double recip_length = 1./sqrt( x*x + y*y + z*z );
            
            x *= recip_length;
            y *= recip_length;
            z *= recip_length;
        }
    }
    
    inline void GenerateStep(I3CLSimStep &newStep,
                             const I3Particle &p,
                             uint32_t identifier,
                             I3RandomServicePtr randomService_,
                             uint32_t photonsPerStep,
                             const double &longitudinalPos)
    {
        const double angularDist_a=0.39;
        const double angularDist_b=2.61;
        const double angularDist_I=1.-exp(-angularDist_b*pow(2., angularDist_a));
        
        const double angular_cos=max(1.-pow(-log(1.-randomService_->Uniform()*angularDist_I)/angularDist_b, 1./angularDist_a), -1.0);
        const double angular_sin=sqrt(1.-angular_cos*angular_cos);
        
        double step_dx = p.GetDir().GetX();
        double step_dy = p.GetDir().GetY();
        double step_dz = p.GetDir().GetZ();
        
        // set all values
        newStep.SetPosX(p.GetX() + longitudinalPos*step_dx);
        newStep.SetPosY(p.GetY() + longitudinalPos*step_dy);
        newStep.SetPosZ(p.GetZ() + longitudinalPos*step_dz);
        newStep.SetTime(p.GetTime() + longitudinalPos/I3Constants::c);
        
        newStep.SetLength(1.*I3Units::mm);
        newStep.SetNumPhotons(photonsPerStep);
        newStep.SetWeight(1.);
        newStep.SetBeta(1.);
        newStep.SetID(identifier);
        
        // rotate in-place
        scatterDirectionByAngle(angular_cos, angular_sin,
                                step_dx, step_dy, step_dz,
                                randomService_);
        
        newStep.SetDir(step_dx, step_dy, step_dz);
        
        
    }
    
    inline void GenerateStepForMuon(I3CLSimStep &newStep,
                                    const I3Particle &p,
                                    uint32_t identifier,
                                    uint32_t photonsPerStep,
                                    double length)
    {
        // set all values
        newStep.SetPosX(p.GetX());
        newStep.SetPosY(p.GetY());
        newStep.SetPosZ(p.GetZ());
        newStep.SetDir(p.GetDir().GetX(), p.GetDir().GetY(), p.GetDir().GetZ());
        newStep.SetTime(p.GetTime());
        
        newStep.SetLength(length);
        newStep.SetNumPhotons(photonsPerStep);
        newStep.SetWeight(1.);
        newStep.SetBeta(1.);
        newStep.SetID(identifier);
    }
    

}


#endif //I3CLSIMPARTICLETOSTEPCONVERTERUTILS_H_INCLUDED
