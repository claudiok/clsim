#include "clsim/I3CLSimParticleToStepConverterCascadeParameterization.h"

#include "clsim/I3CLSimWlenDependentValue.h"
#include "dataclasses/I3Constants.h"

#include <gsl/gsl_integration.h>

const uint32_t I3CLSimParticleToStepConverterCascadeParameterization::default_photonsPerStep=200;


I3CLSimParticleToStepConverterCascadeParameterization::I3CLSimParticleToStepConverterCascadeParameterization
(I3RandomServicePtr randomService,
 uint32_t photonsPerStep)
:
randomService_(randomService),
initialized_(false),
barrier_is_enqueued_(false),
bunchSizeGranularity_(512),
maxBunchSize_(512000),
photonsPerStep_(photonsPerStep)
{
    if (!randomService_)
        throw I3CLSimParticleToStepConverter_exception("No random services was provided!");
    
    if (photonsPerStep_<=0)
        throw I3CLSimParticleToStepConverter_exception("photonsPerStep may not be <= 0!");
        
}

I3CLSimParticleToStepConverterCascadeParameterization::~I3CLSimParticleToStepConverterCascadeParameterization()
{

}

namespace {
    double f(double wlen, void *params)
    {
        const I3CLSimWlenDependentValue *phaseRefIndex = static_cast<const I3CLSimWlenDependentValue *>(params);
        
        const double beta=1.;
        
        const double retval = (2.*M_PI/(137.*( pow(wlen,2.) )))*(1. - 1./( pow(beta*phaseRefIndex->GetValue(wlen), 2.) ));
        
        
        log_trace("r(%fnm)=%g  f(%fnm)=%g",
                  wlen/I3Units::nanometer,
                  phaseRefIndex->GetValue(wlen),
                  wlen/I3Units::nanometer,
                  retval);
        
        return retval;
    }
    
    double NumberOfPhotonsPerMeter(const I3CLSimWlenDependentValue &phaseRefIndex, double fromWlen, double toWlen)
    {
        log_trace("integrating in range [%f;%f]nm",
                  fromWlen/I3Units::nanometer,
                  toWlen/I3Units::nanometer);
        
        gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
        
        gsl_function F;
        F.function = &f;
        F.params = const_cast<void *>(static_cast<const void *>(&phaseRefIndex)); // I'm not proud of this..
        
        double result, error;
        
        gsl_integration_qag(&F,
                            fromWlen,
                            toWlen,
                            0,
                            1e-7,
                            1000,
                            GSL_INTEG_GAUSS61,
                            w,
                            &result,
                            &error); 
        
        gsl_integration_workspace_free(w);
        
        return result/(1./I3Units::meter);
    }
    
    
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

    inline void GenerateStep(I3CLSimStep &newStep, const I3Particle &p,
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

void I3CLSimParticleToStepConverterCascadeParameterization::Initialize()
{
    if (initialized_)
        throw I3CLSimParticleToStepConverter_exception("I3CLSimParticleToStepConverterCascadeParameterization already initialized!");
    
    if (!mediumProperties_)
        throw I3CLSimParticleToStepConverter_exception("MediumProperties not set!");
    
    if (bunchSizeGranularity_ > maxBunchSize_)
        throw I3CLSimParticleToStepConverter_exception("BunchSizeGranularity must not be greater than MaxBunchSize!");
    
    if (maxBunchSize_%bunchSizeGranularity_ != 0)
        throw I3CLSimParticleToStepConverter_exception("MaxBunchSize is not a multiple of BunchSizeGranularity!");
    
    // make a copy of the medium properties
    {
        I3CLSimMediumPropertiesConstPtr copiedMediumProperties(new I3CLSimMediumProperties(*mediumProperties_));
        mediumProperties_ = copiedMediumProperties;
    }

    // calculate the number of photons per meter (at beta==1) for each layer
    meanPhotonsPerMeterInLayer_.clear();
    for (uint32_t i=0;i<mediumProperties_->GetLayersNum();++i)
    {
        const double nPhot = NumberOfPhotonsPerMeter(*(mediumProperties_->GetPhaseRefractiveIndex(i)),
                                                     mediumProperties_->GetMinWavelength(),
                                                     mediumProperties_->GetMaxWavelength());
        
        meanPhotonsPerMeterInLayer_.push_back(nPhot);
        
        log_debug("layer #%u: number of photons per meter=%f (range=[%f;%f]nm)",
                  static_cast<unsigned int>(i),
                  nPhot,
                  mediumProperties_->GetMinWavelength()/I3Units::nanometer,
                  mediumProperties_->GetMaxWavelength()/I3Units::nanometer
                  );
    }

    initialized_=true;
}


bool I3CLSimParticleToStepConverterCascadeParameterization::IsInitialized() const
{
    return initialized_;
}

void I3CLSimParticleToStepConverterCascadeParameterization::SetBunchSizeGranularity(uint64_t num)
{
    if (initialized_)
        throw I3CLSimParticleToStepConverter_exception("I3CLSimParticleToStepConverterCascadeParameterization already initialized!");
    
    if (num<=0)
        throw I3CLSimParticleToStepConverter_exception("BunchSizeGranularity of 0 is invalid!");
    
    bunchSizeGranularity_=num;
}

void I3CLSimParticleToStepConverterCascadeParameterization::SetMaxBunchSize(uint64_t num)
{
    if (initialized_)
        throw I3CLSimParticleToStepConverter_exception("I3CLSimParticleToStepConverterCascadeParameterization already initialized!");

    if (num<=0)
        throw I3CLSimParticleToStepConverter_exception("MaxBunchSize of 0 is invalid!");

    maxBunchSize_=num;
}

void I3CLSimParticleToStepConverterCascadeParameterization::SetMediumProperties(I3CLSimMediumPropertiesConstPtr mediumProperties)
{
    if (initialized_)
        throw I3CLSimParticleToStepConverter_exception("I3CLSimParticleToStepConverterCascadeParameterization already initialized!");

    mediumProperties_=mediumProperties;
}

void I3CLSimParticleToStepConverterCascadeParameterization::EnqueueParticle(const I3Particle &particle, uint32_t identifier)
{
    if (!initialized_)
        throw I3CLSimParticleToStepConverter_exception("I3CLSimParticleToStepConverterCascadeParameterization is not initialized!");

    if (barrier_is_enqueued_)
        throw I3CLSimParticleToStepConverter_exception("A barrier is enqueued! You must receive all steps before enqueuing a new particle.");

    // allocate new step series if necessary
    if (!currentStepSeries_) currentStepSeries_ = I3CLSimStepSeriesPtr(new I3CLSimStepSeries());
    
    // determine current layer
    uint32_t mediumLayer =static_cast<uint32_t>(std::max(0.,(particle.GetPos().GetZ()-mediumProperties_->GetLayersZStart())/(mediumProperties_->GetLayersHeight())));
    if (mediumLayer >= mediumProperties_->GetLayersNum()) mediumLayer=mediumProperties_->GetLayersNum()-1;
    
    const double density = mediumProperties_->GetMediumDensity();
    const double meanPhotonsPerMeter = meanPhotonsPerMeterInLayer_[mediumLayer];
    //const double muonRestMass=0.105658389*I3Units.GeV;  // muon rest mass [GeV]

    log_trace("density=%fg/cm^3 in layer %u", density/(I3Units::g/I3Units::cm3), static_cast<unsigned int>(mediumLayer));
        
    
    const bool isElectron =
    (particle.GetType()==I3Particle::EMinus) ||
    (particle.GetType()==I3Particle::EPlus) ||
    (particle.GetType()==I3Particle::Brems) ||
    (particle.GetType()==I3Particle::DeltaE) ||
    (particle.GetType()==I3Particle::PairProd) ||
    (particle.GetType()==I3Particle::Gamma);

    const bool isHadron =
    (particle.GetType()==I3Particle::Hadrons) ||
    (particle.GetType()==I3Particle::Neutron) ||
    (particle.GetType()==I3Particle::Pi0) ||
    (particle.GetType()==I3Particle::PiPlus) ||
    (particle.GetType()==I3Particle::PiMinus) ||
    (particle.GetType()==I3Particle::K0_Long) ||
    (particle.GetType()==I3Particle::KPlus) ||
    (particle.GetType()==I3Particle::KMinus) ||
    (particle.GetType()==I3Particle::PPlus) ||
    (particle.GetType()==I3Particle::PMinus) ||
    (particle.GetType()==I3Particle::K0_Short) ||
    (particle.GetType()==I3Particle::NuclInt);

    const bool isMuon =
    (particle.GetType()==I3Particle::MuMinus) ||
    (particle.GetType()==I3Particle::MuPlus);

    const double E = particle.GetEnergy()/I3Units::GeV;
    const double logE = log(E);
    const double Lrad=0.358*(I3Units::g/I3Units::cm3)/density;

    if (isElectron) {
        const double pa=2.03+0.604*logE;
        const double pb=Lrad/0.633;
        
        const double nph=5.21*(0.924*I3Units::g/I3Units::cm3)/density;

        const double meanNumPhotons = meanPhotonsPerMeter*nph*E;
        const uint32_t numPhotons = randomService_->Poisson(meanNumPhotons);
        
        const uint32_t numSteps = numPhotons/photonsPerStep_;
        const uint32_t numPhotonsInLastStep = numPhotons%photonsPerStep_;
        
        
        for (uint32_t i=0; i<numSteps; ++i)
        {
            currentStepSeries_->push_back(I3CLSimStep());
            I3CLSimStep &newStep = currentStepSeries_->back();

            const double longitudinalPos = pb*gammaDistributedNumber(pa, randomService_)*I3Units::m;
            GenerateStep(newStep, particle,
                         identifier,
                         randomService_,
                         photonsPerStep_,
                         longitudinalPos);
        }
        
        if (numPhotonsInLastStep > 0)
        {
            currentStepSeries_->push_back(I3CLSimStep());
            I3CLSimStep &newStep = currentStepSeries_->back();

            const double longitudinalPos = pb*gammaDistributedNumber(pa, randomService_)*I3Units::m;
            GenerateStep(newStep, particle,
                         identifier,
                         randomService_,
                         numPhotonsInLastStep,
                         longitudinalPos);
        }
        
        
        log_trace("Generate %u steps for E=%fGeV. (electron)", static_cast<unsigned int>(numSteps+1), E);
    } else if (isHadron) {
        double pa=1.49+0.359*logE;
        double pb=Lrad/0.772;

        const double em=5.21*(0.924*I3Units::g/I3Units::cm3)/density;
        
        double f=1.0;
        const double E0=0.399;
        const double m=0.130;
        const double f0=0.467;
        const double rms0=0.379;
        const double gamma=1.160;
        
        double e=max(10.0, E);
        double F=1.-pow(e/E0, -m)*(1.-f0);
        double dF=F*rms0*pow(log10(e), -gamma);
        do {f=F+dF*randomService_->Gaus(0.,1.);} while((f<0.) || (1.<f));

        const double nph=f*em;

        const double meanNumPhotons = meanPhotonsPerMeter*nph*E;
        const uint32_t numPhotons = randomService_->Poisson(meanNumPhotons);

        const uint32_t numSteps = numPhotons/photonsPerStep_;
        const uint32_t numPhotonsInLastStep = numPhotons%photonsPerStep_;
        
        
        for (uint32_t i=0; i<numSteps; ++i)
        {
            currentStepSeries_->push_back(I3CLSimStep());
            I3CLSimStep &newStep = currentStepSeries_->back();
            
            const double longitudinalPos = pb*gammaDistributedNumber(pa, randomService_)*I3Units::m;
            GenerateStep(newStep, particle,
                         identifier,
                         randomService_,
                         photonsPerStep_,
                         longitudinalPos);
        }
        
        if (numPhotonsInLastStep > 0)
        {
            currentStepSeries_->push_back(I3CLSimStep());
            I3CLSimStep &newStep = currentStepSeries_->back();

            const double longitudinalPos = pb*gammaDistributedNumber(pa, randomService_)*I3Units::m;
            GenerateStep(newStep, particle,
                         identifier,
                         randomService_,
                         numPhotonsInLastStep,
                         longitudinalPos);
        }        
        
        log_trace("Generate %u steps for E=%fGeV. (hadron)", static_cast<unsigned int>(numSteps+1), E);
    } else if (isMuon) {
        const double length = isnan(particle.GetLength())?(2000.*I3Units::m):(particle.GetLength());
        
        if (isnan(particle.GetLength()))
            log_warn("Muon without length found! Assigned a length of 2000m.");

        const double additionalTrackLengthFromCascades=max(0.0, 0.1720+0.0324*logE);

        // calculation from PPC
        const double extr = 1. + max(0.0, 0.1720+0.0324*logE);
        const double cascadeFraction = 1./extr;
        
        const double meanNumPhotonsTotal = meanPhotonsPerMeter*(length/I3Units::m)*extr;
        
        log_warn("meanNumPhotonsTotal=%f", meanNumPhotonsTotal);
        
        const double meanNumPhotonsFromMuon = meanNumPhotonsTotal*(1.-cascadeFraction);
        const uint32_t numPhotonsFromMuon = randomService_->Poisson(meanNumPhotonsFromMuon);
        
        const double meanNumPhotonsFromCascades = meanNumPhotonsTotal*cascadeFraction;
        const uint32_t numPhotonsFromCascades = randomService_->Poisson(meanNumPhotonsFromCascades);


        // steps from muon

        const uint32_t numStepsFromMuon = numPhotonsFromMuon/photonsPerStep_;
        const uint32_t numPhotonsFromMuonInLastStep = numPhotonsFromMuon%photonsPerStep_;
        
        for (uint32_t i=0; i<(numStepsFromMuon+((numPhotonsFromMuonInLastStep>0)?1:0)); ++i)
        {
            currentStepSeries_->push_back(I3CLSimStep());
            I3CLSimStep &newStep = currentStepSeries_->back();

            GenerateStepForMuon(newStep,
                                particle,
                                identifier,
                                photonsPerStep_,
                                length);
        }
        
        if (numPhotonsFromMuonInLastStep>0)
        {
            currentStepSeries_->push_back(I3CLSimStep());
            I3CLSimStep &newStep = currentStepSeries_->back();

            GenerateStepForMuon(newStep,
                                particle,
                                identifier,
                                numPhotonsFromMuonInLastStep,
                                length);
        }
        
        log_warn("Generate %u steps for E=%fGeV. (muon[muon])", static_cast<unsigned int>((numStepsFromMuon+((numPhotonsFromMuonInLastStep>0)?1:0))), E);
        
        
        // steps from cascade
        
        uint32_t photonsPerStepForCascadeLight=photonsPerStep_/10;
        if (photonsPerStepForCascadeLight==0) photonsPerStepForCascadeLight=1;
        
        const uint32_t numStepsFromCascades = numPhotonsFromCascades/photonsPerStepForCascadeLight;
        const uint32_t numPhotonsFromCascadesInLastStep = numStepsFromCascades%photonsPerStepForCascadeLight;

        for (uint32_t i=0; i<numStepsFromCascades; ++i)
        {
            currentStepSeries_->push_back(I3CLSimStep());
            I3CLSimStep &newStep = currentStepSeries_->back();
            
            const double longitudinalPos = randomService_->Uniform()*length*I3Units::m;
            GenerateStep(newStep, particle,
                         identifier,
                         randomService_,
                         photonsPerStep_,
                         longitudinalPos);
        }
        
        if (numPhotonsFromCascadesInLastStep > 0)
        {
            currentStepSeries_->push_back(I3CLSimStep());
            I3CLSimStep &newStep = currentStepSeries_->back();
            
            const double longitudinalPos = randomService_->Uniform()*length*I3Units::m;
            GenerateStep(newStep, particle,
                         identifier,
                         randomService_,
                         numPhotonsFromCascadesInLastStep,
                         longitudinalPos);
        }        
        
        log_warn("Generate %u steps for E=%fGeV. (muon[cascade])", static_cast<unsigned int>((numStepsFromCascades+((numPhotonsFromCascadesInLastStep>0)?1:0))), E);
        
    } else {
        log_fatal("I3CLSimParticleToStepConverterCascadeParameterization cannot handle a %s.",
                  particle.GetTypeString().c_str());
    }
}

void I3CLSimParticleToStepConverterCascadeParameterization::EnqueueBarrier()
{
    if (!initialized_)
        throw I3CLSimParticleToStepConverter_exception("I3CLSimParticleToStepConverterCascadeParameterization is not initialized!");

    if (barrier_is_enqueued_)
        throw I3CLSimParticleToStepConverter_exception("A barrier is already enqueued!");

    barrier_is_enqueued_=true;
}

bool I3CLSimParticleToStepConverterCascadeParameterization::BarrierActive() const
{
    if (!initialized_)
        throw I3CLSimParticleToStepConverter_exception("I3CLSimParticleToStepConverterCascadeParameterization is not initialized!");

    return barrier_is_enqueued_;
}

bool I3CLSimParticleToStepConverterCascadeParameterization::MoreStepsAvailable() const
{
    if (!initialized_)
        throw I3CLSimParticleToStepConverter_exception("I3CLSimParticleToStepConverterCascadeParameterization is not initialized!");

    if (currentStepSeries_) return true;
    return false;
}

I3CLSimStepSeriesConstPtr I3CLSimParticleToStepConverterCascadeParameterization::GetConversionResultWithBarrierInfo(bool &barrierWasReset, double timeout)
{
    if (!initialized_)
        throw I3CLSimParticleToStepConverter_exception("I3CLSimParticleToStepConverterCascadeParameterization is not initialized!");
    
    barrierWasReset=false;
    
    if (!currentStepSeries_) {
        if (barrier_is_enqueued_) {
            if (barrier_is_enqueued_) barrierWasReset=true;
            barrier_is_enqueued_=false;
            return I3CLSimStepSeriesConstPtr();
        }
        throw I3CLSimParticleToStepConverter_exception("I3CLSimParticleToStepConverterCascadeParameterization: no particle is enqueued!");
    }
    if (barrier_is_enqueued_) barrierWasReset=true;
    barrier_is_enqueued_=false;
    
    I3CLSimStepSeriesConstPtr retVal = currentStepSeries_;
    currentStepSeries_.reset();
    
    return retVal;
}
