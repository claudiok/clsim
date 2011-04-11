#include "clsim/I3CLSimParticleToStepConverterCascadeParameterization.h"

#include "clsim/I3CLSimWlenDependentValue.h"


#include "clsim/I3CLSimParticleToStepConverterUtils.h"
using namespace I3CLSimParticleToStepConverterUtils;


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

void I3CLSimParticleToStepConverterCascadeParameterization::Initialize()
{
    if (initialized_)
        throw I3CLSimParticleToStepConverter_exception("I3CLSimParticleToStepConverterCascadeParameterization already initialized!");

    if (!wlenBias_)
        throw I3CLSimParticleToStepConverter_exception("WlenBias not set!");

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
                                                     *(wlenBias_),
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

void I3CLSimParticleToStepConverterCascadeParameterization::SetWlenBias(I3CLSimWlenDependentValueConstPtr wlenBias)
{
    if (initialized_)
        throw I3CLSimParticleToStepConverter_exception("I3CLSimParticleToStepConverterCascadeParameterization already initialized!");
    
    wlenBias_=wlenBias;
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

        // calculation from PPC
        const double extr = 1. + max(0.0, 0.1720+0.0324*logE);
        const double cascadeFraction = 1./extr;
        
        const double meanNumPhotonsTotal = meanPhotonsPerMeter*(length/I3Units::m)*extr;
        
        log_trace("meanNumPhotonsTotal=%f", meanNumPhotonsTotal);
        
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
        
        log_trace("Generate %u steps for E=%fGeV, l=%fm. (muon[muon])", static_cast<unsigned int>((numStepsFromMuon+((numPhotonsFromMuonInLastStep>0)?1:0))), E, length/I3Units::m);
        
        
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
        
        log_trace("Generate %u steps for E=%fGeV, l=%fm. (muon[cascade])", static_cast<unsigned int>((numStepsFromCascades+((numPhotonsFromCascadesInLastStep>0)?1:0))), E, length/I3Units::m);
        
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
