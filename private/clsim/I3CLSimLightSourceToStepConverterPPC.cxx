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
 * @file I3CLSimLightSourceToStepConverterPPC.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>

#include <cmath>

#include "clsim/I3CLSimLightSourceToStepConverterPPC.h"
#include "sim-services/I3SimConstants.h"

#include "clsim/function/I3CLSimFunction.h"

#include "clsim/I3CLSimLightSourceToStepConverterUtils.h"
using namespace I3CLSimLightSourceToStepConverterUtils;


const uint32_t I3CLSimLightSourceToStepConverterPPC::default_photonsPerStep=200;
const uint32_t I3CLSimLightSourceToStepConverterPPC::default_highPhotonsPerStep=0;
const double I3CLSimLightSourceToStepConverterPPC::default_useHighPhotonsPerStepStartingFromNumPhotons=1.0e9;



I3CLSimLightSourceToStepConverterPPC::I3CLSimLightSourceToStepConverterPPC
(uint32_t photonsPerStep,
 uint32_t highPhotonsPerStep,
 double useHighPhotonsPerStepStartingFromNumPhotons)
:
initialized_(false),
barrier_is_enqueued_(false),
bunchSizeGranularity_(1),
maxBunchSize_(512000),
photonsPerStep_(photonsPerStep),
highPhotonsPerStep_(highPhotonsPerStep),
useHighPhotonsPerStepStartingFromNumPhotons_(useHighPhotonsPerStepStartingFromNumPhotons),
useCascadeExtension_(true)
{
    if (photonsPerStep_<=0)
        throw I3CLSimLightSourceToStepConverter_exception("photonsPerStep may not be <= 0!");
    
    if (highPhotonsPerStep_==0) highPhotonsPerStep_=photonsPerStep_;
    
    if (highPhotonsPerStep_<photonsPerStep_)
        throw I3CLSimLightSourceToStepConverter_exception("highPhotonsPerStep may not be < photonsPerStep!");
    
    if (useHighPhotonsPerStepStartingFromNumPhotons_<=0.)
        throw I3CLSimLightSourceToStepConverter_exception("useHighPhotonsPerStepStartingFromNumPhotons may not be <= 0!");
}

I3CLSimLightSourceToStepConverterPPC::~I3CLSimLightSourceToStepConverterPPC()
{

}

void I3CLSimLightSourceToStepConverterPPC::Initialize()
{
    if (initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterPPC already initialized!");

    if (!randomService_)
        throw I3CLSimLightSourceToStepConverter_exception("RandomService not set!");

    if (!wlenBias_)
        throw I3CLSimLightSourceToStepConverter_exception("WlenBias not set!");

    if (!mediumProperties_)
        throw I3CLSimLightSourceToStepConverter_exception("MediumProperties not set!");
    
    if (bunchSizeGranularity_ > maxBunchSize_)
        throw I3CLSimLightSourceToStepConverter_exception("BunchSizeGranularity must not be greater than MaxBunchSize!");

    if (maxBunchSize_%bunchSizeGranularity_ != 0)
        throw I3CLSimLightSourceToStepConverter_exception("MaxBunchSize is not a multiple of BunchSizeGranularity!");
    
    // initialize a fast rng
    rngA_ = 1640531364; // magic number from numerical recipies
    rngState_ = mwcRngInitState(randomService_, rngA_);
    
    // initialize the pre-calculator threads
    preCalc_ = boost::shared_ptr<GenerateStepPreCalculator>(new GenerateStepPreCalculator(randomService_, /*a=*/0.39, /*b=*/2.61));

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


bool I3CLSimLightSourceToStepConverterPPC::IsInitialized() const
{
    return initialized_;
}

void I3CLSimLightSourceToStepConverterPPC::SetBunchSizeGranularity(uint64_t num)
{
    if (initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterPPC already initialized!");
    
    if (num<=0)
        throw I3CLSimLightSourceToStepConverter_exception("BunchSizeGranularity of 0 is invalid!");

    if (num!=1)
        throw I3CLSimLightSourceToStepConverter_exception("A BunchSizeGranularity != 1 is currently not supported!");

    bunchSizeGranularity_=num;
}

void I3CLSimLightSourceToStepConverterPPC::SetMaxBunchSize(uint64_t num)
{
    if (initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterPPC already initialized!");

    if (num<=0)
        throw I3CLSimLightSourceToStepConverter_exception("MaxBunchSize of 0 is invalid!");

    maxBunchSize_=num;
}

void I3CLSimLightSourceToStepConverterPPC::SetWlenBias(I3CLSimFunctionConstPtr wlenBias)
{
    if (initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterPPC already initialized!");
    
    wlenBias_=wlenBias;
}

void I3CLSimLightSourceToStepConverterPPC::SetRandomService(I3RandomServicePtr random)
{
    if (initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterPPC already initialized!");
    
    randomService_=random;
}

void I3CLSimLightSourceToStepConverterPPC::SetMediumProperties(I3CLSimMediumPropertiesConstPtr mediumProperties)
{
    if (initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterPPC already initialized!");

    mediumProperties_=mediumProperties;
}

void I3CLSimLightSourceToStepConverterPPC::EnqueueLightSource(const I3CLSimLightSource &lightSource, uint32_t identifier)
{
    if (!initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterPPC is not initialized!");

    if (barrier_is_enqueued_)
        throw I3CLSimLightSourceToStepConverter_exception("A barrier is enqueued! You must receive all steps before enqueuing a new particle.");

    if (lightSource.GetType() != I3CLSimLightSource::Particle)
        throw I3CLSimLightSourceToStepConverter_exception("The I3CLSimLightSourceToStepConverterPPC parameterization only works on particles.");
    
    const I3Particle &particle = lightSource.GetParticle();
    
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
    (particle.GetType()==I3Particle::Gamma) ||
    (particle.GetType()==I3Particle::Pi0) ; // Pi0 decays to 2 gammas and produce EM showers

    bool isHadron =
    (particle.GetType()==I3Particle::Hadrons) ||
    (particle.GetType()==I3Particle::Neutron) ||
    (particle.GetType()==I3Particle::PiPlus) ||
    (particle.GetType()==I3Particle::PiMinus) ||
    (particle.GetType()==I3Particle::K0_Long) ||
    (particle.GetType()==I3Particle::KPlus) ||
    (particle.GetType()==I3Particle::KMinus) ||
    (particle.GetType()==I3Particle::PPlus) ||
    (particle.GetType()==I3Particle::PMinus) ||
    (particle.GetType()==I3Particle::K0_Short) ||
#ifdef I3PARTICLE_SUPPORTS_PDG_ENCODINGS
    (particle.GetType()==I3Particle::Eta) ||
    (particle.GetType()==I3Particle::Lambda) ||
    (particle.GetType()==I3Particle::SigmaPlus) ||
    (particle.GetType()==I3Particle::Sigma0) ||
    (particle.GetType()==I3Particle::SigmaMinus) ||
    (particle.GetType()==I3Particle::Xi0) ||
    (particle.GetType()==I3Particle::XiMinus) ||
    (particle.GetType()==I3Particle::OmegaMinus) ||
    (particle.GetType()==I3Particle::NeutronBar) ||
    (particle.GetType()==I3Particle::LambdaBar) ||
    (particle.GetType()==I3Particle::SigmaMinusBar) ||
    (particle.GetType()==I3Particle::Sigma0Bar) ||
    (particle.GetType()==I3Particle::SigmaPlusBar) ||
    (particle.GetType()==I3Particle::Xi0Bar) ||
    (particle.GetType()==I3Particle::XiPlusBar) ||
    (particle.GetType()==I3Particle::OmegaPlusBar) ||
    (particle.GetType()==I3Particle::DPlus) ||
    (particle.GetType()==I3Particle::DMinus) ||
    (particle.GetType()==I3Particle::D0) ||
    (particle.GetType()==I3Particle::D0Bar) ||
    (particle.GetType()==I3Particle::DsPlus) ||
    (particle.GetType()==I3Particle::DsMinusBar) ||
    (particle.GetType()==I3Particle::LambdacPlus) ||
    (particle.GetType()==I3Particle::WPlus) ||
    (particle.GetType()==I3Particle::WMinus) ||
    (particle.GetType()==I3Particle::Z0) ||
#endif
    (particle.GetType()==I3Particle::NuclInt);

    const bool isMuon =
    (particle.GetType()==I3Particle::MuMinus) ||
    (particle.GetType()==I3Particle::MuPlus);

    const bool isTau =
    (particle.GetType()==I3Particle::TauMinus) ||
    (particle.GetType()==I3Particle::TauPlus);

#ifdef I3PARTICLE_SUPPORTS_PDG_ENCODINGS
    if ((!isHadron) && (!isElectron) && (!isMuon) && (!isTau))
    {
        // if we don't know it but it has a pdg code,
        // it is probably a hadron..
        isHadron = true;
    }
#endif

    const double E = particle.GetEnergy()/I3Units::GeV;
    const double logE = std::max(0., std::log(E)); // protect against extremely low energies

    if (isElectron || isHadron) {
        const double nph=5.21*(0.924*I3Units::g/I3Units::cm3)/density;
        
        // Get cascade extension and fluctuation parameters
        I3SimConstants::ShowerParameters shower_params(particle.GetType(), E, density);
        double f = 1.;
        if (shower_params.emScaleSigma != 0.) {
            do {
                f=shower_params.emScale +
                    shower_params.emScaleSigma*randomService_->Gaus(0.,1.);
            } while((f<0.) || (1.<f));
        }
        
        const double meanNumPhotons = f*meanPhotonsPerMeter*nph*E;
        
        uint64_t numPhotons;
        if (meanNumPhotons > 1e7)
        {
            log_debug("HUGE EVENT: (cascade) meanNumPhotons=%f, nph=%f, E=%f (approximating possion by gaussian)", meanNumPhotons, nph, E);
            double numPhotonsDouble=0;
            do {
                numPhotonsDouble = randomService_->Gaus(meanNumPhotons, std::sqrt(meanNumPhotons));
            } while (numPhotonsDouble<0.);
            
            if (numPhotonsDouble > static_cast<double>(std::numeric_limits<uint64_t>::max()))
                log_fatal("Too many photons for counter. internal limitation.");
            numPhotons = static_cast<uint64_t>(numPhotonsDouble);
        }
        else
        {
            numPhotons = static_cast<uint64_t>(randomService_->Poisson(meanNumPhotons));
        }
        
        uint64_t usePhotonsPerStep = static_cast<uint64_t>(photonsPerStep_);
        if (static_cast<double>(numPhotons) > useHighPhotonsPerStepStartingFromNumPhotons_)
            usePhotonsPerStep = static_cast<uint64_t>(highPhotonsPerStep_);
        
        const uint64_t numSteps = numPhotons/usePhotonsPerStep;
        const uint64_t numPhotonsInLastStep = numPhotons%usePhotonsPerStep;

        log_trace("Generating %" PRIu64 " steps for cascade", numSteps);
        
        if (particle.GetShape() == I3Particle::CascadeSegment) {
            if (!(particle.GetLength() > 0))
                log_fatal_stream("Found a cascade segment with length "
                    << particle.GetLength() << ". This should not be.");
            MuonStepData_t segmentStepGenInfo;
            segmentStepGenInfo.particle=particle;
            segmentStepGenInfo.particle=particle;
            segmentStepGenInfo.particleIdentifier=identifier;
            segmentStepGenInfo.photonsPerStep=usePhotonsPerStep;
            segmentStepGenInfo.numSteps=numSteps;
            segmentStepGenInfo.numPhotonsInLastStep=numPhotonsInLastStep;
            segmentStepGenInfo.stepIsCascadeLike=true;
            segmentStepGenInfo.length=particle.GetLength();

            log_trace("== enqueue cascade segment");
            stepGenerationQueue_.push_back(segmentStepGenInfo);
        } else {
            CascadeStepData_t cascadeStepGenInfo;
            cascadeStepGenInfo.particle=particle;
            cascadeStepGenInfo.particle=particle;
            cascadeStepGenInfo.particleIdentifier=identifier;
            cascadeStepGenInfo.photonsPerStep=usePhotonsPerStep;
            cascadeStepGenInfo.numSteps=numSteps;
            cascadeStepGenInfo.numPhotonsInLastStep=numPhotonsInLastStep;
            cascadeStepGenInfo.pa=shower_params.a;
            cascadeStepGenInfo.pb=(useCascadeExtension_ ? shower_params.b : 0.);
        
            log_trace("== enqueue cascade");
            stepGenerationQueue_.push_back(cascadeStepGenInfo);
        }
    } else if (isMuon || isTau) {
        // TODO: for now, treat muons and taus (with lengths after MMC) the same.
        // This is compatible to what hit-maker does, but is of course not the right thing
        // to do. (Hit-maker treats all "tracks" the same and uses muon-tables for them.)
        
        const double length = std::isnan(particle.GetLength())?(2000.*I3Units::m):(particle.GetLength());
        
        if (std::isnan(particle.GetLength()))
            log_warn("Muon without length found! Assigned a length of 2000m.");

        log_trace("Parameterizing muon (ID=(%" PRIu64 "/%i)) with E=%fTeV, length=%fm",
                  particle.GetMajorID(), particle.GetMinorID(),
                  E*I3Units::GeV/I3Units::TeV, length/I3Units::m);
        
        // calculation the way it's done by PPC (updated June 2018)
        const double extr=1+std::max(0.0, 0.1880+0.0206*logE);
        const double muonFraction = 1./extr;
        
        const double meanNumPhotonsTotal = meanPhotonsPerMeter*(length/I3Units::m)*extr;
        
        log_trace("meanNumPhotonsTotal=%f", meanNumPhotonsTotal);
        
        const double meanNumPhotonsFromMuon = meanNumPhotonsTotal*muonFraction;
        
        uint64_t numPhotonsFromMuon;
        if (meanNumPhotonsFromMuon > 1e7)
        {
            log_debug("HUGE EVENT: (muon [   muon-like steps]) meanNumPhotons=%f, E=%f, len=%fm (approximating possion by gaussian)", meanNumPhotonsFromMuon, E, length/I3Units::m);
            double numPhotonsDouble=0;
            do {
                numPhotonsDouble = randomService_->Gaus(meanNumPhotonsFromMuon, std::sqrt(meanNumPhotonsFromMuon));
            } while (numPhotonsDouble<0.);
            
            if (numPhotonsDouble > static_cast<double>(std::numeric_limits<uint64_t>::max()))
                log_fatal("Too many photons for counter. internal limitation.");
            numPhotonsFromMuon = static_cast<uint64_t>(numPhotonsDouble);
        }
        else
        {
            numPhotonsFromMuon = static_cast<uint64_t>(randomService_->Poisson(meanNumPhotonsFromMuon));
        }
        log_trace("Generating %" PRIu64 " muon-steps for muon (mean=%f)", numPhotonsFromMuon, meanNumPhotonsFromMuon);

        const double meanNumPhotonsFromCascades = meanNumPhotonsTotal*(1.-muonFraction);

        log_trace("Generating a mean of %f cascade-steps for muon", meanNumPhotonsFromCascades);


        uint64_t numPhotonsFromCascades;
        if (meanNumPhotonsFromCascades > 1e7)
        {
            log_debug("HUGE EVENT: (muon [cascade-like steps]) meanNumPhotons=%f, E=%f, len=%fm (approximating possion by gaussian)", meanNumPhotonsFromCascades, E, length/I3Units::m);
            double numPhotonsDouble=0;
            do {
                numPhotonsDouble = randomService_->Gaus(meanNumPhotonsFromCascades, std::sqrt(meanNumPhotonsFromCascades));
            } while (numPhotonsDouble<0.);
            
            if (numPhotonsDouble > static_cast<double>(std::numeric_limits<uint64_t>::max()))
                log_fatal("Too many photons for counter. internal limitation.");
            numPhotonsFromCascades = static_cast<uint64_t>(numPhotonsDouble);
        }
        else
        {
            numPhotonsFromCascades = static_cast<uint64_t>(randomService_->Poisson(meanNumPhotonsFromCascades));
        }
        log_trace("Generating %" PRIu64 " cascade-steps for muon (mean=%f)", numPhotonsFromCascades, meanNumPhotonsFromCascades);

         
        // steps from muon
        
        uint64_t usePhotonsPerStep = static_cast<uint64_t>(photonsPerStep_);
        if (static_cast<double>(numPhotonsFromMuon) > useHighPhotonsPerStepStartingFromNumPhotons_)
            usePhotonsPerStep = static_cast<uint64_t>(highPhotonsPerStep_);
        
        const uint64_t numStepsFromMuon = numPhotonsFromMuon/usePhotonsPerStep;
        const uint64_t numPhotonsFromMuonInLastStep = numPhotonsFromMuon%usePhotonsPerStep;
        
        log_trace("Generating %" PRIu64 " steps for muon", numStepsFromMuon);

        MuonStepData_t muonStepGenInfo;
        muonStepGenInfo.particle=particle;
        muonStepGenInfo.particleIdentifier=identifier;
        muonStepGenInfo.photonsPerStep=usePhotonsPerStep;
        muonStepGenInfo.numSteps=numStepsFromMuon;
        muonStepGenInfo.numPhotonsInLastStep=numPhotonsFromMuonInLastStep;
        muonStepGenInfo.stepIsCascadeLike=false;
        muonStepGenInfo.length=length;
        log_trace("== enqueue muon (muon-like)");
        stepGenerationQueue_.push_back(muonStepGenInfo);
        
        log_trace("Generate %lu steps for E=%fGeV, l=%fm. (muon[muon])", static_cast<unsigned long>((numStepsFromMuon+((numPhotonsFromMuonInLastStep>0)?1:0))), E, length/I3Units::m);
        
        
        // steps from cascade
        
        usePhotonsPerStep = static_cast<uint64_t>(photonsPerStep_);
        if (static_cast<double>(numPhotonsFromCascades) > useHighPhotonsPerStepStartingFromNumPhotons_)
            usePhotonsPerStep = static_cast<uint64_t>(highPhotonsPerStep_);

        
        const uint64_t numStepsFromCascades = numPhotonsFromCascades/usePhotonsPerStep;
        const uint64_t numPhotonsFromCascadesInLastStep = numStepsFromCascades%usePhotonsPerStep;

        //MuonStepData_t muonStepGenInfo;
        muonStepGenInfo.particle=particle;
        muonStepGenInfo.particleIdentifier=identifier;
        muonStepGenInfo.photonsPerStep=usePhotonsPerStep;
        muonStepGenInfo.numSteps=numStepsFromCascades;
        muonStepGenInfo.numPhotonsInLastStep=numPhotonsFromCascadesInLastStep;
        muonStepGenInfo.stepIsCascadeLike=true;
        muonStepGenInfo.length=length;
        log_trace("== enqueue muon (cascade-like)");
        stepGenerationQueue_.push_back(muonStepGenInfo);
        
        log_trace("Generate %u steps for E=%fGeV, l=%fm. (muon[cascade])", static_cast<unsigned int>((numStepsFromCascades+((numPhotonsFromCascadesInLastStep>0)?1:0))), E, length/I3Units::m);
        
    } else {
#ifdef I3PARTICLE_SUPPORTS_PDG_ENCODINGS
        log_fatal("I3CLSimLightSourceToStepConverterPPC cannot handle a %s. (pdg %u)",
                  particle.GetTypeString().c_str(), particle.GetPdgEncoding());
#else
        log_fatal("I3CLSimLightSourceToStepConverterPPC cannot handle a %s.",
                  particle.GetTypeString().c_str());
#endif
    }
}

void I3CLSimLightSourceToStepConverterPPC::EnqueueBarrier()
{
    if (!initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterPPC is not initialized!");

    if (barrier_is_enqueued_)
        throw I3CLSimLightSourceToStepConverter_exception("A barrier is already enqueued!");

    // actually enqueue the barrier
    log_trace("== enqueue barrier");
    stepGenerationQueue_.push_back(BarrierData_t());
    barrier_is_enqueued_=true;
}

bool I3CLSimLightSourceToStepConverterPPC::BarrierActive() const
{
    if (!initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterPPC is not initialized!");

    return barrier_is_enqueued_;
}

bool I3CLSimLightSourceToStepConverterPPC::MoreStepsAvailable() const
{
    if (!initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterPPC is not initialized!");

    if (stepGenerationQueue_.size() > 0) return true;
    return false;
}

I3CLSimLightSourceToStepConverterPPC::MakeSteps_visitor::MakeSteps_visitor
(uint64_t &rngState, uint32_t rngA,
 uint64_t maxNumStepsPerStepSeries,
 GenerateStepPreCalculator &preCalc)
:rngState_(rngState), rngA_(rngA),
maxNumStepsPerStepSeries_(maxNumStepsPerStepSeries),
preCalc_(preCalc)
{;}

void I3CLSimLightSourceToStepConverterPPC::MakeSteps_visitor::FillStep
(I3CLSimLightSourceToStepConverterPPC::CascadeStepData_t &data,
 I3CLSimStep &newStep,
 uint64_t photonsPerStep,
 double particleDir_x, double particleDir_y, double particleDir_z) const
{
    const double longitudinalPos = data.pb*I3CLSimLightSourceToStepConverterUtils::gammaDistributedNumber(data.pa, rngState_, rngA_)*I3Units::m;
    GenerateStep(newStep,
                 data.particle,
                 particleDir_x, particleDir_y, particleDir_z,
                 data.particleIdentifier,
                 photonsPerStep,
                 longitudinalPos,
                 preCalc_);
}

void I3CLSimLightSourceToStepConverterPPC::MakeSteps_visitor::FillStep
(I3CLSimLightSourceToStepConverterPPC::MuonStepData_t &data,
 I3CLSimStep &newStep,
 uint64_t photonsPerStep,
 double particleDir_x, double particleDir_y, double particleDir_z) const
{
    if (data.stepIsCascadeLike) {
        const double longitudinalPos = mwcRngRandomNumber_co(rngState_, rngA_)*data.length;
        GenerateStep(newStep,
                     data.particle,
                     particleDir_x, particleDir_y, particleDir_z,
                     data.particleIdentifier,
                     photonsPerStep,
                     longitudinalPos,
                     preCalc_);
    } else {
        GenerateStepForMuon(newStep,
                            data.particle,
                            particleDir_x, particleDir_y, particleDir_z,
                            data.particleIdentifier,
                            photonsPerStep,
                            data.length);
    }
}

template <typename T>
std::pair<I3CLSimStepSeriesConstPtr, bool>
I3CLSimLightSourceToStepConverterPPC::MakeSteps_visitor::operator()
(T &data) const
{
    I3CLSimStepSeriesPtr currentStepSeries(new I3CLSimStepSeries());
    
    uint64_t useNumSteps = data.numSteps;
    if (useNumSteps > maxNumStepsPerStepSeries_) useNumSteps=maxNumStepsPerStepSeries_;
    
    const double particleDir_x = data.particle.GetDir().GetX();
    const double particleDir_y = data.particle.GetDir().GetY();
    const double particleDir_z = data.particle.GetDir().GetZ();
    
    // make useNumSteps steps
    for (uint64_t i=0; i<useNumSteps; ++i)
    {
        currentStepSeries->push_back(I3CLSimStep());
        I3CLSimStep &newStep = currentStepSeries->back();
        FillStep(data, newStep, data.photonsPerStep, particleDir_x, particleDir_y, particleDir_z);
    }
    
    // reduce the number of steps that have still to be processed
    data.numSteps -= useNumSteps;
    
    if ((data.numSteps==0) && (useNumSteps<maxNumStepsPerStepSeries_))
    {
        // make the last step (with a possible different number of photons than all the others)
        
        if (data.numPhotonsInLastStep > 0)
        {
            currentStepSeries->push_back(I3CLSimStep());
            I3CLSimStep &newStep = currentStepSeries->back();
            FillStep(data, newStep, data.numPhotonsInLastStep, particleDir_x, particleDir_y, particleDir_z);
        }
        
        // we are finished with this entry, it can be removed (signal this using the return value's .second entry)
        return std::make_pair(currentStepSeries, true);
    }
    else
    {
        // we are not finished with this entry, there are more steps to come..
        return std::make_pair(currentStepSeries, false);
    }
}

// specialization for BarrierData_t
template <>
std::pair<I3CLSimStepSeriesConstPtr, bool>
I3CLSimLightSourceToStepConverterPPC::MakeSteps_visitor::operator()
(I3CLSimLightSourceToStepConverterPPC::BarrierData_t &data) const
{
    log_trace("== end barrier in visitor::operator()");
    return make_pair(I3CLSimStepSeriesConstPtr(), true); // true=="entry can be removed from the queue"
}


I3CLSimStepSeriesConstPtr I3CLSimLightSourceToStepConverterPPC::MakeSteps(bool &barrierWasReset)
{
    barrierWasReset=false;
    
    if (stepGenerationQueue_.empty()) return I3CLSimStepSeriesConstPtr(); // queue is empty
    
    StepData_t &currentElement = stepGenerationQueue_.front();
    
    //  Let the visitor convert it into steps (the step pointer will be NULL if it is a barrier)
    std::pair<I3CLSimStepSeriesConstPtr, bool> retval =
    boost::apply_visitor(MakeSteps_visitor(rngState_, rngA_, maxBunchSize_, *preCalc_), currentElement);
    
    I3CLSimStepSeriesConstPtr &steps = retval.first;
    const bool entryCanBeRemoved = retval.second;
    
    if (entryCanBeRemoved) stepGenerationQueue_.pop_front();
    
    // steps==NULL means a barrier was reset. Return an empty list of 
    if (!steps) {
        barrierWasReset=true;
        return I3CLSimStepSeriesConstPtr(new I3CLSimStepSeries());
    } else {
        return steps;
    }
}
    
    

I3CLSimStepSeriesConstPtr I3CLSimLightSourceToStepConverterPPC::GetConversionResultWithBarrierInfo(bool &barrierWasReset, double timeout)
{
    if (!initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterPPC is not initialized!");
    
    barrierWasReset=false;
    
    if (stepGenerationQueue_.empty())
    {
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterPPC: no particle is enqueued!");
        return I3CLSimStepSeriesConstPtr();
    }
    
    I3CLSimStepSeriesConstPtr returnSteps = MakeSteps(barrierWasReset);
    if (!returnSteps) log_fatal("logic error. returnSteps==NULL");

    if (barrierWasReset) {
        if (!barrier_is_enqueued_)
            log_fatal("logic error: barrier encountered, but enqueued flag is false.");
        
        barrier_is_enqueued_=false;
    }
    
    return returnSteps;
}




/////// HELPERS

I3CLSimLightSourceToStepConverterPPC::GenerateStepPreCalculator::GenerateStepPreCalculator(I3RandomServicePtr randomService,
                                                     double angularDist_a,
                                                     double angularDist_b,
                                                     std::size_t numberOfValues)
:
one_over_angularDist_a_(1./angularDist_a),
angularDist_b_(angularDist_b),
angularDist_I_(1.-std::exp(-angularDist_b*std::pow(2., angularDist_a)) ),
numberOfValues_(numberOfValues),
index_(numberOfValues),
queueFromFeederThreads_(10)  // size 10 for 5 threads
{
    const unsigned int numFeederThreads = 4;
    const uint32_t rngAs[8] = { // numbers taken from Numerical Recipies
        3874257210,
        2936881968,
        2811536238,
        2654432763,
        4294957665,
        4294963023,
        4162943475,
        3947008974,
    };
    
    // start threads
    for (unsigned int i=0;i<numFeederThreads;++i)
    {
        const uint64_t rngState = mwcRngInitState(randomService, rngAs[i]);
        boost::shared_ptr<boost::thread> newThread(new boost::thread(boost::bind(&I3CLSimLightSourceToStepConverterPPC::GenerateStepPreCalculator::FeederThread, this, i, rngState, rngAs[i])));
        feederThreads_.push_back(newThread);
    }

}

I3CLSimLightSourceToStepConverterPPC::GenerateStepPreCalculator::~GenerateStepPreCalculator()
{
    // terminate threads
    
    for (std::size_t i=0;i<feederThreads_.size();++i)
    {
        if (!feederThreads_[i]) continue;
        if (!feederThreads_[i]->joinable()) continue;

        log_debug("Stopping the worker thread #%zu", i);
        feederThreads_[i]->interrupt();
        
        // wait for the thread to stop (not indefinitely, just to be sure)
        const bool did_join = feederThreads_[i]->timed_join(boost::posix_time::seconds(10));
        
        if (did_join) {
            log_debug("Worker thread #%zu stopped.", i);
        } else {
            log_warn("Worker thread #%zu did not stop. leaking memory.", i);
        }
    }

    feederThreads_.clear();
}

void I3CLSimLightSourceToStepConverterPPC::GenerateStepPreCalculator::FeederThread(unsigned int threadId,
                                                                                uint64_t initialRngState,
                                                                                uint32_t rngA)
{
    // set up storage
    uint64_t rngState = initialRngState;
    
    for (;;)
    {
        // make a bunch of steps
        boost::shared_ptr<queueVector_t> outputVector(new queueVector_t());
        outputVector->reserve(numberOfValues_);
        
        // calculate values
        for (std::size_t i=0;i<numberOfValues_;++i)
        {
            const double cos_val=std::max(1.-std::pow(-std::log(1.-mwcRngRandomNumber_co(rngState, rngA)*angularDist_I_)/angularDist_b_, one_over_angularDist_a_), -1.);
            const double sin_val=std::sqrt(1.-cos_val*cos_val);
            const double random_value = mwcRngRandomNumber_co(rngState, rngA);
            
            outputVector->push_back(std::make_pair(std::make_pair(sin_val, cos_val), random_value));
        }

        try 
        {
            // this blocks in case the queue is full
            queueFromFeederThreads_.Put(outputVector);
            log_debug("thread %u just refilled the queue", threadId);
        }
        catch(boost::thread_interrupted &i)
        {
            break;
        }
    }
}

void I3CLSimLightSourceToStepConverterPPC::GenerateStepPreCalculator::RegenerateValues()
{
    currentVector_ = queueFromFeederThreads_.Get();
    log_trace("queueSize=%zu", queueFromFeederThreads_.size());
    index_=0;
}




void I3CLSimLightSourceToStepConverterPPC::GenerateStep(I3CLSimStep &newStep,
                                                        const I3Particle &p,
                                                        double particleDir_x, double particleDir_y, double particleDir_z,
                                                        uint32_t identifier,
                                                        uint32_t photonsPerStep,
                                                        const double &longitudinalPos,
                                                        GenerateStepPreCalculator &preCalc)
{
    double angular_cos, angular_sin, random_value;
    preCalc.GetAngularCosSinValue(angular_cos, angular_sin, random_value);
    
    double step_dx = particleDir_x;
    double step_dy = particleDir_y;
    double step_dz = particleDir_z;
    
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
    newStep.SetSourceType(0); // cherenkov emission
    
    // rotate in-place
    scatterDirectionByAngle(angular_cos, angular_sin,
                            step_dx, step_dy, step_dz,
                            random_value);
    
    newStep.SetDir(step_dx, step_dy, step_dz);
}

void I3CLSimLightSourceToStepConverterPPC::GenerateStepForMuon(I3CLSimStep &newStep,
                                                               const I3Particle &p,
                                                               double particleDir_x, double particleDir_y, double particleDir_z,
                                                               uint32_t identifier,
                                                               uint32_t photonsPerStep,
                                                               double length)
{
    // set all values
    newStep.SetPosX(p.GetX());
    newStep.SetPosY(p.GetY());
    newStep.SetPosZ(p.GetZ());
    newStep.SetDir(particleDir_x, particleDir_y, particleDir_z);
    newStep.SetTime(p.GetTime());
    
    newStep.SetLength(length);
    newStep.SetNumPhotons(photonsPerStep);
    newStep.SetWeight(1.);
    newStep.SetBeta(1.);
    newStep.SetID(identifier);
    newStep.SetSourceType(0); // cherenkov emission
}


