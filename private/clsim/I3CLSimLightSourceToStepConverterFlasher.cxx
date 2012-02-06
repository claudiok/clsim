#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>

#include <cmath>
#include <boost/math/special_functions/fpclassify.hpp>
#include <limits>

#include "clsim/I3CLSimLightSourceToStepConverterFlasher.h"

#include "clsim/I3CLSimWlenDependentValue.h"

#include "clsim/I3CLSimLightSourceToStepConverterUtils.h"
using namespace I3CLSimLightSourceToStepConverterUtils;


const uint32_t I3CLSimLightSourceToStepConverterFlasher::default_photonsPerStep=2000;


I3CLSimLightSourceToStepConverterFlasher::I3CLSimLightSourceToStepConverterFlasher
(I3CLSimWlenDependentValueConstPtr flasherSpectrumNoBias,
 I3CLSimSpectrumTablePtr spectrumTable,
 uint32_t photonsPerStep)
:
initialized_(false),
barrier_is_enqueued_(false),
bunchSizeGranularity_(1),
maxBunchSize_(512000),
photonsPerStep_(photonsPerStep),
flasherSpectrumNoBias_(flasherSpectrumNoBias)
{
    if (photonsPerStep_<=0)
        throw I3CLSimLightSourceToStepConverter_exception("photonsPerStep may not be <= 0!");
    
    if (!flasherSpectrumNoBias_)
        throw I3CLSimLightSourceToStepConverter_exception("You have to provide a non-NULL spectrum!");
    
    if (!spectrumTable)
        throw I3CLSimLightSourceToStepConverter_exception("You have to provide a non-NULL spectrum table!");
    
    // register the spectrum in the table and retain the index / "source type"
    std::size_t spectrumSourceTypeIndex = spectrumTable->append(flasherSpectrumNoBias_);
    
    if (spectrumSourceTypeIndex >= 256)
        throw I3CLSimLightSourceToStepConverter_exception("No more than 255 distinct spectra are allowed!");
    
    spectrumSourceTypeIndex_ = static_cast<uint8_t>(spectrumSourceTypeIndex);
    
    log_debug("flasher spectrum registered in table. got index %zu.", spectrumSourceTypeIndex);
}

I3CLSimLightSourceToStepConverterFlasher::~I3CLSimLightSourceToStepConverterFlasher()
{

}

void I3CLSimLightSourceToStepConverterFlasher::Initialize()
{
    if (initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterFlasher already initialized!");

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
    
    // make a copy of the medium properties
    {
        I3CLSimMediumPropertiesConstPtr copiedMediumProperties(new I3CLSimMediumProperties(*mediumProperties_));
        mediumProperties_ = copiedMediumProperties;
    }

    double minWlen=flasherSpectrumNoBias_->GetMinWlen();
    double maxWlen=flasherSpectrumNoBias_->GetMaxWlen();
    if (boost::math::isinf(minWlen)) minWlen=mediumProperties_->GetMinWavelength();
    if (boost::math::isinf(maxWlen)) maxWlen=mediumProperties_->GetMaxWavelength();
    
    // calculate the photon number correction factor for the selected spectrum and bias
    photonNumberCorrectionFactorForBias_ = 
    PhotonNumberCorrectionFactorAfterBias(*(flasherSpectrumNoBias_),
                                          *(wlenBias_),
                                          minWlen,
                                          maxWlen);
    log_debug("photon number correction factor for spectrum and bias is %f", photonNumberCorrectionFactorForBias_);

    initialized_=true;
}


bool I3CLSimLightSourceToStepConverterFlasher::IsInitialized() const
{
    return initialized_;
}

void I3CLSimLightSourceToStepConverterFlasher::SetBunchSizeGranularity(uint64_t num)
{
    if (initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterFlasher already initialized!");
    
    if (num<=0)
        throw I3CLSimLightSourceToStepConverter_exception("BunchSizeGranularity of 0 is invalid!");

    if (num!=1)
        throw I3CLSimLightSourceToStepConverter_exception("A BunchSizeGranularity != 1 is currently not supported!");

    bunchSizeGranularity_=num;
}

void I3CLSimLightSourceToStepConverterFlasher::SetMaxBunchSize(uint64_t num)
{
    if (initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterFlasher already initialized!");

    if (num<=0)
        throw I3CLSimLightSourceToStepConverter_exception("MaxBunchSize of 0 is invalid!");

    maxBunchSize_=num;
}

void I3CLSimLightSourceToStepConverterFlasher::SetWlenBias(I3CLSimWlenDependentValueConstPtr wlenBias)
{
    if (initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterFlasher already initialized!");
    
    wlenBias_=wlenBias;
}

void I3CLSimLightSourceToStepConverterFlasher::SetRandomService(I3RandomServicePtr random)
{
    if (initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterFlasher already initialized!");
    
    randomService_=random;
}

void I3CLSimLightSourceToStepConverterFlasher::SetMediumProperties(I3CLSimMediumPropertiesConstPtr mediumProperties)
{
    if (initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterFlasher already initialized!");

    mediumProperties_=mediumProperties;
}

void I3CLSimLightSourceToStepConverterFlasher::EnqueueLightSource(const I3CLSimLightSource &lightSource, uint32_t identifier)
{
    if (!initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterFlasher is not initialized!");

    if (barrier_is_enqueued_)
        throw I3CLSimLightSourceToStepConverter_exception("A barrier is enqueued! You must receive all steps before enqueuing a new particle.");

    if (lightSource.GetType() != I3CLSimLightSource::Flasher)
        throw I3CLSimLightSourceToStepConverter_exception("The I3CLSimLightSourceToStepConverterFlasher parameterization only works on flashers.");
    
    const I3CLSimFlasherPulse &flasherPulse = lightSource.GetFlasherPulse();

    // just skip the entry if there are no photons to generate
    if (flasherPulse.GetNumberOfPhotonsNoBias() <= 0.) return;
    
    const double numPhotonsWithBias = flasherPulse.GetNumberOfPhotonsNoBias()*photonNumberCorrectionFactorForBias_;
    if (numPhotonsWithBias <= 0.) return;
    
    uint64_t numPhotons;
    
    if (numPhotonsWithBias > 1e6)
    {
        log_debug("HUGE EVENT: (e-m) numPhotonsWithBias=%f (approximating possion by gaussian)", numPhotonsWithBias);
        double numPhotonsDouble=0;
        do {
            numPhotonsDouble = randomService_->Gaus(numPhotonsWithBias, std::sqrt(numPhotonsWithBias));
        } while (numPhotonsDouble<0.);
        
        if (numPhotonsDouble > static_cast<double>(std::numeric_limits<uint64_t>::max()))
            log_fatal("Too many photons for counter. hit an internal limitation.");
        numPhotons = static_cast<uint64_t>(numPhotonsDouble);
    }
    else
    {
        numPhotons = static_cast<uint64_t>(randomService_->Poisson(numPhotonsWithBias));
    }
    if (numPhotons==0) return;
    

    log_warn("Generating %zu photons for flasher (after bias). Requested number was %f, with bias %f",
             static_cast<std::size_t>(numPhotons),
             flasherPulse.GetNumberOfPhotonsNoBias(),
             numPhotonsWithBias);
    
    LightSourceData_t newEntry;
    newEntry.isBarrier = false;
    newEntry.flasherPulse = flasherPulse;
    newEntry.identifier = identifier;
    newEntry.numPhotonsWithBias = numPhotons;
    inputQueue_.push_back(newEntry);
}

void I3CLSimLightSourceToStepConverterFlasher::EnqueueBarrier()
{
    if (!initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterFlasher is not initialized!");

    if (barrier_is_enqueued_)
        throw I3CLSimLightSourceToStepConverter_exception("A barrier is already enqueued!");

    // actually enqueue the barrier
    log_trace("== enqueue barrier");
    LightSourceData_t newEntry;
    newEntry.isBarrier = true;
    inputQueue_.push_back(newEntry);

    barrier_is_enqueued_=true;
}

bool I3CLSimLightSourceToStepConverterFlasher::BarrierActive() const
{
    if (!initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterFlasher is not initialized!");

    return barrier_is_enqueued_;
}

bool I3CLSimLightSourceToStepConverterFlasher::MoreStepsAvailable() const
{
    if (!initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterFlasher is not initialized!");

    if (inputQueue_.size() > 0) return true;
    return false;
}


I3CLSimStepSeriesConstPtr I3CLSimLightSourceToStepConverterFlasher::GetConversionResultWithBarrierInfo(bool &barrierWasReset, double timeout)
{
    if (!initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterFlasher is not initialized!");
    
    barrierWasReset=false;
    
    if (inputQueue_.empty())
    {
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterFlasher: no flasher pulse is in queue!");
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



// the actual work is done here:
I3CLSimStepSeriesConstPtr I3CLSimLightSourceToStepConverterFlasher::MakeSteps(bool &barrierWasReset)
{
    barrierWasReset=false;
    
    if (inputQueue_.empty()) return I3CLSimStepSeriesConstPtr(); // queue is empty
    
    LightSourceData_t &currentElement = inputQueue_.front();

    // barrier?
    if (currentElement.isBarrier) {
        barrierWasReset=true;
        return I3CLSimStepSeriesConstPtr(new I3CLSimStepSeries());
    }

    I3CLSimStepSeriesPtr outputSteps(new I3CLSimStepSeries());
    
    bool entryCanBeRemoved;
    
    uint64_t numSteps;                  // number of steps to generate
    uint64_t numAppendedDummySteps;     // number of dummy steps to append in order to achieve the correct steps vector size granularity
    uint32_t numPhotonsInLastStep;      // number of photons in the last step
    // (the number of photons in all other steps is photonsPerStep_)
    
    const uint64_t maxPhotonsPerResult = maxBunchSize_*static_cast<uint64_t>(photonsPerStep_);
    if (currentElement.numPhotonsWithBias >= maxPhotonsPerResult) {
        numSteps = maxBunchSize_;
        numPhotonsInLastStep = photonsPerStep_; // nothing special for the last step
        
        numAppendedDummySteps=0;
        
        currentElement.numPhotonsWithBias -= numSteps*static_cast<uint64_t>(photonsPerStep_);
        entryCanBeRemoved = (currentElement.numPhotonsWithBias == 0);
        
        outputSteps->reserve(maxBunchSize_);
    } else {
        if (currentElement.numPhotonsWithBias <= static_cast<uint64_t>(photonsPerStep_)) {
            // only a single step
            numSteps = 1;
            numPhotonsInLastStep=static_cast<uint32_t>(currentElement.numPhotonsWithBias);
        } else {
            numSteps = currentElement.numPhotonsWithBias / static_cast<uint64_t>(photonsPerStep_);
            numPhotonsInLastStep = static_cast<uint32_t>(currentElement.numPhotonsWithBias % static_cast<uint64_t>(photonsPerStep_));
            if (numPhotonsInLastStep>0) {
                // one step more with less photons
                ++numSteps;
            }
        }

        entryCanBeRemoved=true;
        currentElement.numPhotonsWithBias = 0;

        const uint64_t modulo = numSteps % bunchSizeGranularity_;
        numAppendedDummySteps=0;
        if (modulo>0) {
            // we need some dummy steps
            
            numAppendedDummySteps = bunchSizeGranularity_-modulo;
        }
        
        outputSteps->reserve(numSteps+numAppendedDummySteps);
    }
    
    // now make the steps
    for (uint64_t i=0;i<numSteps;++i)
    {
        const uint32_t numberOfPhotonsForThisStep = (i==numSteps-1)?numPhotonsInLastStep:photonsPerStep_;
        if (numberOfPhotonsForThisStep==0) {++numAppendedDummySteps; continue;}
        
        outputSteps->push_back(I3CLSimStep());
        I3CLSimStep &newStep = outputSteps->back();

        FillStep(newStep,
                 numberOfPhotonsForThisStep,
                 currentElement.flasherPulse,
                 currentElement.identifier
                );
    }

    // and the dummy steps
    for (uint64_t i=0;i<numAppendedDummySteps;++i)
    {
        outputSteps->push_back(I3CLSimStep());
        I3CLSimStep &newStep = outputSteps->back();

        newStep.SetPosX(0.); newStep.SetPosY(0.); newStep.SetPosZ(0.);
        newStep.SetTime(0.);
        newStep.SetDirTheta(0.); newStep.SetDirPhi(0.);
        newStep.SetLength(0.);
        newStep.SetBeta(1.);
        newStep.SetNumPhotons(0);
        newStep.SetWeight(0);
        newStep.SetSourceType(0);
        newStep.SetID(currentElement.identifier);
    }
    
    if (entryCanBeRemoved) inputQueue_.pop_front();
    
    return outputSteps;
}


void I3CLSimLightSourceToStepConverterFlasher::FillStep(I3CLSimStep &step,
                                                        uint32_t numberOfPhotons,
                                                        const I3CLSimFlasherPulse &flasherPulse,
                                                        uint32_t identifier)
{
    step.SetPos(flasherPulse.GetPos());
    step.SetTime(flasherPulse.GetTime());
    step.SetDir(flasherPulse.GetDir());
    step.SetLength(0.);
    step.SetBeta(1.);
    step.SetNumPhotons(numberOfPhotons);
    step.SetWeight(1.);
    step.SetID(identifier);
    step.SetSourceType(spectrumSourceTypeIndex_);
}


