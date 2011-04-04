#include "TrkUserEventInformation.hh"
#include "TrkDetectorConstruction.hh"

TrkUserEventInformation::TrkUserEventInformation(uint64_t maxBunchSize_,
                                                 I3CLSimStepStorePtr stepStore_,
                                                 const I3CLSimParticleParameterizationSeries &parameterizationAvailable_,
                                                 boost::shared_ptr<I3CLSimQueue<I3CLSimParticleToStepConverterGeant4::FromGeant4Pair_t> > queueFromGeant4_,
                                                 boost::this_thread::disable_interruption &threadDisabledInterruptionState_,
                                                 uint32_t currentExternalParticleID_,
                                                 double maxRefractiveIndex_)
:
abortRequested(false),
maxBunchSize(maxBunchSize_),
stepStore(stepStore_),
parameterizationAvailable(parameterizationAvailable_),
queueFromGeant4(queueFromGeant4_),
threadDisabledInterruptionState(threadDisabledInterruptionState_),
currentExternalParticleID(currentExternalParticleID_),
maxRefractiveIndex(maxRefractiveIndex_)
{
}

TrkUserEventInformation::~TrkUserEventInformation()
{
}

