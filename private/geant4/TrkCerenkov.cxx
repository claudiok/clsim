//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
////////////////////////////////////////////////////////////////////////
// Modified Cerenkov Radiation Class Definition 
////////////////////////////////////////////////////////////////////////
//
// File:        TrkCerenkov.hh  
// Description:	Discrete Process - Generation of Cerenkov Photons
// Based on:    G4Cerenkov.hh from the original Geant4 distribution
//              (v4.9.4)
//
////////////////////////////////////////////////////////////////////////

#include "G4ios.hh"
#include "G4Poisson.hh"
#include "G4EmProcessSubType.hh"

#include "G4LossTableManager.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ParticleDefinition.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4RunManager.hh"
#include "TrkUserEventInformation.hh"

#include "TrkCerenkov.hh"

#include "icetray/I3Units.h"

#include <boost/thread.hpp>


/////////////////////////
// Class Implementation  
/////////////////////////


/////////////////
// Constructors
/////////////////

TrkCerenkov::TrkCerenkov(const G4String& processName, G4ProcessType type)
: G4VProcess(processName, type)
{
    SetProcessSubType(fCerenkov);
    
    //fTrackSecondariesFirst = false;
    fMaxBetaChange = 0.;
    fMaxPhotons = 0;
    
    thePhysicsTable1 = NULL;
    thePhysicsTable2 = NULL;
    
    if (verboseLevel>0) {
        G4cout << GetProcessName() << " is created " << G4endl;
    }
    
    BuildThePhysicsTable();
}

// TrkCerenkov::TrkCerenkov(const TrkCerenkov &right)
// {
// }

////////////////
// Destructors
////////////////

TrkCerenkov::~TrkCerenkov() 
{
    if (thePhysicsTable1 != NULL) {
        thePhysicsTable1->clearAndDestroy();
        delete thePhysicsTable1;
    }

    if (thePhysicsTable2 != NULL) {
        thePhysicsTable2->clearAndDestroy();
        delete thePhysicsTable2;
    }
}

////////////
// Methods
////////////

// PostStepDoIt
// -------------
//
G4VParticleChange*
TrkCerenkov::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)

// This routine is called for each tracking Step of a charged particle
// in a radiator. A Poisson-distributed number of photons is generated
// according to the Cerenkov formula, distributed evenly along the track
// segment and uniformly azimuth w.r.t. the particle direction. The 
// parameters are then transformed into the Master Reference System, and 
// they are added to the particle change. 

{
    //////////////////////////////////////////////////////
    // Should we ensure that the material is dispersive?
    //////////////////////////////////////////////////////
    
    aParticleChange.Initialize(aTrack);
    
    const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
    const G4Material* aMaterial = aTrack.GetMaterial();
    
    G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
    G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();
    
    G4ThreeVector x0 = pPreStepPoint->GetPosition();
    G4ThreeVector p0 = aStep.GetDeltaPosition().unit();
    G4double t0 = pPreStepPoint->GetGlobalTime();
    
    G4MaterialPropertiesTable* aMaterialPropertiesTable =
    aMaterial->GetMaterialPropertiesTable();
    if (!aMaterialPropertiesTable) return pParticleChange;
    
    const G4MaterialPropertyVector* Rindex = 
    aMaterialPropertiesTable->GetProperty("RINDEX"); 
    if (!Rindex) return pParticleChange;
    
    // particle charge
    const G4double charge = aParticle->GetDefinition()->GetPDGCharge();
    
    // particle beta
    const G4double beta2 = pPostStepPoint->GetBeta();
    const G4double beta = (pPreStepPoint ->GetBeta() +
                           beta2)/2.;
    
    G4double MeanNumberOfPhotons = 
    GetAverageNumberOfPhotons(charge,beta,aMaterial,Rindex);
    
    if (MeanNumberOfPhotons <= 0.0) {
        
        // return unchanged particle and no secondaries
        
        aParticleChange.SetNumberOfSecondaries(0);
        
        return pParticleChange;
        
    }
    
    const G4double step_length = aStep.GetStepLength();
    
    MeanNumberOfPhotons = MeanNumberOfPhotons * step_length;
    
    aParticleChange.SetNumberOfSecondaries(0);

    G4int NumPhotons = (G4int) G4Poisson(MeanNumberOfPhotons);

    
    
    const double canHeight = 1000.*m;
    const double canRadius = 750.*m;
    if (fabs(x0.z()) > canHeight/2.) {
        // we are above or below the can
        NumPhotons=0;
    } else {
        G4double posRadius = std::sqrt(x0.x()*x0.x() + x0.y()*x0.y());
        if (posRadius > canRadius) {
            NumPhotons=0;
        }
    }
    
    
    if (NumPhotons <= 0) 
    {
        // return unchanged particle and no secondaries  
        aParticleChange.SetNumberOfSecondaries(0);
        return pParticleChange;
    }

	TrkUserEventInformation* eventInformation
	=(TrkUserEventInformation*)G4EventManager::GetEventManager()
	->GetConstCurrentEvent()->GetUserInformation();
    
    
    // add this step to the step store
    I3CLSimStepStorePtr stepStore = eventInformation->stepStore;

    {
        // insert a new step
        I3CLSimStep &newStep = stepStore->insert_new(NumPhotons); // insert @ NumPhotons
    
        // set all values
        newStep.SetPosX(x0.x()*I3Units::m/m);
        newStep.SetPosY(x0.y()*I3Units::m/m);
        newStep.SetPosZ(x0.z()*I3Units::m/m);
        newStep.SetTime(t0*I3Units::ns/ns);
        newStep.SetDir(p0.x(), p0.y(), p0.z());

        newStep.SetLength(step_length*I3Units::m/m);
        newStep.SetNumPhotons(NumPhotons);
        newStep.SetWeight(1.);
        newStep.SetBeta(beta);
        newStep.SetID(eventInformation->currentExternalParticleID);
    }
    
    // if the store size is large enough, flush some events to the external queue
    if (stepStore->size() >= eventInformation->maxBunchSize*2)
    {
        I3CLSimStepSeriesPtr steps(new I3CLSimStepSeries());

        stepStore->pop_bunch_to_vector(eventInformation->maxBunchSize, *steps);
        
        eventInformation->StopClock();
        
        {
            boost::this_thread::restore_interruption ri(eventInformation->threadDisabledInterruptionState);
            try {
                eventInformation->queueFromGeant4->Put(std::make_pair(steps, false));
            } catch(boost::thread_interrupted &i) {
                G4cout << "G4 thread was interrupted. shutting down Geant4!" << G4endl;
                
                eventInformation->abortRequested = true;
                G4RunManager::GetRunManager()->AbortRun();
            }
        }
        
        //const double stepsPerTime = static_cast<double>(eventInformation->maxBunchSize)/eventInformation->GetElapsedWallTime();
        //G4cout << "Geant4 just sent " << eventInformation->maxBunchSize<< " steps. => " << stepsPerTime/(1./I3Units::s) << " steps/second" << G4endl;
        
        eventInformation->StartClock();
    }
    else 
    {
        // check if the thread was interrupted
        {
            boost::this_thread::restore_interruption ri(eventInformation->threadDisabledInterruptionState);
            try {
                boost::this_thread::interruption_point();
            } catch(boost::thread_interrupted &i) {
                G4cout << "G4 thread was interrupted. shutting down Geant4!" << G4endl;

                eventInformation->abortRequested = true;
                G4RunManager::GetRunManager()->AbortRun();
            }
        }
    }
    stepStore.reset();
    
    
    // kill the particle if it has fallen below the threshold
    if (1./beta2 > eventInformation->maxRefractiveIndex)
    {
        aParticleChange.ProposeTrackStatus(fStopAndKill);
    }
    
    /*
    ////////////////////////////////////////////////////////////////
    
    aParticleChange.SetNumberOfSecondaries(NumPhotons);
    
    //if (fTrackSecondariesFirst) {
    //    if (aTrack.GetTrackStatus() == fAlive )
    //        aParticleChange.ProposeTrackStatus(fSuspend);
    //}
    
    ////////////////////////////////////////////////////////////////
    
    G4double Pmin = Rindex->GetMinPhotonEnergy();
    G4double Pmax = Rindex->GetMaxPhotonEnergy();
    G4double dp = Pmax - Pmin;
    
    G4double nMax = Rindex->GetMaxProperty();
    
    G4double BetaInverse = 1./beta;
    
    G4double maxCos = BetaInverse / nMax; 
    G4double maxSin2 = (1.0 - maxCos) * (1.0 + maxCos);
    
    const G4double beta1 = pPreStepPoint ->GetBeta();
    const G4double beta2 = pPostStepPoint->GetBeta();
    
    G4double MeanNumberOfPhotons1 =
    GetAverageNumberOfPhotons(charge,beta1,aMaterial,Rindex);
    G4double MeanNumberOfPhotons2 =
    GetAverageNumberOfPhotons(charge,beta2,aMaterial,Rindex);
    
    for (G4int i = 0; i < NumPhotons; i++) {
        
        // Determine photon energy
        
        G4double rand;
        G4double sampledEnergy, sampledRI; 
        G4double cosTheta, sin2Theta;
        
        // sample an energy
        
        do {
            rand = G4UniformRand();	
            sampledEnergy = Pmin + rand * dp; 
            sampledRI = Rindex->GetProperty(sampledEnergy);
            cosTheta = BetaInverse / sampledRI;  
            
            sin2Theta = (1.0 - cosTheta)*(1.0 + cosTheta);
            rand = G4UniformRand();	
            
        } while (rand*maxSin2 > sin2Theta);
        
        // Generate random position of photon on cone surface 
        // defined by Theta 
        
        rand = G4UniformRand();
        
        G4double phi = twopi*rand;
        G4double sinPhi = std::sin(phi);
        G4double cosPhi = std::cos(phi);
        
        // calculate x,y, and z components of photon energy
        // (in coord system with primary particle direction 
        //  aligned with the z axis)
        
        G4double sinTheta = std::sqrt(sin2Theta); 
        G4double px = sinTheta*cosPhi;
        G4double py = sinTheta*sinPhi;
        G4double pz = cosTheta;
        
        // Create photon momentum direction vector 
        // The momentum direction is still with respect
        // to the coordinate system where the primary
        // particle direction is aligned with the z axis  
        
        G4ParticleMomentum photonMomentum(px, py, pz);
        
        // Rotate momentum direction back to global reference
        // system 
        
        photonMomentum.rotateUz(p0);
        
        // Determine polarization of new photon 
        
        G4double sx = cosTheta*cosPhi;
        G4double sy = cosTheta*sinPhi; 
        G4double sz = -sinTheta;
        
        G4ThreeVector photonPolarization(sx, sy, sz);
        
        // Rotate back to original coord system 
        
        photonPolarization.rotateUz(p0);
        
        // Generate a new photon:
        
        G4DynamicParticle* aCerenkovPhoton =
        new G4DynamicParticle(G4OpticalPhoton::OpticalPhoton(), 
                              photonMomentum);
        aCerenkovPhoton->SetPolarization
        (photonPolarization.x(),
         photonPolarization.y(),
         photonPolarization.z());
        
        aCerenkovPhoton->SetKineticEnergy(sampledEnergy);
        
        // Generate new G4Track object:
        
        G4double delta, NumberOfPhotons, N;
        
        do {
            rand = G4UniformRand();
            delta = rand * aStep.GetStepLength();
            NumberOfPhotons = MeanNumberOfPhotons1 - delta *
            (MeanNumberOfPhotons1-MeanNumberOfPhotons2)/
            aStep.GetStepLength();
            N = G4UniformRand() *
            std::max(MeanNumberOfPhotons1,MeanNumberOfPhotons2);
        } while (N > NumberOfPhotons);
        
        G4double deltaTime = delta /
        ((pPreStepPoint->GetVelocity()+
          pPostStepPoint->GetVelocity())/2.);
        
        G4double aSecondaryTime = t0 + deltaTime;
        
        G4ThreeVector aSecondaryPosition =
        x0 + rand * aStep.GetDeltaPosition();
        
        G4Track* aSecondaryTrack = 
        new G4Track(aCerenkovPhoton,aSecondaryTime,aSecondaryPosition);
        
        aSecondaryTrack->SetTouchableHandle(
                                            aStep.GetPreStepPoint()->GetTouchableHandle());
        
        aSecondaryTrack->SetParentID(aTrack.GetTrackID());
        
        aParticleChange.AddSecondary(aSecondaryTrack);
    }
    
    if (verboseLevel>0) {
        G4cout << "\n Exiting from TrkCerenkov::DoIt -- NumberOfSecondaries = " 
        << aParticleChange.GetNumberOfSecondaries() << G4endl;
    }
    */
    
    return pParticleChange;
}

// BuildThePhysicsTable for the Cerenkov process
// ---------------------------------------------
//

void TrkCerenkov::BuildThePhysicsTable()
{
    if (thePhysicsTable1 != NULL) {
        thePhysicsTable1->clearAndDestroy();
        delete thePhysicsTable1;
        thePhysicsTable1=NULL;
    }
    
    if (thePhysicsTable2 != NULL) {
        thePhysicsTable2->clearAndDestroy();
        delete thePhysicsTable2;
        thePhysicsTable2=NULL;
    }

    
    const G4MaterialTable* theMaterialTable=
    G4Material::GetMaterialTable();
    G4int numOfMaterials = G4Material::GetNumberOfMaterials();
    
    // create new physics table
    
    thePhysicsTable1 = new G4PhysicsTable(numOfMaterials);
    thePhysicsTable2 = new G4PhysicsTable(numOfMaterials);
    
    // loop for materials
    
    for (G4int i=0 ; i < numOfMaterials; i++)
    {
        G4PhysicsOrderedFreeVector* aPhysicsOrderedFreeVector1 =
        new G4PhysicsOrderedFreeVector();

        G4PhysicsOrderedFreeVector* aPhysicsOrderedFreeVector2 =
        new G4PhysicsOrderedFreeVector();
        
        // Retrieve vector of refraction indices for the material
        // from the material's optical properties table 
        
        G4Material* aMaterial = (*theMaterialTable)[i];
        
        G4MaterialPropertiesTable* aMaterialPropertiesTable =
        aMaterial->GetMaterialPropertiesTable();
        
        if (aMaterialPropertiesTable) {
            
            G4MaterialPropertyVector* theRefractionIndexVector = 
            aMaterialPropertiesTable->GetProperty("RINDEX");
            
            if (theRefractionIndexVector) {
                
                // Retrieve the first refraction index in vector
                // of (photon energy, refraction index) pairs 
                
                theRefractionIndexVector->ResetIterator();
                ++(*theRefractionIndexVector);	// advance to 1st entry 
                
                G4double currentRI = theRefractionIndexVector->
                GetProperty();
                
                if (currentRI > 1.0) {
                    
                    // Create first (photon energy, Cerenkov Integral)
                    // pair  
                    
                    G4double currentPM = theRefractionIndexVector->
                    GetPhotonEnergy();
                    G4double currentCAI1 = 0.0;
                    G4double currentCAI2 = 0.0;
                    
                    aPhysicsOrderedFreeVector1->
                    InsertValues(currentPM , currentCAI1);
                    aPhysicsOrderedFreeVector2->
                    InsertValues(currentPM , currentCAI2);
                    
                    // Set previous values to current ones prior to loop
                    
                    G4double prevPM  = currentPM;
                    G4double prevCAI1 = currentCAI1;
                    G4double prevCAI2 = currentCAI2;
                    G4double prevRI  = currentRI;
                    
                    // loop over all (photon energy, refraction index)
                    // pairs stored for this material  
                    
                    while(++(*theRefractionIndexVector))
                    {
                        currentRI=theRefractionIndexVector->	
                        GetProperty();
                        
                        currentPM = theRefractionIndexVector->
                        GetPhotonEnergy();
                        
                        double currentBiasFactor=1.;
                        double prevBiasFactor=1.;
                        if (fWlenBias) {
                            currentBiasFactor = fWlenBias->GetValue(((h_Planck*c_light/currentPM)/nm)*I3Units::nanometer);
                            prevBiasFactor = fWlenBias->GetValue(((h_Planck*c_light/prevPM)/nm)*I3Units::nanometer);
                        }
                        
                        currentCAI1 = 0.5*(prevBiasFactor + currentBiasFactor);
                        currentCAI2 = 0.5*(prevBiasFactor/(prevRI*prevRI) +
                                           currentBiasFactor/(currentRI*currentRI));

                        currentCAI1 = prevCAI1 + 
                        (currentPM - prevPM) * currentCAI1;
                        currentCAI2 = prevCAI2 + 
                        (currentPM - prevPM) * currentCAI2;
                        
                        aPhysicsOrderedFreeVector1->
                        InsertValues(currentPM, currentCAI1);
                        aPhysicsOrderedFreeVector2->
                        InsertValues(currentPM, currentCAI2);
                        
                        prevPM  = currentPM;
                        prevCAI1 = currentCAI1;
                        prevCAI2 = currentCAI2;
                        prevRI  = currentRI;
                    }
                    
                }
            }
        }
        
        // The Cerenkov integral for a given material
        // will be inserted in thePhysicsTable2
        // according to the position of the material in
        // the material table. 
        
        thePhysicsTable1->insertAt(i,aPhysicsOrderedFreeVector1); 
        thePhysicsTable2->insertAt(i,aPhysicsOrderedFreeVector2); 
        
    }
}

// GetMeanFreePath
// ---------------
//

G4double TrkCerenkov::GetMeanFreePath(const G4Track&,
                                     G4double,
                                     G4ForceCondition*)
{
    return 1.;
}

G4double TrkCerenkov::PostStepGetPhysicalInteractionLength(
                                                          const G4Track& aTrack,
                                                          G4double,
                                                          G4ForceCondition* condition)
{
    *condition = NotForced;
    G4double StepLimit = DBL_MAX;
    
    const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
    const G4Material* aMaterial = aTrack.GetMaterial();
    const G4MaterialCutsCouple* couple = aTrack.GetMaterialCutsCouple();
    
    const G4double kineticEnergy = aParticle->GetKineticEnergy();
    const G4ParticleDefinition* particleType = aParticle->GetDefinition();
    const G4double mass = particleType->GetPDGMass();
    
    // particle beta
    const G4double beta = aParticle->GetTotalMomentum() /
    aParticle->GetTotalEnergy();
    // particle gamma
    const G4double gamma = 1./std::sqrt(1.-beta*beta);
    
    G4MaterialPropertiesTable* aMaterialPropertiesTable =
    aMaterial->GetMaterialPropertiesTable();
    
    const G4MaterialPropertyVector* Rindex = NULL;
    
    if (aMaterialPropertiesTable)
        Rindex = aMaterialPropertiesTable->GetProperty("RINDEX");
    
    G4double nMax;
    if (Rindex) {
        nMax = Rindex->GetMaxProperty();
    } else {
        return StepLimit;
    }
    
    G4double BetaMin = 1./nMax;
    if ( BetaMin >= 1. ) return StepLimit;
    
    G4double GammaMin = 1./std::sqrt(1.-BetaMin*BetaMin);
    
    if (gamma < GammaMin ) return StepLimit;
    
    G4double kinEmin = mass*(GammaMin-1.);
    
    G4double RangeMin = G4LossTableManager::Instance()->
    GetRange(particleType,
             kinEmin,
             couple);
    G4double Range    = G4LossTableManager::Instance()->
    GetRange(particleType,
             kineticEnergy,
             couple);
    
    G4double Step = Range - RangeMin;
    if (Step < 1.*um ) return StepLimit;
    
    if (Step > 0. && Step < StepLimit) StepLimit = Step; 
    
    // If user has defined an average maximum number of photons to
    // be generated in a Step, then calculate the Step length for
    // that number of photons. 
    
    if (fMaxPhotons > 0) {
        
        // particle charge
        const G4double charge = aParticle->
        GetDefinition()->GetPDGCharge();
        
        G4double MeanNumberOfPhotons = 
        GetAverageNumberOfPhotons(charge,beta,aMaterial,Rindex);
        
        G4double Step = 0.;
        if (MeanNumberOfPhotons > 0.0) Step = fMaxPhotons /
            MeanNumberOfPhotons;
        
        if (Step > 0. && Step < StepLimit) StepLimit = Step;
    }
    
    // If user has defined an maximum allowed change in beta per step
    if (fMaxBetaChange > 0.) {
        
        G4double dedx = G4LossTableManager::Instance()->
        GetDEDX(particleType,
                kineticEnergy,
                couple);
        
        G4double deltaGamma = gamma - 
        1./std::sqrt(1.-beta*beta*
                     (1.-fMaxBetaChange)*
                     (1.-fMaxBetaChange));
        
        G4double Step = mass * deltaGamma / dedx;
        
        if (Step > 0. && Step < StepLimit) StepLimit = Step;
        
    }
    
    *condition = StronglyForced;
    return StepLimit;
}

// GetAverageNumberOfPhotons
// -------------------------
// This routine computes the number of Cerenkov photons produced per
// GEANT-unit (millimeter) in the current medium. 
//             ^^^^^^^^^^

G4double 
TrkCerenkov::GetAverageNumberOfPhotons(const G4double charge,
                                      const G4double beta, 
                                      const G4Material* aMaterial,
                                      const G4MaterialPropertyVector* Rindex) const
{
    const G4double Rfact = 369.81/(eV * cm);
    
    if(beta <= 0.0)return 0.0;
    
    G4double BetaInverse = 1./beta;
    
    // Vectors used in computation of Cerenkov Angle Integral:
    // 	- Refraction Indices for the current material
    //	- new G4PhysicsOrderedFreeVector allocated to hold CAI's
    
    G4int materialIndex = aMaterial->GetIndex();
    
    // Retrieve the Cerenkov Angle Integrals for this material  
    
    G4PhysicsOrderedFreeVector* CerenkovAngleIntegrals1 =
    (G4PhysicsOrderedFreeVector*)((*thePhysicsTable1)(materialIndex));
    G4PhysicsOrderedFreeVector* CerenkovAngleIntegrals2 =
    (G4PhysicsOrderedFreeVector*)((*thePhysicsTable2)(materialIndex));
    
    if(!(CerenkovAngleIntegrals1->IsFilledVectorExist()))return 0.0;
    if(!(CerenkovAngleIntegrals2->IsFilledVectorExist()))return 0.0;
    
    // Min and Max photon energies 
    G4double Pmin = Rindex->GetMinPhotonEnergy();
    G4double Pmax = Rindex->GetMaxPhotonEnergy();
    
    // Min and Max Refraction Indices 
    G4double nMin = Rindex->GetMinProperty();	
    G4double nMax = Rindex->GetMaxProperty();
    
    // Max Cerenkov Angle Integral 
    G4double CAImax1 = CerenkovAngleIntegrals1->GetMaxValue();
    G4double CAImax2 = CerenkovAngleIntegrals2->GetMaxValue();
    
    G4double dp, ge;
    
    // If n(Pmax) < 1/Beta -- no photons generated 
    
    int kind=-1;
    
    if (nMax < BetaInverse) {
        kind=0;

        dp = 0;
        ge = 0;
    } 
    
    // otherwise if n(Pmin) >= 1/Beta -- photons generated  
    
    else if (nMin > BetaInverse) {
        kind=1;
        
        dp = CAImax1; // Pmax-Pmin for unbiased production
        ge = CAImax2; 
    } 
    
    // If n(Pmin) < 1/Beta, and n(Pmax) >= 1/Beta, then
    // we need to find a P such that the value of n(P) == 1/Beta.
    // Interpolation is performed by the GetPhotonEnergy() and
    // GetProperty() methods of the G4MaterialPropertiesTable and
    // the GetValue() method of G4PhysicsVector.  
    
    else {
        kind=2;
        
        Pmin = Rindex->GetPhotonEnergy(BetaInverse);
        dp = Pmax - Pmin;
        
        // need boolean for current implementation of G4PhysicsVector
        // ==> being phased out
        G4bool isOutRange;

        G4double CAImin1 = CerenkovAngleIntegrals1->
        GetValue(Pmin, isOutRange);
        dp = CAImax1 - CAImin1;

        G4double CAImin2 = CerenkovAngleIntegrals2->
        GetValue(Pmin, isOutRange);
        ge = CAImax2 - CAImin2;
        
        if (verboseLevel>0) {
            G4cout << "CAImin2 = " << CAImin2 << G4endl;
            G4cout << "ge = " << ge << G4endl;
        }
    }
    
    // Calculate number of photons 
    G4double NumPhotons = Rfact * charge/eplus * charge/eplus *
    (dp - ge * BetaInverse*BetaInverse);

    //std::cout << "NumPhotons/m=" << NumPhotons/(charge/eplus * charge/eplus)/(1./m) << " @ beta=" << beta << ", kind=" << kind << std::endl;
    //std::cout << "CAImax1=" << CAImax1 << ", CAImax2=" << CAImax2 << ", Pmax=" << Pmax << ", Pmin=" << Pmin << std::endl;
    //std::cout << "Wmax=" << (h_Planck*c_light/Pmin)/nm << ", Wmin=" << (h_Planck*c_light/Pmax)/nm << std::endl;
    
    return NumPhotons;		
}
