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
// Description: Discrete Process - Generation of Cerenkov Photons
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
#include "clsim/I3CLSimStep.h"

#include <boost/thread.hpp>

#include "G4Version.hh"

#if G4VERSION_NUMBER >= 950
// The G4MaterialPropertyVector is gone since 4.9.5.
// It has been typedef'd to a G4UnorderedPhysicsVector
// with a different interface. Try to support both
// versions with an ifdef.
#define MATERIAL_PROPERTY_VECTOR_IS_PHYSICS_VECTOR

#if G4VERSION_NUMBER >= 1000
using CLHEP::h_Planck;
using CLHEP::c_light;
using CLHEP::eplus;
using CLHEP::m;
using CLHEP::cm;
using CLHEP::um;
using CLHEP::nm;
using CLHEP::ns;
using CLHEP::eV;
#endif
#endif

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
    
    G4MaterialPropertyVector* Rindex = 
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
    {
        I3CLSimStep newStep;
    
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
        newStep.SetSourceType(0); // cherenkov emission
        
        // emit the step, possibly blocking
        eventInformation->StopClock();
        eventInformation->emitStep(newStep);
        eventInformation->StartClock();
    }
    
    // changed Apr 19, 2012: do not automatically assume a particle
    // below the Cherenkov threshold should be killed. It might decay
    // later while we can still se it.
    //
    // kill the particle if it has fallen below the threshold
    //if (1./beta2 > eventInformation->maxRefractiveIndex)
    //{
    //    aParticleChange.ProposeTrackStatus(fStopAndKill);
    //}
    
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
                
#ifdef MATERIAL_PROPERTY_VECTOR_IS_PHYSICS_VECTOR
                G4double currentRI = (*theRefractionIndexVector)[0];
#else
                theRefractionIndexVector->ResetIterator();
                ++(*theRefractionIndexVector);  // advance to 1st entry 
                
                G4double currentRI = theRefractionIndexVector->
                GetProperty();
#endif
                
                if (currentRI > 1.0) {
                    
                    // Create first (photon energy, Cerenkov Integral)
                    // pair  
                    
#ifdef MATERIAL_PROPERTY_VECTOR_IS_PHYSICS_VECTOR
                    G4double currentPM = theRefractionIndexVector->
                                           Energy(0);
#else
                    G4double currentPM = theRefractionIndexVector->
                    GetPhotonEnergy();
#endif
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
                    
#ifdef MATERIAL_PROPERTY_VECTOR_IS_PHYSICS_VECTOR
                    for (size_t i = 1;
                         i < theRefractionIndexVector->GetVectorLength();
                         i++)
                    {
                        currentRI=(*theRefractionIndexVector)[i];
                        currentPM = theRefractionIndexVector->Energy(i);
#else
                    while(++(*theRefractionIndexVector))
                    {
                        currentRI=theRefractionIndexVector->
                        GetProperty();
                        
                        currentPM = theRefractionIndexVector->
                        GetPhotonEnergy();
#endif
                        
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
    
    G4MaterialPropertyVector* Rindex = NULL;
    
    if (aMaterialPropertiesTable)
        Rindex = aMaterialPropertiesTable->GetProperty("RINDEX");
    
    G4double nMax;
    if (Rindex) {
#ifdef MATERIAL_PROPERTY_VECTOR_IS_PHYSICS_VECTOR
        nMax = Rindex->GetMaxValue();
#else
        nMax = Rindex->GetMaxProperty();
#endif
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
                                       G4MaterialPropertyVector* Rindex) const
{
    const G4double Rfact = 369.81/(eV * cm);
    
    if(beta <= 0.0)return 0.0;
    
    G4double BetaInverse = 1./beta;
    
    // Vectors used in computation of Cerenkov Angle Integral:
    //  - Refraction Indices for the current material
    //  - new G4PhysicsOrderedFreeVector allocated to hold CAI's
    
    G4int materialIndex = aMaterial->GetIndex();
    
    // Retrieve the Cerenkov Angle Integrals for this material  
    
    G4PhysicsOrderedFreeVector* CerenkovAngleIntegrals1 =
    (G4PhysicsOrderedFreeVector*)((*thePhysicsTable1)(materialIndex));
    G4PhysicsOrderedFreeVector* CerenkovAngleIntegrals2 =
    (G4PhysicsOrderedFreeVector*)((*thePhysicsTable2)(materialIndex));
    
    if(!(CerenkovAngleIntegrals1->IsFilledVectorExist()))return 0.0;
    if(!(CerenkovAngleIntegrals2->IsFilledVectorExist()))return 0.0;
    
    // Min and Max Refraction Indices 
#ifdef MATERIAL_PROPERTY_VECTOR_IS_PHYSICS_VECTOR
    G4double nMin = Rindex->GetMinValue();
    G4double nMax = Rindex->GetMaxValue();
#else
    G4double nMin = Rindex->GetMinProperty();
    G4double nMax = Rindex->GetMaxProperty();
#endif
    
    // Max Cerenkov Angle Integral 
    G4double CAImax1 = CerenkovAngleIntegrals1->GetMaxValue();
    G4double CAImax2 = CerenkovAngleIntegrals2->GetMaxValue();
    
    G4double dp, ge;
    
    // If n(Pmax) < 1/Beta -- no photons generated 
    
    if (nMax < BetaInverse) {
        dp = 0;
        ge = 0;
    } 
    
    // otherwise if n(Pmin) >= 1/Beta -- photons generated  
    
    else if (nMin > BetaInverse) {
        dp = CAImax1; // Pmax-Pmin for unbiased production
        ge = CAImax2; 
    } 
    
    // If n(Pmin) < 1/Beta, and n(Pmax) >= 1/Beta, then
    // we need to find a P such that the value of n(P) == 1/Beta.
    // Interpolation is performed by the GetPhotonEnergy() and
    // GetProperty() methods of the G4MaterialPropertiesTable and
    // the GetValue() method of G4PhysicsVector.  
    
    else {
#ifdef MATERIAL_PROPERTY_VECTOR_IS_PHYSICS_VECTOR
      G4double Pmin = Rindex->GetEnergy(BetaInverse);
#else
      G4double Pmin = Rindex->GetPhotonEnergy(BetaInverse);
#endif

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

    return NumPhotons;      
}
