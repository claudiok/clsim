#ifndef TrkDetectorConstruction_H
#define TrkDetectorConstruction_H 1

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Box;

#include "G4Material.hh"
#include "G4RotationMatrix.hh"
#include "G4VUserDetectorConstruction.hh"

#include "clsim/I3CLSimMediumProperties.h"

#include <map>

class TrkDetectorConstruction : public G4VUserDetectorConstruction
{
public:
    TrkDetectorConstruction(I3CLSimMediumPropertiesConstPtr mediumProperties);
    virtual ~TrkDetectorConstruction();
    
    G4VPhysicalVolume* Construct();
    
    //Functions to modify the geometry
    
    //rebuild the geometry based on changes. must be called
    void UpdateGeometry();
    G4bool GetUpdated() {return updated;}

private:
    void DefineMaterials();
    G4VPhysicalVolume* ConstructDetector();
    
    G4bool updated;
    
    G4Box* world_sol;
    G4LogicalVolume* world_log;
    G4VPhysicalVolume* world_phys;
    
    G4Box* rock_sol;
    G4LogicalVolume* rock_log;
    G4VPhysicalVolume* rock_phys;

    G4Box* air_sol;
    G4LogicalVolume* air_log;
    G4VPhysicalVolume* air_phys;
    
    //Materials & Elements
    G4Element* O;
    G4Element* N;
    G4Element* C;
    G4Material* Air;
    G4Material* Vacuum;
    G4Material* H2O;
    G4Material* StdRock;
    G4Element* H;
    G4Element* Na;
    G4Element* Mg;
    G4Element* Cl;
    G4Element* Ca;
    
    G4MaterialPropertiesTable* MPT;
    
    
    // external medium properties
    I3CLSimMediumPropertiesConstPtr mediumProperties_;
};

#endif
