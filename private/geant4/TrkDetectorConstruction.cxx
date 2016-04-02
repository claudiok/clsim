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
 * @file TrkDetectorConstruction.cxx
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#include "TrkDetectorConstruction.hh"

#include <vector>
#include <map>

#include "G4RunManager.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4MaterialTable.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4UserLimits.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4GeometryManager.hh"

#include "icetray/I3Logging.h"
#include "icetray/I3Units.h"

TrkDetectorConstruction::TrkDetectorConstruction(I3CLSimMediumPropertiesConstPtr mediumProperties)
:
mediumProperties_(mediumProperties)
{
    if (!mediumProperties_)
        log_fatal("MediumProperties are NULL!");
}

TrkDetectorConstruction::~TrkDetectorConstruction()
{
}

void TrkDetectorConstruction::DefineMaterials(){
    G4double a;  // atomic mass
    G4double z;  // atomic number
    G4double density;
    
    //***Elements
    H = new G4Element("Hydrogen",   "H",  z=1.,  a=  1.00794*g/mole);
    C = new G4Element("C",          "C",  z=6.,  a= 12.01   *g/mole);
    N = new G4Element("Nitrogen",   "N",  z=7.,  a= 14.0067 *g/mole);
    O = new G4Element("Oxygen",     "O",  z=8.,  a= 15.9994 *g/mole);
    
    Na = new G4Element("Sodium",    "Na", z=11., a= 22.98976928*g/mole);
    Mg = new G4Element("Magnesium", "Mg", z=12., a= 24.3050 *g/mole);
    Cl = new G4Element("Chlorine",  "Cl", z=17., a= 35.4527 *g/mole);
    Ca = new G4Element("Calcium",   "Ca", z=20., a= 40.0784*g/mole);
    
    //***Materials
    //Vacuum
    Vacuum = new G4Material("Vacuum",z=1.,a=1.01*g/mole,
                            density=universe_mean_density,kStateGas,0.1*kelvin,
                            1.e-19*pascal); 

    //Standard Rock (approximation with a composition of individual elements.
    //Geant4 does not play well with average atomic weights, etc.)
    //This is not too different from Gran Sasso rock, which
    //in turn should not be too different from standard rock. The density of 2.65g/cm^3
    //is the original "standard rock" density given by the full PDG report.
    StdRock = new G4Material("StdRock",density=2.65*g/cm3,4, kStateSolid );
    StdRock->AddElement(O,  52.*perCent);
    StdRock->AddElement(Ca, 27.*perCent);
    StdRock->AddElement(C,  12.*perCent);
    StdRock->AddElement(Mg,  9.*perCent);

    //Air
    Air = new G4Material("Air", density= 1.29*mg/cm3, 2);
    Air->AddElement(N, 70*perCent);
    Air->AddElement(O, 30*perCent);

    // Sea Water (with salt)
    const double mediumDensity = (mediumProperties_->GetMediumDensity()/( I3Units::g/I3Units::cm3 )) * (g/cm3);
    
    H2O = new G4Material("H2O", density=mediumDensity, 5, kStateLiquid,286.35*kelvin,1.*atmosphere);
    H2O->AddElement(H,  10.74*perCent);
    H2O->AddElement(O,  85.41*perCent);
    H2O->AddElement(Na,  1.32*perCent);
    H2O->AddElement(Mg,  0.15*perCent);
    H2O->AddElement(Cl,  2.37*perCent);
    
    
    //***Material properties tables
    
    // refractive index. We get it as a function, convert it to a table.
    
    // the refractive index has to be constant for all layers at the moment.
    // TODO: improve this, implement layers in Geant4!
    {
        const uint32_t numLayers = mediumProperties_->GetLayersNum();
        if (numLayers==0) log_fatal("numLayers==0");
        
        I3CLSimFunctionConstPtr firstLayer = mediumProperties_->GetPhaseRefractiveIndex(0);
        
        for (uint32_t i=1;i<numLayers;++i)
        {
            I3CLSimFunctionConstPtr thisLayer = mediumProperties_->GetPhaseRefractiveIndex(i);
            
            if (*firstLayer != *thisLayer)
                log_fatal("The current version of the Geant4 detector constructor can only handle media with un-layered refractie indices..");
        }
        
    }
    
    // fill the refractive index map for a fixed number of wavelengths
    I3CLSimFunctionConstPtr refractiveIndexFunc = mediumProperties_->GetPhaseRefractiveIndex(0);
    
    const unsigned int initialPoints=20;
    const double fromWlen = mediumProperties_->GetMinWavelength();
    const double toWlen = mediumProperties_->GetMaxWavelength();
    
    // sanity check
    if ((fromWlen < refractiveIndexFunc->GetMinWlen()) ||
        (toWlen > refractiveIndexFunc->GetMaxWlen()))
        log_fatal("Global medium wavelength range is larger than the refractive index wavelength range!");
    
    std::map<double, double> wlenToRindexMap;
    
    //G4cout << "map<wlen,rindex> (before):" << G4endl;
    for (unsigned int i=0;i<initialPoints;++i)
    {
        const double wlen = fromWlen + static_cast<double>(i)*(toWlen-fromWlen)/static_cast<double>(initialPoints-1);
        const double rindex = refractiveIndexFunc->GetValue(wlen);
        
        wlenToRindexMap.insert(std::make_pair(wlen, rindex));
        
        //G4cout << " " << wlen/I3Units::nanometer << "nm -> " << rindex << G4endl;
    }
    //G4cout << G4endl;
    
    // check if any intermediate points are needed
    const double maxRelativeError = 1e-5; // no more than 1% deviation of interpolated result from real function

    for (;;)
    {
        bool doItAgain=false;

        for (std::map<double, double>::iterator it=wlenToRindexMap.begin();
             it!=wlenToRindexMap.end();++it)
        {
            std::map<double, double>::iterator next_it=it; ++next_it;
            if (next_it==wlenToRindexMap.end()) break;
            
            const double inbetween_wlen = (next_it->first + it->first)/2.;
            const double interpolated_rindex = (next_it->second + it->second)/2.;
            const double real_rindex = refractiveIndexFunc->GetValue(inbetween_wlen);
            const double relativeError = (real_rindex-interpolated_rindex)/real_rindex;
            
            // do we need to insert a new data point?
            if (fabs(relativeError) > maxRelativeError)
            {
                //G4cout << "relative error @ " << inbetween_wlen/I3Units::nanometer << "nm: " << relativeError << ", interp=" << interpolated_rindex << ", real=" << real_rindex << G4endl;

                wlenToRindexMap.insert(std::make_pair(inbetween_wlen, real_rindex));
                
                // start from scratch
                doItAgain=true;
                break;
            }
        }
        
        if (!doItAgain) break;
    }

    std::vector<double> rindex;
    std::vector<double> rindex_ppckov;

    //G4cout << "map<wlen,rindex>:" << G4endl;
    for (std::map<double, double>::reverse_iterator it=wlenToRindexMap.rbegin();
         it!=wlenToRindexMap.rend();++it)
    {
        //G4cout << " " << it->first/I3Units::nanometer << "nm -> " << it->second << G4endl;
        
        const double wlenInG4Units = (it->first/I3Units::nanometer)*nanometer;
        
        rindex_ppckov.push_back((h_Planck*c_light)/wlenInG4Units);
        rindex.push_back(it->second);
    }
    G4cout << G4endl;

    
    // set up material properties table
    MPT = new G4MaterialPropertiesTable();
    MPT->AddProperty("RINDEX", &(rindex_ppckov[0]), &(rindex[0]), rindex_ppckov.size());
    
    H2O->SetMaterialPropertiesTable(MPT);
}

G4VPhysicalVolume* TrkDetectorConstruction::Construct()
{
    DefineMaterials();
    return ConstructDetector();
}

G4VPhysicalVolume* TrkDetectorConstruction::ConstructDetector()
{
    G4cout << " ===---***---=== CONSTRUCTING DETECTOR" << G4endl;
    
    const double seaFloorZCoordinate = (mediumProperties_->GetRockZCoord()/I3Units::m)*m;
    const double airZCoordinate = (mediumProperties_->GetAirZCoord()/I3Units::m)*m;

    // Set up the world volume. It is a box filled with water.
    const double simvolume_x = 5.*km; // a rather arbitrary safety margin
    const double simvolume_y = 5.*km; // a rather arbitrary safety margin
    // add 1km of rock and 1km of air // this, together with the simulation volume should
    // be enough for high-energy input that ran through MMC/PROPOSAL
    const double simvolume_z = std::max(fabs(seaFloorZCoordinate)+1.*km, fabs(airZCoordinate)+1.*km);
    
    
    G4cout << "world_r_x=+/-" << simvolume_x/km << "km, " 
           << "world_r_y=+/-" << simvolume_y/km << "km, " 
           << "world_r_z=+/-" << simvolume_z/km << "km"
           << G4endl; 
    
    world_sol = new G4Box("world_sol",
                          simvolume_x, // half-lengths!
                          simvolume_y,
                          simvolume_z);
    world_log = new G4LogicalVolume(world_sol,H2O,
                                    "world_log"
                                    //,0,0,0,false
                                    );
    
    world_phys = new G4PVPlacement(0,G4ThreeVector(),
                                  world_log,"world_phys",NULL,false,0,false);

    
    
    G4cout << "Placing rock at z=" << seaFloorZCoordinate/m << "m." << G4endl;

    // place the sea floor rock box
    {
        const double seafloor_rock_height = simvolume_z+seaFloorZCoordinate;
        rock_sol = new G4Box("rock_sol",simvolume_x, simvolume_y, seafloor_rock_height/2.);
        rock_log = new G4LogicalVolume(rock_sol,StdRock,"rock_log");
        G4RotationMatrix* rot = new G4RotationMatrix();
        rock_phys = new G4PVPlacement(rot,G4ThreeVector(0.,0.,-simvolume_z+seafloor_rock_height/2.),rock_log,"rock_phys",world_log,false,0);
    }
    
    G4cout << "Placing air at z=" << airZCoordinate/m << "m." << G4endl;

    // place the air box
    {
        const double air_height = simvolume_z-airZCoordinate;
        air_sol = new G4Box("air_sol",simvolume_x, simvolume_y, air_height/2.);
        air_log = new G4LogicalVolume(air_sol,Air,"air_log");
        G4RotationMatrix* rot = new G4RotationMatrix();
        air_phys = new G4PVPlacement(rot,G4ThreeVector(0.,0.,simvolume_z-air_height/2.),air_log,"air_phys",world_log,false,0);
    }    
    
    
    
    
    return world_phys;
}

void TrkDetectorConstruction::UpdateGeometry(){
    
    // clean-up previous geometry
    G4GeometryManager::GetInstance()->OpenGeometry();
    
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();
    G4LogicalSkinSurface::CleanSurfaceTable();
    G4LogicalBorderSurface::CleanSurfaceTable();
    
    //define new one
    G4RunManager::GetRunManager()->DefineWorldVolume(ConstructDetector());
    G4RunManager::GetRunManager()->GeometryHasBeenModified();

    updated=false;
}



