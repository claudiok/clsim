#ifndef TRKDETECTORCONSTRUCTION_H_INCLUDED
#define TRKDETECTORCONSTRUCTION_H_INCLUDED
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
 * @file TrkDetectorConstruction.hh
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef TrkDetectorConstruction_H
#define TrkDetectorConstruction_H 1

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Box;
class G4Material;

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
    G4Material*
    MaterialWithSingleIsotope(G4String, G4String, G4double, G4int, G4int);
    
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
    G4Element* D;
    G4Element* O18;
    G4Material* normal_ice;
    G4Material* semiheavy_ice;
    G4Material* H2O18;
    G4Material* iso_ice;
    G4Material* Ice;

    G4MaterialPropertiesTable* MPT;
    
    
    // external medium properties
    I3CLSimMediumPropertiesConstPtr mediumProperties_;
};

#endif

#endif  // TRKDETECTORCONSTRUCTION_H_INCLUDED
