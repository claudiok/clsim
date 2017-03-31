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
 * @file I3CLSimI3ParticleGeantConverter.hh
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMI3PARTCILEGEANTCONVERTER_H_INCLUDED
#define I3CLSIMI3PARTCILEGEANTCONVERTER_H_INCLUDED

/**
 * Helper function: Converts I3Particles from/to G4Particles.
 *
 * It either takes an G4Track object and sets the properties
 * of an I3Particle to correspond to it takes an I3Particle
 * and sets the properties of an G4ParticleGun to shoot paticles
 * of this type.
 */

#include <string>
#include <cmath>

#include "G4ParticleGun.hh"
#include "G4Track.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "G4Version.hh"

namespace I3CLSimI3ParticleGeantConverter
{
#ifndef I3PARTICLE_SUPPORTS_PDG_ENCODINGS
    inline I3Particle::ParticleType ConvertPDGEncodingToI3ParticleType(int pdg_code)
    {
        switch (pdg_code)
        {
            //case -3122:  return I3Particle::Hadrons; // DeltaBar
            case -2212:  return I3Particle::PMinus;
            //case -2112:  return I3Particle::Hadrons; // NeutronBar
            case -321:   return I3Particle::KMinus;
            case -311:   return I3Particle::K0_Short; // not really right.. (K0)
            case -211:   return I3Particle::PiMinus;
            case -16:    return I3Particle::NuTauBar;
            case -15:    return I3Particle::TauPlus;
            case -14:    return I3Particle::NuMuBar;
            case -13:    return I3Particle::MuPlus;
            case -12:    return I3Particle::NuEBar;
            case -11:    return I3Particle::EPlus;
            case 0:      return I3Particle::unknown; // invalid
            case 11:     return I3Particle::EMinus;
            case 12:     return I3Particle::NuE;
            case 13:     return I3Particle::MuMinus;
            case 14:     return I3Particle::NuMu;
            case 15:     return I3Particle::TauMinus;
            case 16:     return I3Particle::NuTau;
            case 22:     return I3Particle::Gamma; 
            case 111:    return I3Particle::Pi0;
            case 130:    return I3Particle::K0_Long;
            case 211:    return I3Particle::PiPlus;
            case 310:    return I3Particle::K0_Short;
            case 311:    return I3Particle::K0_Long; // not really right.. (K0Bar)
            case 321:    return I3Particle::KPlus;
            //case 411:    return I3Particle::Hadrons; // D+
            //case 421:    return I3Particle::Hadrons; // D0
            //case 431:    return I3Particle::Hadrons; // (something else.. we should store PDG codes directly..)
            case 2112:   return I3Particle::Neutron;
            case 2212:   return I3Particle::PPlus;
            //case 3112:   return I3Particle::Hadrons; // Sigma-
            //case 3122:   return I3Particle::Hadrons; // Delta
            //case 3212:   return I3Particle::Hadrons; // (something else.. we should store PDG codes directly..)
            //case 3222:   return I3Particle::Hadrons; // Sigma+
            //case 4122:   return I3Particle::Hadrons; // (something else.. we should store PDG codes directly..)
            //case 4212:   return I3Particle::Hadrons; // (something else.. we should store PDG codes directly..)
            //case 4222:   return I3Particle::Hadrons; // (something else.. we should store PDG codes directly..)
                
            case 1000020040:  return I3Particle::He4Nucleus;
            case 1000030070:  return I3Particle::Li7Nucleus;
            case 1000040090:  return I3Particle::Be9Nucleus;
            case 1000050110:  return I3Particle::B11Nucleus;
            case 1000060120:  return I3Particle::C12Nucleus;
            case 1000070140:  return I3Particle::N14Nucleus;
            case 1000080160:  return I3Particle::O16Nucleus;
            case 1000090190:  return I3Particle::F19Nucleus;
            case 1000100200:  return I3Particle::Ne20Nucleus;
            case 1000110230:  return I3Particle::Na23Nucleus;
            case 1000120240:  return I3Particle::Mg24Nucleus;
            case 1000130270:  return I3Particle::Al27Nucleus;
            case 1000140280:  return I3Particle::Si28Nucleus;
            case 1000150310:  return I3Particle::P31Nucleus;
            case 1000160320:  return I3Particle::S32Nucleus;
            case 1000170350:  return I3Particle::Cl35Nucleus;
            case 1000180400:  return I3Particle::Ar40Nucleus;
            case 1000190390:  return I3Particle::K39Nucleus;
            case 1000200400:  return I3Particle::Ca40Nucleus;
            case 1000210450:  return I3Particle::Sc45Nucleus;
            case 1000220480:  return I3Particle::Ti48Nucleus;
            case 1000230510:  return I3Particle::V51Nucleus;
            case 1000240520:  return I3Particle::Cr52Nucleus;
            case 1000250550:  return I3Particle::Mn55Nucleus;
            case 1000260560:  return I3Particle::Fe56Nucleus;
                
            default: return I3Particle::unknown; // invalid/unknown
        }
    }

    
    inline int ConvertI3ParticleTypeToPDGEncoding(I3Particle::ParticleType type)
    {
        switch (type)
        {
            case I3Particle::unknown: return 0; // invalid
            case I3Particle::Gamma: return 22;
            case I3Particle::EPlus: return -11;
            case I3Particle::EMinus: return 11;
            case I3Particle::MuPlus: return -13;
            case I3Particle::MuMinus: return 13;
            case I3Particle::Pi0: return 111;
            case I3Particle::PiPlus: return 211;
            case I3Particle::PiMinus: return -211;
            case I3Particle::K0_Long: return 130;
            case I3Particle::KPlus: return 321;
            case I3Particle::KMinus: return -321;
            case I3Particle::Neutron: return 2112;
            case I3Particle::PPlus: return 2212;
            case I3Particle::PMinus: return -2212;
            case I3Particle::K0_Short: return 310;
            case I3Particle::NuE: return 12;
            case I3Particle::NuEBar: return -12;
            case I3Particle::NuMu: return 14;
            case I3Particle::NuMuBar: return -14;
            case I3Particle::TauPlus: return -15;
            case I3Particle::TauMinus: return 15;
            case I3Particle::NuTau: return 16;
            case I3Particle::NuTauBar: return -16;
                
            case I3Particle::He4Nucleus: return 1000020040;
            case I3Particle::Li7Nucleus: return 1000030070;
            case I3Particle::Be9Nucleus: return 1000040090;
            case I3Particle::B11Nucleus: return 1000050110;
            case I3Particle::C12Nucleus: return 1000060120;
            case I3Particle::N14Nucleus: return 1000070140;
            case I3Particle::O16Nucleus: return 1000080160;
            case I3Particle::F19Nucleus: return 1000090190;
            case I3Particle::Ne20Nucleus: return 1000100200;
            case I3Particle::Na23Nucleus: return 1000110230;
            case I3Particle::Mg24Nucleus: return 1000120240;
            case I3Particle::Al27Nucleus: return 1000130270;
            case I3Particle::Si28Nucleus: return 1000140280;
            case I3Particle::P31Nucleus: return 1000150310;
            case I3Particle::S32Nucleus: return 1000160320;
            case I3Particle::Cl35Nucleus: return 1000170350;
            case I3Particle::Ar40Nucleus: return 1000180400;
            case I3Particle::K39Nucleus: return 1000190390;
            case I3Particle::Ca40Nucleus: return 1000200400;
            case I3Particle::Sc45Nucleus: return 1000210450;
            case I3Particle::Ti48Nucleus: return 1000220480;
            case I3Particle::V51Nucleus: return 1000230510;
            case I3Particle::Cr52Nucleus: return 1000240520;
            case I3Particle::Mn55Nucleus: return 1000250550;
            case I3Particle::Fe56Nucleus: return 1000260560;
            case I3Particle::CherenkovPhoton: return 0; // invalid
                
            case I3Particle::Nu: return 0;          // invalid
            case I3Particle::Monopole: return 0;    //invalid
            case I3Particle::Brems: return 11;      // NOTE: check these
            case I3Particle::DeltaE: return 11;     // NOTE: check these
            case I3Particle::PairProd: return 22;   // NOTE: check these
            case I3Particle::NuclInt: return 211;   // we just assume these are pions
            case I3Particle::MuPair: return 22;     // NOTE: check these
            case I3Particle::Hadrons: return 211;   // we just assume these are pions
            case I3Particle::FiberLaser: return 0;  // invalid
            case I3Particle::N2Laser: return 0;     // invalid
            case I3Particle::YAGLaser: return 0;    // invalid
            case I3Particle::STauPlus: return 0;    // invalid
            case I3Particle::STauMinus: return 0;   // invalid
            default: return 0;          // invalid
        }
    }
#endif
    
    
    /**
     * Configures a G4ParticleGun to shoot particles as configured
     * by an I3Particle. Returns "false" if the I3Particle does not contain
     * enough information (i.e., some of the required fields are NaN).
     */
    inline bool SetParticleGun(G4ParticleGun *particleGun, const I3Particle &particle)
    {
        // check the I3Particle
        
        if (std::isnan(particle.GetEnergy())) {
            log_debug("Particle has no energy. Cannot shoot.");
            return false;
        }

        if (std::isnan(particle.GetZenith())) {
            log_debug("Particle has no zenith direction. Cannot shoot.");
            return false;
        }

        if (std::isnan(particle.GetAzimuth())) {
            log_debug("Particle has no zenith direction. Cannot shoot.");
            return false;
        }

        if (std::isnan(particle.GetX()) || std::isnan(particle.GetY()) || std::isnan(particle.GetZ())) {
            log_debug("Particle has no position. Cannot shoot.");
            return false;
        }

        if (std::isnan(particle.GetTime())) {
            log_debug("Particle has no time. Cannot shoot.");
            return false;
        }
        
        if (particle.GetType() == I3Particle::unknown)
        {
            log_debug("Particle has no valid type. Cannot shoot.");
            return false;
        }
        
#ifndef I3PARTICLE_SUPPORTS_PDG_ENCODINGS
        int pdgcode = ConvertI3ParticleTypeToPDGEncoding(particle.GetType());
        if (pdgcode==0) 
#else
        int pdgcode = particle.GetPdgEncoding();
        
        // simple replacements for a few special codes
        if      (pdgcode==-2000001001) pdgcode=11;  // Brems    -> EMinus // NOTE: check these
        else if (pdgcode==-2000001002) pdgcode=11;  // DeltaE   -> EMinus // NOTE: check these
        else if (pdgcode==-2000001003) pdgcode=22;  // PairProd -> Gamma  // NOTE: check these
        else if (pdgcode==-2000001004) pdgcode=211; // NuclInt  -> PiPlus (just assume this is a pion..)
        else if (pdgcode==-2000001005) pdgcode=22;  // MuPair   -> EMinus // NOTE: check these
        else if (pdgcode==-2000001006) pdgcode=211; // Hadrons  -> PiPlus (just assume this is a pion..)
        
        if ((pdgcode==0) || (pdgcode <= -2000000000) || (pdgcode >= 2000000000)) // no special IceTray codes allowed in Geant4..
#endif
        {
            log_warn("Particle with type \"%s\" is not supported by Geant4. Ignoring.",
                     particle.GetTypeString().c_str());
            return false;
        }
        
        // look it up in the particle table
        G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
        G4ParticleDefinition* particleDef = particleTable->FindParticle(pdgcode);
        
        if (!particleDef)
        {
            log_warn("Particle with type I3Particle type \"%s\" (pdgcode=%i) was not found in the Geant4 particle table. Ignoring.",
                     particle.GetTypeString().c_str(), pdgcode);
            return false;
        }
        
        // configure the particle gun
        particleGun->SetParticleDefinition(particleDef);
        
        // configure energy, position, momentum and time. Convert from I3Units to G4Units
#if G4VERSION_NUMBER >= 1000
        particleGun->SetParticleEnergy((particle.GetEnergy()/I3Units::GeV)*CLHEP::GeV);
        particleGun->SetParticlePosition(G4ThreeVector(particle.GetX()/I3Units::m, particle.GetY()/I3Units::m, particle.GetZ()/I3Units::m)*CLHEP::m);
        particleGun->SetParticleTime((particle.GetTime()/I3Units::ns)*CLHEP::ns);
#else
        particleGun->SetParticleEnergy((particle.GetEnergy()/I3Units::GeV)*GeV);
        particleGun->SetParticlePosition(G4ThreeVector(particle.GetX()/I3Units::m, particle.GetY()/I3Units::m, particle.GetZ()/I3Units::m)*m);
        particleGun->SetParticleTime((particle.GetTime()/I3Units::ns)*ns);
#endif
        {
            const I3Direction &dir=particle.GetDir();
            particleGun->SetParticleMomentumDirection(G4ThreeVector(dir.GetX(), dir.GetY(), dir.GetZ()));
        }

        
        particleGun->SetNumberOfParticles(1);
        
        return true;
    }
    
    
    
    
};

#endif //I3CLSIMI3PARTCILEGEANTCONVERTER_H_INCLUDED
