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
/// \file PrimaryGenerator.hh
/// \brief Definition of the PrimaryGenerator class
//
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef GarfieldParticleGenerator_h
#define GarfieldParticleGenerator_h 1

#include <G4VPrimaryGenerator.hh>
#include "G4IonTable.hh"
#include "G4GenericMessenger.hh"

#include "G4BetaMinusDecay.hh"
#include "FileHandling.h"
#include "GeometryBase.h"

class G4GenericMessenger;
class G4Event;
class G4ParticleDefinition;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace nexus{
    class GeometryBase;
    class GarfieldParticleGenerator : public G4VPrimaryGenerator
    {
      public:
        GarfieldParticleGenerator();
       ~GarfieldParticleGenerator();

      public:
        void GeneratePrimaryVertex(G4Event*) ;
        void GeneratePrimaryVertexIon(G4Event*);

        void GenerateSingleParticle(G4Event*);



    private:
        void SetParticleDefinition(G4String);
        const GeometryBase* geom_; ///< Pointer to the detector geometry
        bool fFSNeutrino;
        G4IonTable* fIonTable;
        G4ParticleDefinition* particle_definition_;

        G4ParticleDefinition *fParentNucleus;
        G4BetaMinusDecay *fBetaMinusChannel;

         G4GenericMessenger* msg_;

         G4ThreeVector momentum_;


        // For reading in files
        FileHandling::FileHandling FileHandler;

        std::vector<std::vector<G4double>> electron_data;
        std::vector<std::vector<G4double>> alpha_data;
        G4String ParticleType_;
        G4ThreeVector Position_;
        G4double energy_;
        G4bool Iso_;
        G4bool useNeedle;
        G4String GeneratorMode_;



    };
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
