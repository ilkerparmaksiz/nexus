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

#ifndef PrimaryGenerator_h
#define PrimaryGenerator_h 1

#include "G4VPrimaryGenerator.hh"
#include "G4IonTable.hh"
#include "G4GenericMessenger.hh"

#include "G4BetaMinusDecay.hh"

#include "FileHandling.hh"


class G4Event;
//class G4DecayProducts;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGenerator : public G4VPrimaryGenerator
{
  public:
    PrimaryGenerator();    
   ~PrimaryGenerator();

  public:
  virtual void GeneratePrimaryVertex(G4Event*) {};
  virtual void GeneratePrimaryVertexIon(G4Event*,std::vector<double> &);
  virtual void GenerateSingleParticle(G4Event*);
  virtual void Generate(G4Event* event, std::vector<double> &xyz); // Choose which generator to launch


private:
    bool fFSNeutrino;
    G4IonTable* fIonTable;
    G4ParticleDefinition *fParentNucleus;
    G4BetaMinusDecay *fBetaMinusChannel; 

     G4GenericMessenger* msg_;

     G4ThreeVector momentum_;


    // For reading in files
    filehandler::FileHandling FileHandler;

    std::vector<std::vector<G4double>> electron_data;
    std::vector<std::vector<G4double>> alpha_data;
    G4String ParticleType_;
    G4ThreeVector Position_;
    G4double energy_;
    G4bool Iso_;
    G4bool useNeedle;
    G4String GeneratorMode_;



};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
