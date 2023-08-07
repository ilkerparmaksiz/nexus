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
// $Id: GarfieldFastSimulationModel.hh 9999994 2015-12-11 14:47:43Z dpfeiffe $
//
/// \file GarfieldFastSimulationModel.hh
/// \brief Definition of the GarfieldFastSimulationModel class

#ifndef GarfieldVUVPhotonModel_h
#define GarfieldVUVPhotonModel_h 1

#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpWLS.hh"
#include "G4OpBoundaryProcess.hh"
#include "FileHandling.hh"

#include "G4VFastSimulationModel.hh"
#include "Medium.hh"
#include "GasBoxSD.hh"
#include "MediumMagboltz.hh"
#include "AvalancheMicroscopic.hh"
#include "AvalancheMC.hh"
#include "ComponentUser.hh"

#include "TrackHeed.hh"

class GasModelParameters;
class DetectorConstruction;
class GasBoxSD;
namespace nexus {
    class GarfieldVUVPhotonModel : public G4VFastSimulationModel
    {
    public:
      //-------------------------
      // Constructor, destructor
      //-------------------------
        GarfieldVUVPhotonModel(GasModelParameters*, G4String, G4Region*,DetectorConstruction*,GasBoxSD*);
        ~GarfieldVUVPhotonModel (){};

        //void SetPhysics(degradPhysics* fdegradPhysics);
        //void WriteGeometryToGDML(G4VPhysicalVolume* physicalVolume);

        virtual G4bool IsApplicable(const G4ParticleDefinition&);
        virtual G4bool ModelTrigger(const G4FastTrack &);
        virtual void DoIt(const G4FastTrack&, G4FastStep&);
        void GenerateVUVPhotons(const G4FastTrack& fastTrack, G4FastStep& fastStep,G4ThreeVector garfPos,G4double garfTime);
            void Reset();
        G4ThreeVector garfPos;
        G4double garfTime;

        // Function to create simple geometry for field
        Garfield::ComponentUser* CreateSimpleGeometry();

        // Generate EL photons in the gap according to Garfield Microphysical Model
        void MakeELPhotonsFromFile(G4FastStep& fastStep, G4double xi, G4double yi, G4double zi, G4double ti);

        // Generate EL photons in the gap according to a simple model
        void MakeELPhotonsSimple(G4FastStep& fastStep, G4double xi, G4double yi, G4double zi, G4double ti);


    private:

        void InitialisePhysics();
        void S1Fill(const G4FastTrack& );

        G4String gasFile;
        G4String ionMobFile;

        DetectorConstruction* detCon;
        G4ThreeVector myPoint;
        G4double time;
        G4double thermalE;
        //degradPhysics* fdegradPhysics;

        Garfield::MediumMagboltz* fMediumMagboltz;
        Garfield::AvalancheMicroscopic* fAvalanche;
        Garfield::AvalancheMC* fAvalancheMC;

        Garfield::Sensor* fSensor;

        GasBoxSD* fGasBoxSD;
        Garfield::TrackHeed* fTrackHeed;
        std::vector<uint> counter {0,0,0,0};

        // Variable to store the EL timing profiles to sample from
        // <event> <photon> <x,y,z,t of photon>
        std::vector<std::vector<std::vector<G4double>>> EL_profiles;

        // Vector consisting of the event numbers from the simulated Garfield using
        // COMSOL geometry
        std::vector<G4double> EL_events;

        // For reading in files
        filehandler::FileHandling FileHandler;

        GasModelParameters* fGasModelParameters;


    };
    void userHandle(double x, double y, double z, double t, int type, int level,Garfield::Medium * m);

}
#endif /* GarfieldVUVPhotonModel_H_ */
