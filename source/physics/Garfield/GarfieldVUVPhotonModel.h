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
#include "IOUtils.h"
#include "G4GenericMessenger.hh"

#include "G4VFastSimulationModel.hh"
#include "Garfield/Medium.hh"
#include "IonizationSD.h"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/ComponentUser.hh"
#include "GarfieldHelper.h"
#include "GeometryBase.h"
#include "CRAB0.h"
#include "IonizationSD.h"

class G4GenericMessenger;
namespace nexus {
    class GarfieldVUVPhotonModel : public G4VFastSimulationModel
    {
    public:
      //-------------------------
      // Constructor, destructor
      //-------------------------
        GarfieldVUVPhotonModel(G4String modelName, G4Region* envelope, GarfieldHelper GH,  IonizationSD* ionisd);
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

        G4String gasFile;
        G4String ionMobFile;

        Garfield::MediumMagboltz* fMediumMagboltz;
        Garfield::AvalancheMC* fAvalancheMC;
        Garfield::Sensor* fSensor;

        std::vector<uint> counter {0,0,0,0};

        // Variable to store the EL timing profiles to sample from
        // <event> <photon> <x,y,z,t of photon>
        std::vector<std::vector<std::vector<G4double>>> EL_profiles;

        // Vector consisting of the event numbers from the simulated Garfield using
        // COMSOL geometry
        std::vector<G4double> EL_events;

        GarfieldHelper GH_;


        G4double ELPos_; // cm
        G4double FCTop_; // cm

        /// The sensitive detector to fill hits into
        IonizationSD* fGarfieldSD;

    };
    void userHandle(double x, double y, double z, double t, int type, int level,Garfield::Medium * m);

}
#endif /* GarfieldVUVPhotonModel_H_ */
