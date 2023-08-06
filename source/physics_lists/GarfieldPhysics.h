//
// Created by ilker on 8/6/23.
//

#ifndef NEXUS_GARFIELDPHYSICS_H
#define NEXUS_GARFIELDPHYSICS_H

#include "G4VModularPhysicsList.hh"
#include "G4SystemOfUnits.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4GenericIon.hh"

#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
// particles
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

#include "G4/NESTProc.hh"
class G4GenericMessenger;
class PhysicsListMessenger;
class DetectorConstruction;
class G4FastSimulationPhysics;

namespace nexus {

    class GarfieldPhysics: public G4VPhysicsConstructor
    {
    public:
        /// Constructor
        GarfieldPhysics();
        /// Destructor
        ~GarfieldPhysics();

        /// Construct all required particles (Geant4 mandatory method)
        virtual void ConstructParticle();
        /// Construct all required physics processes (Geant4 mandatory method)
        virtual void ConstructProcess();

        void AddParametrisation();



    private:
        G4GenericMessenger* msg_;
        void AddIonGasModels();
        bool IonGasModels_;



        PhysicsListMessenger* pMessenger;
        G4FastSimulationPhysics* fastSimulationPhysics;
    };

} // end namespace nexus
#endif //NEXUS_GARFIELDPHYSICS_H
