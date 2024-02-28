// ----------------------------------------------------------------------------
// nexus | NexusPhysics.cc
//
// This class registers any new physics process defined in nexus.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "NexusPhysics.h"

#include "IonizationElectron.h"
#include "IonizationClustering.h"
#include "IonizationDrift.h"
#include "Electroluminescence.h"
#include "WavelengthShifting.h"
#include "OpPhotoelectricEffect.h"

// Krishan: this needs cleaning up
#include "G4OpAbsorption.hh"
#include "G4OpBoundaryProcess.hh"
#include "gasNESTdet.h"
#include "NEST/G4/NESTProc.hh"
#include "NEST/G4/NESTS1Photon.hh"

#include <G4GenericMessenger.hh>
#include <G4OpticalPhoton.hh>
#include <G4ProcessManager.hh>
#include <G4ProcessTable.hh>
#include <G4StepLimiter.hh>
#include <G4FastSimulationManagerProcess.hh>
#include <G4PhysicsConstructorFactory.hh>

#include "G4FastSimulationPhysics.hh"
namespace nexus {

  /// Macro that allows the use of this physics constructor
  /// with the generic physics list
  G4_DECLARE_PHYSCONSTR_FACTORY(NexusPhysics);



  NexusPhysics::NexusPhysics():
    G4VPhysicsConstructor("NexusPhysics"),
    clustering_(true), drift_(true), electroluminescence_(true), photoelectric_(false), fastsim_(false), nest_(false)
  {
    msg_ = new G4GenericMessenger(this, "/PhysicsList/Nexus/",
      "Control commands of the nexus physics list.");

    msg_->DeclareProperty("clustering", clustering_,
      "Switch on/off the ionization clustering");

    msg_->DeclareProperty("drift", drift_,
      "Switch on/off the ionization drift.");

    msg_->DeclareProperty("electroluminescence", electroluminescence_,
      "Switch on/off the electroluminescence.");

    msg_->DeclareProperty("photoelectric", photoelectric_,
      "Switch on/off the photoelectric effect.");

    msg_->DeclareProperty("fastsim", fastsim_,
      "Switch on/off NEST/Garfield/Degrad physics.");

    msg_->DeclareProperty("nest", nest_,
      "Switch on/off NEST physics.");

  }



  NexusPhysics::~NexusPhysics()
  {
    delete msg_;
  }



  void NexusPhysics::ConstructParticle()
  {
    IonizationElectron::Definition();
    G4OpticalPhoton::Definition();
    // G4OpticalPhoton::OpticalPhotonDefinition();
    NESTS1Photon::Definition();
  }



  void NexusPhysics::ConstructProcess()
  {
    G4ProcessManager* pmanager = 0;

    // Add our own wavelength shifting process for the optical photon
    pmanager = G4OpticalPhoton::Definition()->GetProcessManager();
    if (!pmanager) {
      G4Exception("[NexusPhysics]", "ConstructProcess()", FatalException,
        "G4OpticalPhoton without a process manager.");
    }
    WavelengthShifting* wls = new WavelengthShifting();
    pmanager->AddDiscreteProcess(wls);

    pmanager = IonizationElectron::Definition()->GetProcessManager();
    if (!pmanager) {
      G4Exception("[NexusPhysics]", "ConstructProcess()", FatalException,
        "Ionization electron without a process manager.");
    }

    // Add drift and electroluminescence to the process table of the ie-

    if (drift_) {
      // First, we remove the standard transportation from the
      // process table of the ionization electron
      G4VProcess* transportation = G4ProcessTable::GetProcessTable()->
        FindProcess("Transportation", IonizationElectron::Definition());
      pmanager->RemoveProcess(transportation);

      IonizationDrift* drift = new IonizationDrift();
      pmanager->AddContinuousProcess(drift);
      pmanager->AddDiscreteProcess(drift);
    }

    if (electroluminescence_) {
      Electroluminescence* el = new Electroluminescence();
      pmanager->AddDiscreteProcess(el);
    }


    // Add clustering to all pertinent particles

    if (clustering_) {

      IonizationClustering* clust = new IonizationClustering();

      auto aParticleIterator = GetParticleIterator();
      aParticleIterator->reset();
      while ((*aParticleIterator)()) {
        G4ParticleDefinition* particle = aParticleIterator->value();
        pmanager = particle->GetProcessManager();

        if (clust->IsApplicable(*particle)) {
          pmanager->AddDiscreteProcess(clust);
          pmanager->AddRestProcess(clust);
        }
      }
    }

    // Add photoelectric effect to optical photons

    if (photoelectric_) {
      OpPhotoelectricEffect* photoe = new OpPhotoelectricEffect();

      auto aParticleIterator = GetParticleIterator();
      aParticleIterator->reset();
      while ((*aParticleIterator)()) {
        G4ParticleDefinition* particle = aParticleIterator->value();

        if (photoe->IsApplicable(*particle)){
          pmanager = particle->GetProcessManager();
          pmanager->AddDiscreteProcess(photoe);
        }
      }
    }

    // Use NEST/Garfield/Degrad physics

    if (fastsim_) {
      
      // This is needed to notify Geant4 that the G4FastSimulationModel is to be used as a possible physics process
      auto fastSimProcess_garfield = new G4FastSimulationManagerProcess("fastSimPhys");

      gasNESTdet* gndet        = new gasNESTdet();
      NEST::NESTcalc* calcNEST = new NEST::NESTcalc(gndet);  

      NEST::NESTProc* theNEST2ScintillationProcess = new NEST::NESTProc("S1",fElectromagnetic, calcNEST, gndet);
      theNEST2ScintillationProcess->SetDetailedSecondaries(true);
      theNEST2ScintillationProcess->SetStackElectrons(true);

      G4OpBoundaryProcess* fBoundaryProcess = new G4OpBoundaryProcess();
      G4OpAbsorption* fAbsorptionProcess    = new G4OpAbsorption();

       // Also switch on clustering for the G4 generated electrons
      IonizationClustering* clust = new IonizationClustering();

      auto aParticleIterator = GetParticleIterator();
      aParticleIterator->reset();

      while ((*aParticleIterator)()) {
        G4ParticleDefinition* particle     = aParticleIterator->value();
        G4ProcessManager*     pmanager     = particle->GetProcessManager();
        G4String              particleName = particle->GetParticleName();

        if (pmanager) {

          // Add the fast sim process for the other particles
          pmanager->AddDiscreteProcess(fastSimProcess_garfield);

          // Add NEST physics
          if (theNEST2ScintillationProcess->IsApplicable(*particle)) {
            // std::cout << "PhysicsList::InitialisePhysics(): particleName, pmanager  " << particleName << ", " << pmanager << "." << std::endl;
            
            // Ony use NEST for alphas
            if (nest_ && particleName == "alpha")
              pmanager->AddProcess(theNEST2ScintillationProcess, ordDefault + 1, ordInActive, ordDefault + 1);
          }

          if (particleName == "opticalphoton") {
            G4cout << " AddDiscreteProcess to OpticalPhoton " << G4endl;
            pmanager->AddDiscreteProcess(fAbsorptionProcess);
            pmanager->AddDiscreteProcess(fBoundaryProcess);
          }

          // Clustering
          if (clust->IsApplicable(*particle) && particleName != "alpha") {
            pmanager->AddDiscreteProcess(clust);
            pmanager->AddRestProcess(clust);
          }

        }
      }

      // Manually add the fastsim physics for ie- created in NEST
      G4ParticleTable*      pTable     = G4ParticleTable::GetParticleTable();
      G4ParticleDefinition* anInstance = pTable->FindParticle("ie-");
      G4ProcessManager*     pmanager   = anInstance->GetProcessManager();
      pmanager->AddDiscreteProcess(fastSimProcess_garfield);

      if (nest_){
        // S1 photons created in NEST
        anInstance = pTable->FindParticle("NESTS1Photon");
        pmanager = anInstance->GetProcessManager();
        pmanager->AddDiscreteProcess(fAbsorptionProcess);
        pmanager->AddDiscreteProcess(fBoundaryProcess);
        pmanager->AddDiscreteProcess(wls);
        pmanager->AddDiscreteProcess(fastSimProcess_garfield);
      }
    }
  
  }

} // end namespace nexus
