//
// Created by ilker on 8/6/23.
//

#include "GarfieldPhysics.h"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4PhysicsListHelper.hh"
#include "G4SystemOfUnits.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4RayleighScattering.hh"
#include "G4KleinNishinaModel.hh"

#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuMultipleScattering.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"

#include "G4ionIonisation.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4NuclearStopping.hh"

#include "G4EmParameters.hh"
#include "G4MscStepLimitType.hh"

#include "G4LossTableManager.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4SystemOfUnits.hh"

#include "G4EmModelActivator.hh"

#include "G4FastSimulationManagerProcess.hh"
#include "G4FastSimulationPhysics.hh"

#include "gasNESTdet.hh"
#include "G4/NESTProc.hh"
#include "G4OpRayleigh.hh"
#include "OpAbsorption.hh"
#include "OpWLS.hh"
#include "OpBoundaryProcess.hh"
#include "S2Photon.hh"
#include <G4PhysicsConstructorFactory.hh>
#include <G4GenericMessenger.hh>

#ifdef theParticleIterator
#undef theParticleIterator
#endif

namespace nexus {
    /// Macro that allows the use of this physics constructor
    /// with the generic physics list
    G4_DECLARE_PHYSCONSTR_FACTORY(GarfieldPhysics);

    GarfieldPhysics::GarfieldPhysics()
            : IonGasModels_(false),G4VPhysicsConstructor("GarfieldPhysics"){
        msg_ = new G4GenericMessenger(this, "/PhysicsList/Garfield/",
                                      "Control commands of the nexus physics list.");
        msg_->DeclareProperty("IonGasModels", IonGasModels_,
                              "Switch on/off the ionization drift.");
        fastSimulationPhysics = new G4FastSimulationPhysics("fastSimPhys");
        if(IonGasModels_) AddIonGasModels();
    }
    void GarfieldPhysics::ConstructParticle() {
        NEST::NESTThermalElectron::Definition();
        G4OpticalPhoton::Definition();
        S2Photon::Definition();
    }

    void GarfieldPhysics::ConstructProcess() {

        ///
        G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

        // Add standard EM Processes
        //
        G4ParticleTable::G4PTblDicIterator* theParticleIterator = theParticleTable->GetIterator();
        theParticleIterator->reset();
        while ((*theParticleIterator)() ) {
            G4ParticleDefinition* particle = theParticleIterator->value();
            G4String particleName = particle->GetParticleName();

            if (particleName == "gamma") {
                ////ph->RegisterProcess(new G4RayleighScattering, particle);
                ph->RegisterProcess(new G4PhotoElectricEffect, particle);
                G4ComptonScattering* cs = new G4ComptonScattering;
                cs->SetEmModel(new G4KleinNishinaModel());
                ph->RegisterProcess(cs, particle);
                ph->RegisterProcess(new G4GammaConversion, particle);

            } else if (particleName == "e-") {
                ph->RegisterProcess(new G4eMultipleScattering(), particle);
                G4eIonisation* eIoni = new G4eIonisation();
                eIoni->SetStepFunction(0.1, 100 * um);
                eIoni->SetDEDXBinning(12 * 10);
                ph->RegisterProcess(eIoni, particle);
                ph->RegisterProcess(new G4eBremsstrahlung(), particle);

            } else if (particleName == "e+") {
                ph->RegisterProcess(new G4eMultipleScattering(), particle);
                G4eIonisation* eIoni = new G4eIonisation();
                eIoni->SetStepFunction(0.1, 100 * um);
                eIoni->SetDEDXBinning(12 * 10);
                ph->RegisterProcess(eIoni, particle);
                ph->RegisterProcess(new G4eBremsstrahlung(), particle);
                ph->RegisterProcess(new G4eplusAnnihilation(), particle);

            } else if (particleName == "mu+" || particleName == "mu-") {
                ph->RegisterProcess(new G4MuMultipleScattering(), particle);
                G4MuIonisation* muIoni = new G4MuIonisation();
                muIoni->SetStepFunction(0.1, 50 * um);
                ph->RegisterProcess(muIoni, particle);
                ph->RegisterProcess(new G4MuBremsstrahlung(), particle);
                ph->RegisterProcess(new G4MuPairProduction(), particle);

            } else if (particleName == "proton" || particleName == "pi-" ||
                       particleName == "pi+") {
                ph->RegisterProcess(new G4hMultipleScattering(), particle);
                G4hIonisation* hIoni = new G4hIonisation();
                hIoni->SetStepFunction(0.1, 20 * um);
                ph->RegisterProcess(hIoni, particle);
                ph->RegisterProcess(new G4hBremsstrahlung(), particle);
                ph->RegisterProcess(new G4hPairProduction(), particle);

            } else if (particleName == "alpha" || particleName == "He3") {
                ph->RegisterProcess(new G4hMultipleScattering(), particle);
                G4ionIonisation* ionIoni = new G4ionIonisation();
                ionIoni->SetStepFunction(0.1, 1 * um);
                ionIoni->SetDEDXBinning(12 * 10);
                ph->RegisterProcess(ionIoni, particle);
                ph->RegisterProcess(new G4NuclearStopping(), particle);

            } else if (particleName == "GenericIon") {
                ph->RegisterProcess(new G4hMultipleScattering(), particle);
                G4ionIonisation* ionIoni = new G4ionIonisation();
                ionIoni->SetEmModel(new G4IonParametrisedLossModel());
                ionIoni->SetStepFunction(0.1, 1 * um);
                ph->RegisterProcess(ionIoni, particle);
                ph->RegisterProcess(new G4NuclearStopping(), particle);

            } else if ((!particle->IsShortLived()) &&
                       (particle->GetPDGCharge() != 0.0) &&
                       (particle->GetParticleName() != "chargedgeantino")) {
                // all others charged particles except geantino
                ph->RegisterProcess(new G4hMultipleScattering(), particle);
                ph->RegisterProcess(new G4hIonisation(), particle);
            }else if (particleName=="S2Photon"){
                ph->RegisterProcess(new OpAbsorption(),particle);
                ph->RegisterProcess(new OpWLS(),particle);
                ph->RegisterProcess(new OpBoundaryProcess(),particle);
            }
            //AddParametrisation();
        }

        // Em options
        //
        // Main options and setting parameters are shown here.
        // Several of them have default values.
        //
        G4EmParameters *emOptions= G4EmParameters::Instance();

        // physics tables
        //
        emOptions->SetMinEnergy(10 * eV);      // default 100 eV
        emOptions->SetMaxEnergy(10 * TeV);     // default 100 TeV
        // emOptions.SetDEDXBinning(12 * 10);    // default=12*7
        // emOptions.SetLambdaBinning(12 * 10);  // default=12*7

        // multiple coulomb scattering
        //
        emOptions->SetMscStepLimitType(fUseSafety); // default  // default

        // Deexcitation
        //
        G4VAtomDeexcitation* de = new G4UAtomicDeexcitation();
        de->SetFluo(true);
        de->SetAuger(false);
        de->SetPIXE(true);
        G4LossTableManager::Instance()->SetAtomDeexcitation(de);

        G4EmModelActivator mact(GetPhysicsName());


        auto particleIteratorP=GetParticleIterator();
        particleIteratorP->reset();
        std::cout  << "PhysicsListEMStandard::ConstructProcess() pit is "  << particleIteratorP << std::endl;

        gasNESTdet* gndet = new gasNESTdet();
        // std::shared_ptr<gasNESTdet> gndet(new gasNESTdet());
        NEST::NESTcalc* calcNEST = new NEST::NESTcalc(gndet);

        NEST::NESTProc* theNEST2ScintillationProcess = new NEST::NESTProc("S1",fElectromagnetic, calcNEST, gndet); //gndet);
        theNEST2ScintillationProcess->SetDetailedSecondaries(true);
        theNEST2ScintillationProcess->SetStackElectrons(true);

        while( (*particleIteratorP)() ){
            G4ParticleDefinition* particle = particleIteratorP->value();
            G4ProcessManager* pmanager = particle->GetProcessManager();
            G4String particleName = particle->GetParticleName();

            // td::cout << "PhysicsListEMStandard::ConstructProcess(): pname, pmanager are " << particleName << ", " << pmanager << std::endl;
            // if ( !( particleName.find("e-")!=std::string::npos || particleName.find("alpha")!=std::string::npos  || particleName.find("opticalphoton")!=std::string::npos ) )
            //   continue;
            if (pmanager) {
                if (theNEST2ScintillationProcess->IsApplicable(*particle) && pmanager) {
                    std::cout << "PhysicsList::InitialisePhysics(): particleName, pmanager  " << particleName << ", " << pmanager << "." << std::endl;
                    std::cout << "ordDefault, ordInActive " << ordDefault << ", " << ordInActive  << std::endl;
                    pmanager->AddProcess(theNEST2ScintillationProcess, ordDefault + 1, ordInActive, ordDefault + 1);
                }

                G4OpBoundaryProcess* fBoundaryProcess = new OpBoundaryProcess();
                G4OpAbsorption* fAbsorptionProcess = new OpAbsorption();
                G4OpWLS* fTheWLSProcess = new OpWLS();

                if (particleName == "opticalphoton" && pmanager) {
                    G4cout << " AddDiscreteProcess to OpticalPhoton " << G4endl;
                    pmanager->AddDiscreteProcess(fAbsorptionProcess);
                    //      pmanager->AddDiscreteProcess(fRayleighScatteringProcess);
                    pmanager->AddDiscreteProcess(fTheWLSProcess);
                    pmanager->AddDiscreteProcess(fBoundaryProcess);
                }


                //      std::cout << "PhysicsList::InitialisePhysicsList(): e-/optphot in particleIteratorP" << std::endl;
                std::cout << "PhysicsListEMStandard()::ConstructProcess() process list (of length) " << pmanager->GetProcessList()->size() << " for " << particleName << " is: " << std::endl;
                //    for (const auto  proc : pmanager->GetProcessList())
                for (short int ii = 0; ii<(short int)(pmanager->GetProcessList()->size());++ii)
                {
                    std::cout << (*pmanager->GetProcessList())[ii]->GetProcessName() << std::endl;
                }
            }
        }



    }



    void GarfieldPhysics::AddIonGasModels() {
        G4EmConfigurator* em_config =
                G4LossTableManager::Instance()->EmConfigurator();
        G4ParticleTable::G4PTblDicIterator* theParticleIterator = theParticleTable->GetIterator();
        theParticleIterator->reset();
        while ((*theParticleIterator)()) {
            G4ParticleDefinition* particle = theParticleIterator->value();
            G4String partname = particle->GetParticleName();
            if (partname == "alpha" || partname == "He3" || partname == "GenericIon") {
                G4BraggIonGasModel* mod1 = new G4BraggIonGasModel();
                G4BetheBlochIonGasModel* mod2 = new G4BetheBlochIonGasModel();
                G4double eth = 2. * MeV * particle->GetPDGMass() / proton_mass_c2;
                em_config->SetExtraEmModel(partname, "ionIoni", mod1, "", 0.0, eth,
                                           new G4IonFluctuations());
                em_config->SetExtraEmModel(partname, "ionIoni", mod2, "", eth, 100 * TeV,
                                           new G4UniversalFluctuation());
            }
        }
    }

    void GarfieldPhysics::AddParametrisation() {
        theParticleTable->GetIterator()->reset();
        while ((*theParticleTable->GetIterator())()) {
            G4String particleName = theParticleTable->GetIterator()->value()->GetParticleName();
            fastSimulationPhysics->ActivateFastSimulation(particleName);
        }
    }
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



} // nexus