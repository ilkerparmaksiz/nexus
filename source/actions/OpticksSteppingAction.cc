//
// Created by ilker parmaksiz on 2/29/24.
// This Class handles S1 photons for
//

#include "OpticksSteppingAction.hh"
#include "FactoryBase.h"

#include <G4Step.hh>
#include <G4SteppingManager.hh>
#include <G4ProcessManager.hh>
#include <G4OpticalPhoton.hh>
#include <G4OpBoundaryProcess.hh>
#include "config.h"
#include <G4VPhysicalVolume.hh>
#ifdef With_GarField
#include "NESTS1Photon.hh"
#endif
#include "G4EventManager.hh"
#include "G4Scintillation.hh"
#include "PersistencyManager.h"

#ifdef With_Opticks
    #include "G4CXOpticks.hh"
    #include "SEvt.hh"
    #include "U4.hh"
#endif

using namespace nexus;
REGISTER_CLASS(OpticksSteppingAction, G4UserSteppingAction)

OpticksSteppingAction::OpticksSteppingAction(): G4UserSteppingAction()
{
}

OpticksSteppingAction::~OpticksSteppingAction()
{

}
void OpticksSteppingAction::UserSteppingAction(const G4Step* step) {
    G4ParticleDefinition *pdef = step->GetTrack()->GetDefinition();
    G4Track *track = step->GetTrack();
    PersistencyManager *pManger= dynamic_cast<PersistencyManager *>(PersistencyManager::GetPersistencyManager());

    //Check whether the track is an optical photon
    if (pdef != G4OpticalPhoton::Definition()) {
        G4SteppingManager *sMg = G4EventManager::GetEventManager()->GetTrackingManager()->GetSteppingManager();
        G4StepStatus stepStatus = sMg->GetfStepStatus();
        if (stepStatus != fAtRestDoItProc) {
            G4ProcessVector *PostStepProc = sMg->GetfPostStepDoItVector();
            size_t MaxSteps = sMg->GetMAXofPostStepLoops();
            for (int stp = 0; stp < MaxSteps; stp++) {
                if ((*PostStepProc)[stp]->GetProcessName() == "Scintillation") {
                    G4Scintillation *ScintProc = (G4Scintillation *) (*PostStepProc)[stp];
                    G4int num_photons = ScintProc->GetNumPhotons();
                    //std::cout << "Scintilation "<< num_photons <<std::endl;

                    if (num_photons > 0) {
                        G4MaterialPropertiesTable *MPT = track->GetMaterial()->GetMaterialPropertiesTable();
                        G4double t1, t2 = 0;
                        G4int singlets, triplets = 0;
                        t1 = MPT->GetConstProperty(kSCINTILLATIONTIMECONSTANT1);
                        t2 = MPT->GetConstProperty(kSCINTILLATIONTIMECONSTANT2);
                        singlets = floor(MPT->GetConstProperty(kSCINTILLATIONYIELD1) * num_photons);
                        triplets = ceil(MPT->GetConstProperty(kSCINTILLATIONYIELD2) * num_photons);
                        pManger->AddPhotons((singlets+triplets)); // S1

                        //std::cout << "Scintilation "<< num_photons <<" Amount of Singlets " <<singlets <<" Triplets " << triplets <<std::endl;
                        #ifdef With_Opticks
                        if (singlets > 0)
                            U4::CollectGenstep_DsG4Scintillation_r4695(track, step, singlets, 0, t1);
                        if (triplets > 0)
                            U4::CollectGenstep_DsG4Scintillation_r4695(track, step, triplets, 1, t2);

                        #endif

                    }

                }
            }
        }
    }else{
#ifdef With_Opticks
        if(step->GetTrack()->GetDefinition()==G4OpticalPhoton::Definition()) step->GetTrack()->SetTrackStatus(fStopAndKill);
#endif
#ifndef With_Opticks
    static G4OpBoundaryProcess* boundary = 0;

    if (!boundary) { // the pointer is not defined yet
        // Get the list of processes defined for the optical photon
        // and loop through it to find the optical boundary process.
        G4ProcessVector* pv = pdef->GetProcessManager()->GetProcessList();
        for (size_t i=0; i<pv->size(); i++) {

            if ((*pv)[i]->GetProcessName() == "OpBoundary") {
                boundary = (G4OpBoundaryProcess*) (*pv)[i];
                break;
            }
        }
    }

    if (step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {
        if (boundary and boundary->GetStatus() == Detection ){
            G4String detector_name = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName();


            //G4cout << "##### Sensitive Volume: " << detector_name << G4endl;
            pManger->AddOpticalHit(detector_name,step->GetPostStepPoint()->GetPosition(),step->GetPostStepPoint()->GetGlobalTime());
        }
    }


#endif
        return;
    }
}