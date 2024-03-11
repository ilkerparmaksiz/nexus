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

#ifdef With_Opticks
#include "NESTS1Photon.hh"
    #include "G4CXOpticks.hh"
    #include "SEvt.hh"
    #include "U4.hh"
    #include "G4Scintillation.hh"
    #include "G4EventManager.hh"
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
#ifdef With_Opticks
    if ( pdef==NESTS1Photon::Definition())
        track->SetTrackStatus(fStopAndKill);
#endif
    //Check whether the track is an optical photon
    if (pdef != G4OpticalPhoton::Definition()) {
#ifdef With_Opticks
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
                        //std::cout << "Scintilation "<< num_photons <<" Amount of Singlets " <<singlets <<" Triplets " << triplets <<std::endl;
                        if (singlets > 0)
                            U4::CollectGenstep_DsG4Scintillation_r4695(track, step, singlets, 0, t1);
                        if (triplets > 0)
                            U4::CollectGenstep_DsG4Scintillation_r4695(track, step, triplets, 1, t2);
                    }

                }
            }
        }
#endif
    }

}