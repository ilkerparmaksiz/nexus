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
    #include "omp.h"
#endif

using namespace nexus;
REGISTER_CLASS(OpticksSteppingAction, G4UserSteppingAction)

OpticksSteppingAction::OpticksSteppingAction(): G4UserSteppingAction()
{}

OpticksSteppingAction::~OpticksSteppingAction()
{}
void OpticksSteppingAction::UserSteppingAction(const G4Step* step) {
    G4ParticleDefinition *pdef = step->GetTrack()->GetDefinition();
    G4Track *track = step->GetTrack();
    PersistencyManager *pManger= dynamic_cast<PersistencyManager *>(PersistencyManager::GetPersistencyManager());

#ifdef With_GarField
     if(pdef==NESTS1Photon::Definition()) {
         step->GetTrack()->SetTrackStatus(fStopAndKill);
         return;
     }
#endif
    //Check whether the track is an optical photon
    if (pdef != G4OpticalPhoton::Definition()) {
        G4SteppingManager *sMg = G4EventManager::GetEventManager()->GetTrackingManager()->GetSteppingManager();
        G4StepStatus stepStatus = sMg->GetfStepStatus();
        G4double t1, t2 = 0;
        G4int singlets=0, triplets=0,TotalPhotns = 0;
        G4MaterialPropertiesTable *MPT = track->GetMaterial()->GetMaterialPropertiesTable();

        if (stepStatus != fAtRestDoItProc) {
            G4ProcessVector *PostStepProc = sMg->GetfPostStepDoItVector();
            size_t MaxSteps = sMg->GetMAXofPostStepLoops();
            //#pragma omp parallel for num_threads(15)
            for (int stp = 0; stp < MaxSteps; stp++) {
                if ((*PostStepProc)[stp]->GetProcessName() == "Scintillation") {
                    G4Scintillation *ScintProc = (G4Scintillation *) (*PostStepProc)[stp];
                    G4int num_photons = ScintProc->GetNumPhotons();

                    //std::cout << "Scintilation "<< num_photons <<std::endl;

                    if (num_photons > 0) {
                        TotalPhotns+=num_photons;
                        //std::cout << "Scintilation PreStep "<< step->GetPreStepPoint()->GetPosition() << " PostStep " << step->GetPostStepPoint()->GetPosition() <<" TotalNumber " << num_photons <<std::endl;
                    }

                }
            }

        #ifdef With_Opticks
            t1 = MPT->GetConstProperty(kSCINTILLATIONTIMECONSTANT1);
            t2 = MPT->GetConstProperty(kSCINTILLATIONTIMECONSTANT2);
            singlets = floor(MPT->GetConstProperty(kSCINTILLATIONYIELD1) * TotalPhotns);
            triplets = ceil(MPT->GetConstProperty(kSCINTILLATIONYIELD2) * TotalPhotns);
            pManger->AddPhotons((singlets+triplets)); // S1
            #ifdef With_G4OpticksTest
            pManger->fOpticksPhotonCounter+=(singlets+triplets);
            #endif
            //std::unique_ptr<G4Step> newStep= std::make_unique<G4Step>(G4Step(*step)) ;
           // G4Step *newStep= new G4Step(*step) ;
            //std::cout <<"New Step" <<std::endl;
            //G4StepPoint newStepPoint=*newStep->GetPreStepPoint();
            //newStepPoint.SetPosition(newStep->GetPreStepPoint()->GetPosition()+(newStep->GetPostStepPoint()->GetPosition()-newStep->GetPreStepPoint()->GetPosition())*G4UniformRand());
            //newStepPoint.SetGlobalTime(newStep->GetPreStepPoint()->GetGlobalTime()+(newStep->GetPostStepPoint()->GetGlobalTime()-newStep->GetPreStepPoint()->GetGlobalTime())*G4UniformRand());
            //newStep->SetPreStepPoint(&newStepPoint);
            //std::cout << "Triplet " <<triplets << " Singlets " << singlets <<std::endl;
            if (singlets > 0)
                U4::CollectGenstep_DsG4Scintillation_r4695(track, step, singlets, 0, t1);
            if (triplets > 0)
                U4::CollectGenstep_DsG4Scintillation_r4695(track, step, triplets, 1, t2);
        #endif
        }
        //std::cout << "PreStep "<< step->GetPreStepPoint()->GetPosition() << " PostStep " << step->GetPostStepPoint()->GetPosition() << std::endl;

    }else{
    #ifdef With_G4OpticksTest
            if( G4OpticalPhoton::Definition() and track->GetTrackStatus()==fStopAndKill) pManger->fG4PhotonCounter+=1;
    #endif

#if defined(With_Opticks)  and not defined(With_G4OpticksTest)

        if(step->GetTrack()->GetDefinition()==G4OpticalPhoton::Definition()) step->GetTrack()->SetTrackStatus(fStopAndKill);
#endif
#if not defined(With_Opticks) || defined(With_G4OpticksTest)
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
            auto PostStep=step->GetPostStepPoint();
            G4double wavelength=1239.8/PostStep->GetTotalEnergy()*CLHEP::eV; //nm
            pManger->AddOpticalHit(detector_name, PostStep->GetPosition(),PostStep->GetGlobalTime(),PostStep->GetMomentum(),PostStep->GetPolarization(),wavelength);
        }
    }


#endif
        return;
    }
}