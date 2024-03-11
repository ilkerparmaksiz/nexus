// ----------------------------------------------------------------------------
// nexus | AnalysisSteppingAction.cc
//
// This class allows the user to print the total number of photons detected by
// all kinds of photosensors at the end of the run.
// It also shows examples of information that can be accessed at the stepping
// level, so it is useful for debugging.
//
// The  NEXT Collaboration
// ----------------------------------------------------------------------------

#include "AnalysisSteppingAction.h"
#include "FactoryBase.h"

#include <G4Step.hh>
#include <G4SteppingManager.hh>
#include <G4ProcessManager.hh>
#include <G4OpticalPhoton.hh>
#include <G4OpBoundaryProcess.hh>
#include "config.h"
#include <G4VPhysicalVolume.hh>
#ifdef With_Opticks
#include "G4CXOpticks.hh"
#include "SEvt.hh"
#include "U4.hh"
#include "G4Scintillation.hh"
#include "G4EventManager.hh"
#endif
using namespace nexus;

REGISTER_CLASS(AnalysisSteppingAction, G4UserSteppingAction)

AnalysisSteppingAction::AnalysisSteppingAction(): G4UserSteppingAction()
{
}



AnalysisSteppingAction::~AnalysisSteppingAction()
{
  G4double total_counts = 0;
  detectorCounts::iterator it = my_counts_.begin();
  while (it != my_counts_.end()) {
    G4cout << "Detector " << it->first << ": " << it->second << " counts" << G4endl;
    total_counts += it->second;
    it ++;
  }
  G4cout << "TOTAL COUNTS: " << total_counts << G4endl;
}



void AnalysisSteppingAction::UserSteppingAction(const G4Step* step)
{
  G4ParticleDefinition* pdef = step->GetTrack()->GetDefinition();
    G4Track* track = step->GetTrack();
#ifdef With_Opticks
    if (pdef==G4OpticalPhoton::Definition())
        track->SetTrackStatus(fStopAndKill);
#endif
  //Check whether the track is an optical photon
  if (pdef != G4OpticalPhoton::Definition()) {
    #ifdef With_Opticks
    G4SteppingManager * sMg=G4EventManager::GetEventManager()->GetTrackingManager()->GetSteppingManager();
    G4StepStatus stepStatus=sMg->GetfStepStatus();
    if(stepStatus!=fAtRestDoItProc){
        G4ProcessVector * PostStepProc=sMg->GetfPostStepDoItVector();
        size_t MaxSteps=sMg->GetMAXofPostStepLoops();
        for (int stp=0;stp<MaxSteps;stp++){
            if((*PostStepProc)[stp]->GetProcessName()=="Scintillation"){
                G4Scintillation *ScintProc= (G4Scintillation*) (*PostStepProc)[stp];
                G4int num_photons=ScintProc->GetNumPhotons();
                //std::cout << "Scintilation "<< num_photons <<std::endl;

                if(num_photons>0){

                    G4MaterialPropertiesTable * MPT = track->GetMaterial()->GetMaterialPropertiesTable();
                    G4double t1,t2=0;
                    G4int singlets,triplets=0;
                    t1=MPT->GetConstProperty(kSCINTILLATIONTIMECONSTANT1);
                    t2=MPT->GetConstProperty(kSCINTILLATIONTIMECONSTANT2);
                    singlets= floor(MPT->GetConstProperty(kSCINTILLATIONYIELD1)*num_photons);
                    triplets= ceil(MPT->GetConstProperty(kSCINTILLATIONYIELD2)*num_photons);
                    //std::cout << "Scintilation "<< num_photons <<" Amount of Singlets " <<singlets <<" Triplets " << triplets <<std::endl;
                    if(singlets>0)
                        U4::CollectGenstep_DsG4Scintillation_r4695(track,step,singlets,0,t1);
                    if(triplets>0)
                        U4::CollectGenstep_DsG4Scintillation_r4695(track,step,triplets,1,t2);
                }

            }
        }
    }
    #endif
    return ;
  }



  /*
  // example of information one can access about optical photons

  G4Track* track = step->GetTrack();
  G4int pid = track->GetParentID();
  G4int tid = track->GetTrackID();
  G4StepPoint* point1 = step->GetPreStepPoint();
  G4StepPoint* point2 = step->GetPostStepPoint();
  G4TouchableHandle touch1 = point1->GetTouchableHandle();
  G4TouchableHandle touch2 = point2->GetTouchableHandle();
  G4String vol1name = touch1->GetVolume()->GetName();
  G4String vol2name = touch2->GetVolume()->GetName();

  G4String proc_name = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
  G4int copy_no = step->GetPostStepPoint()->GetTouchable()->GetReplicaNumber(1);
  */

  // Retrieve the pointer to the optical boundary process.
  // We do this only once per run defining our local pointer as static.
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
    if (boundary->GetStatus() == Detection ){
      G4String detector_name = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName();
      //G4cout << "##### Sensitive Volume: " << detector_name << G4endl;

      detectorCounts::iterator it = my_counts_.find(detector_name);
      if (it != my_counts_.end()) my_counts_[it->first] += 1;
      else my_counts_[detector_name] = 1;
    }
  }

}
