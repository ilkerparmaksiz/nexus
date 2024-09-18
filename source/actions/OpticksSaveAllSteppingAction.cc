// ----------------------------------------------------------------------------
// nexus | OpticksSaveAllSteppingAction.cc
//
// This class adds a new group and table to the output file, "/DEBUG/steps".
// This table contains information (position and volume of both the
// pre- and post-step points, average time, process name and other identifiers)
// of some steps of the simulation. By default all steps are stored. However,
// a subset of them can be selected by cherry-picking the volumes and particles
// involved in the step. This can be achieved with the commands
// /Actions/SaveAllSteppingAction/select_particle
// and
// /Actions/OpticksSaveAllSteppingAction/select_volume
// without the need for re-compilation.
// It must be noted that the files produced with this action become large
// very quickly. Therefore, strict filtering and small number of events are
// encouraged.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "OpticksSaveAllSteppingAction.h"
#include "PersistencyManager.h"
#include "FactoryBase.h"

#include <G4Step.hh>
#include <G4VPersistencyManager.hh>
#include <G4ProcessManager.hh>
#include <G4ParticleTable.hh>

#include "G4EventManager.hh"
#include "G4Scintillation.hh"
#include <G4SteppingManager.hh>
#include <G4OpticalPhoton.hh>
#include <G4OpBoundaryProcess.hh>
#include "config.h"
#ifdef With_Opticks
#include "G4CXOpticks.hh"
    #include "SEvt.hh"
    #include "U4.hh"
#endif
using namespace nexus;

REGISTER_CLASS(OpticksSaveAllSteppingAction, G4UserSteppingAction)

OpticksSaveAllSteppingAction::OpticksSaveAllSteppingAction():
G4UserSteppingAction(),
msg_(0),
selected_volumes_(),
selected_particles_(),
initial_volumes_(),
final_volumes_(),
proc_names_(),
initial_poss_(),
final_poss_(),
times_(),
kill_after_selection_(false)
{
  msg_ = new G4GenericMessenger(this, "/Actions/OpticksSaveAllSteppingAction/");

  msg_->DeclareMethod("select_particle",
                      &OpticksSaveAllSteppingAction::AddSelectedParticle,
                      "add a new particle to select");

  msg_->DeclareMethod("select_volume",
                      &OpticksSaveAllSteppingAction::AddSelectedVolume,
                      "add a new volume to select");

  msg_->DeclareProperty("kill_after_selection",
                        kill_after_selection_,
                        "Whether to kill a particle after a step has been selected");

  PersistencyManager* pm = dynamic_cast<PersistencyManager*>
        (G4VPersistencyManager::GetPersistencyManager());

  pm->StoreSteps(true);

}



OpticksSaveAllSteppingAction::~OpticksSaveAllSteppingAction()
{
}



void OpticksSaveAllSteppingAction::UserSteppingAction(const G4Step* step)
{
  G4ParticleDefinition* pdef          = step->GetTrack()->GetDefinition();
  G4int                 track_id      = step->GetTrack()->GetTrackID();
  G4String              particle_name = pdef->GetParticleName();
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
                G4double wavelength=1239.8/step->GetPostStepPoint()->GetTotalEnergy()*CLHEP::eV; //nm

                //G4cout << "##### Sensitive Volume: " << detector_name << G4endl;
                pManger->AddOpticalHit(detector_name, PostStep->GetPosition(),PostStep->GetGlobalTime(),PostStep->GetMomentum(),PostStep->GetPolarization(),wavelength);
            }
        }


#endif
        //return;
    }

  if (!KeepParticle(pdef)) return;

  G4StepPoint* pre  = step->GetPreStepPoint();
  G4StepPoint* post = step->GetPostStepPoint();

  G4ThreeVector initial_pos = pre ->GetPosition();
  G4ThreeVector   final_pos = post->GetPosition();
  G4double        step_time = (pre->GetGlobalTime()  +
                              post->GetGlobalTime()) / 2.;

  if (! post->GetTouchableHandle()->GetVolume()) return; // Particle exits the world

  G4String initial_volume = pre ->GetTouchableHandle()->GetVolume()->GetName();
  G4String   final_volume = post->GetTouchableHandle()->GetVolume()->GetName();
  G4String      proc_name = post->GetProcessDefinedStep()->GetProcessName();

  if (!KeepVolume(initial_volume, final_volume))
    return;

  std::pair<G4int, G4String> key = std::make_pair(track_id, particle_name);

  initial_volumes_[key].push_back(initial_volume);
    final_volumes_[key].push_back(  final_volume);
       proc_names_[key].push_back(     proc_name);

  initial_poss_   [key].push_back(initial_pos);
    final_poss_   [key].push_back(  final_pos);
         times_   [key].push_back(  step_time);

  if (kill_after_selection_)
    step->GetTrack()->SetTrackStatus(fStopAndKill);
}


void OpticksSaveAllSteppingAction::AddSelectedParticle(G4String particle_name)
{
  G4ParticleDefinition* pdef = G4ParticleTable::GetParticleTable()->FindParticle(particle_name);
  if (!pdef) {
    G4String msg = "No particle description was found for particle name " + particle_name;
    G4Exception("[OpticksSaveAllSteppingAction]", "AddSelectedParticle()", FatalException, msg);
  }
  selected_particles_.push_back(pdef);
}


void OpticksSaveAllSteppingAction::AddSelectedVolume(G4String volume_name)
{
  selected_volumes_.push_back(volume_name);
}


G4bool OpticksSaveAllSteppingAction::KeepParticle(G4ParticleDefinition* pdef)
{
  if (!selected_particles_.size()) return true;

  auto it = std::find(selected_particles_.begin(), selected_particles_.end(), pdef);
  return it != selected_particles_.end();
}


G4bool OpticksSaveAllSteppingAction::KeepVolume(G4String& initial_volume, G4String& final_volume)
{
  if (!selected_volumes_.size()) return true;

  for (auto volume=selected_volumes_.begin(); volume != selected_volumes_.end(); volume++)
  {
    if (G4StrUtil::contains(initial_volume, *volume)) return true;
    if (G4StrUtil::contains(  final_volume, *volume)) return true;
  }

  return false;
}



void OpticksSaveAllSteppingAction::Reset()
{
  initial_volumes_.clear();
    final_volumes_.clear();
       proc_names_.clear();

  initial_poss_   .clear();
    final_poss_   .clear();
         times_   .clear();
}
