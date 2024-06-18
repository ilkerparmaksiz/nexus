// ----------------------------------------------------------------------------
// nexus | DefaultEventAction.cc
//
// This is the default event action of the NEXT simulations. Only events with
// deposited energy larger than 0 are saved in the nexus output file.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "DefaultEventAction.h"
#include "Trajectory.h"
#include "PersistencyManager.h"
#include "IonizationHit.h"
#include "FactoryBase.h"

#include <G4Event.hh>
#include <G4VVisManager.hh>
#include <G4Trajectory.hh>
#include <G4GenericMessenger.hh>

#include <globals.hh>
// Krishan Garfield specific code
#include "config.h"
#ifdef With_GarField
#include "DegradModel.h"
#endif
#ifdef With_Opticks
    #include "G4GlobalFastSimulationManager.hh"
    #  include "SEvt.hh"
    #  include "G4CXOpticks.hh"
namespace {G4Mutex opticks_mt =G4MUTEX_INITIALIZER;}
#endif

namespace nexus {

REGISTER_CLASS(DefaultEventAction, G4UserEventAction)

  DefaultEventAction::DefaultEventAction():
    G4UserEventAction(), nevt_(0), nupdate_(10), energy_min_(0.), energy_max_(DBL_MAX)
  {
    msg_ = new G4GenericMessenger(this, "/Actions/DefaultEventAction/");

    G4GenericMessenger::Command& thresh_cmd =
       msg_->DeclareProperty("min_energy", energy_min_,
                             "Minimum deposited energy to save the event to file.");
    thresh_cmd.SetParameterName("min_energy", true);
    thresh_cmd.SetUnitCategory("Energy");
    thresh_cmd.SetRange("min_energy>0.");

    G4GenericMessenger::Command& max_energy_cmd =
      msg_->DeclareProperty("max_energy", energy_max_,
                            "Maximum deposited energy to save the event to file.");
    max_energy_cmd.SetParameterName("max_energy", true);
    max_energy_cmd.SetUnitCategory("Energy");
    max_energy_cmd.SetRange("max_energy>0.");

    PersistencyManager* pm = dynamic_cast<PersistencyManager*>
      (G4VPersistencyManager::GetPersistencyManager());

    pm->SaveNumbOfInteractingEvents(true);
  }



  DefaultEventAction::~DefaultEventAction()
  {
  }



  void DefaultEventAction::BeginOfEventAction(const G4Event* event)
  {
    startTime=std::chrono::high_resolution_clock::now();
    // Print out event number info
    if ((nevt_ % nupdate_) == 0) {
      G4cout << " >> Event no. " << nevt_  << G4endl;
      if (nevt_  == (10 * nupdate_)) nupdate_ *= 10;
    }

    // Krishan: Garfield specific code
#ifdef With_GarField
    DegradModel* dm = (DegradModel*)(G4GlobalFastSimulationManager::GetInstance()->GetFastSimulationModel("DegradModel"));
    if(dm)
        dm->Reset();
#endif
    G4PrimaryVertex* pVtx;
    pVtx = event->GetPrimaryVertex();
    if (pVtx)
      {
        G4double PKE = pVtx->GetPrimary(0)->GetKineticEnergy();
        #ifdef With_GarField
            dm->SetPrimaryKE(PKE/eV);
        #endif
      }



  }



  void DefaultEventAction::EndOfEventAction(const G4Event* event)
  {
    nevt_++;


    // Determine whether total energy deposit in ionization sensitive
    // detectors is above threshold
    if (energy_min_ >= 0.) {

        #ifdef With_Opticks

            G4AutoLock lock(&opticks_mt);
            G4CXOpticks * g4cx=G4CXOpticks::Get();

            G4int eventID=event->GetEventID();
            G4int ngenstep=SEvt::GetNumGenstepFromGenstep(0);
            G4int nphotons=SEvt::GetNumPhotonCollected(0);



            // Simulate the photons
            if(nphotons>0 and ngenstep>0){
                //std::cout<<g4cx->desc()<<std::endl;
                //std::cout<<"--- G4Optickx ---" << g4cx->descSimulate() <<std::endl;
                g4cx->simulate(eventID,0); // For Simulation
                cudaDeviceSynchronize();

               // std::cout<<"Event " <<eventID <<" Simulating with Opticks nphotons "<< nphotons << " nsteps " << ngenstep << " Hits " <<SEvt::GetNumHit(0) << std::endl;

            }
        //G4cout<<" Opticks End of Event Action" <<G4endl;



        #endif

      // Get the trajectories stored for this event and loop through them
      // to calculate the total energy deposit

      G4double edep = 0.;

      G4TrajectoryContainer* tc = event->GetTrajectoryContainer();
      if (tc) {
        // in interactive mode, a G4TrajectoryContainer would exist
        // but the trajectories will not cast to Trajectory
        Trajectory* trj = dynamic_cast<Trajectory*>((*tc)[0]);
        if (trj == nullptr){
          G4Exception("[DefaultEventAction]", "EndOfEventAction()", FatalException,
                      "DefaultTrackingAction is required when using DefaultEventAction");
        }
        for (unsigned int i=0; i<tc->size(); ++i) {
          Trajectory* tr = dynamic_cast<Trajectory*>((*tc)[i]);
          edep += tr->GetEnergyDeposit();
        }
      }
      else {
        G4Exception("[DefaultEventAction]", "EndOfEventAction()", FatalException,"DefaultTrackingAction is required when using DefaultEventAction");
      }

      PersistencyManager* pm = dynamic_cast<PersistencyManager*>
        (G4VPersistencyManager::GetPersistencyManager());

      if (!event->IsAborted() && edep>0) {
        pm->InteractingEvent(true);
      } else {
        pm->InteractingEvent(false);
      }
      if (!event->IsAborted() && edep > energy_min_ && edep < energy_max_) {
        pm->StoreCurrentEvent(true);
      } else {
        pm->StoreCurrentEvent(false);
      }
        endTime=std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = endTime - startTime;
        pm->SetCompletionTime(duration.count());
        pm->SaveTimeInfo();
    }



  }


} // end namespace nexus
