#include "GarfieldEventAction.h"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"
#include "GarfieldRunAction.h"

#include "G4SDManager.hh"
#include "G4Threading.hh"
#include "GarfieldSteppingAction.h"
#include "G4VPhysicalVolume.hh"
#include "G4GlobalFastSimulationManager.hh"
#include "GarfieldVUVPhotonModel.h"
#include "FactoryBase.h"
#include "Trajectory.h"
#include <G4Trajectory.hh>
#include "DegradModel.h"

#include "PersistencyManager.h"

namespace nexus {
    REGISTER_CLASS(GarfieldEventAction,G4UserEventAction)
    GarfieldEventAction::GarfieldEventAction():G4UserEventAction() , nevt_(0), nupdate_(10), energy_min_(0.), energy_max_(DBL_MAX)
    {
        msg_ = new G4GenericMessenger(this, "/Actions/GarfieldEventAction/");

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

    GarfieldEventAction::~GarfieldEventAction() {
        G4cout << "Deleting EventAction" << G4endl;
    }


    void GarfieldEventAction::BeginOfEventAction(const G4Event *ev) {
      G4cout << " GarfieldEventAction::BeginOfEventAction()  0 " << G4endl;
        DegradModel* dm = (DegradModel*)(G4GlobalFastSimulationManager::GetInstance()->GetFastSimulationModel("DegradModel"));
        if(dm)
            dm->Reset();
        // Print out event number info
        if ((nevt_ % nupdate_) == 0) {
            G4cout << " >> Event no. " << nevt_  << G4endl;
            if (nevt_  == (10 * nupdate_)) nupdate_ *= 10;
        }
    }

    void GarfieldEventAction::EndOfEventAction(const G4Event *evt) {

        GarfieldVUVPhotonModel* gvm = (GarfieldVUVPhotonModel*)(G4GlobalFastSimulationManager::GetInstance()->GetFastSimulationModel("GarfieldVUVPhotonModel"));
        if(gvm)
          gvm->Reset(); // zero out the sensor: meaning reset the nexcitations, which is cumulative.

        G4cout << " EventAction::EndOfEventAction()  " << G4endl;



        // Get the trajectories stored for this event and loop through them
        // to calculate the total energy deposit

        G4double fEDepPrim = 0.;

        G4TrajectoryContainer* tc = evt->GetTrajectoryContainer();
        if (tc) {
            // in interactive mode, a G4TrajectoryContainer would exist
            // but the trajectories will not cast to Trajectory
            Trajectory* trj = dynamic_cast<Trajectory*>((*tc)[0]);
            if (trj == nullptr){
                G4Exception("[GarfieldEventAction]", "EndOfEventAction()", FatalException,
                            "GarfieldEventAction is required when using GarfieldEventAction");
            }
            for (unsigned int i=0; i<tc->size(); ++i) {
                Trajectory* tr = dynamic_cast<Trajectory*>((*tc)[i]);
                fEDepPrim += tr->GetEnergyDeposit();
            }
        }
        else {
            G4Exception("[DefaultEventAction]", "EndOfEventAction()", FatalException,
                        "GarfieldEventAction is required when using GarfieldEventAction");
        }

        PersistencyManager* pm = dynamic_cast<PersistencyManager*>
        (G4VPersistencyManager::GetPersistencyManager());

        if (!evt->IsAborted() && fEDepPrim>0) {
            pm->InteractingEvent(true);
        } else {
            pm->InteractingEvent(false);
        }
        if (!evt->IsAborted() && fEDepPrim > energy_min_ && fEDepPrim < energy_max_) {
            pm->StoreCurrentEvent(true);
        } else {
            pm->StoreCurrentEvent(false);
        }

    }

    void GarfieldEventAction::EDepPrim(const G4double &Ed)
    {
        fEDepPrim+=Ed;
    }

}
