#include "GarfieldTrackingAction.h"
#include "G4VPhysicalVolume.hh"
#include "G4VTouchable.hh"

#include "G4Track.hh"
#include "G4OpticalPhoton.hh"

#include "G4SDManager.hh"
#include "G4Run.hh"
#include "G4EventManager.hh"
#include "S2Photon.hh"
#include <G4Trajectory.hh>
#include "Trajectory.h"
#include "TrajectoryMap.h"
#include "FactoryBase.h"

using namespace nexus;
REGISTER_CLASS(GarfieldTrackingAction,G4UserTrackingAction)

GarfieldTrackingAction::GarfieldTrackingAction():G4UserTrackingAction(),fPPID(0.0),fPKE(0.0), fFPKE(0.0)
{
}

GarfieldTrackingAction::~GarfieldTrackingAction()
{
}
void GarfieldTrackingAction::PreUserTrackingAction(const G4Track *Track) {


    // have to include  NEST ThermalElectrons here
    // Do not Store Trajectories to help on memory
    if (Track->GetDefinition() == G4OpticalPhoton::Definition() ||
        Track->GetDefinition() == S2Photon::Definition()) {
        fpTrackingManager->SetStoreTrajectory(false);
        return;
    }
    // Create a new trajectory associated to the track.
    // N.B. If the processesing of a track is interrupted to be resumed
    // later on (to process, for instance, its secondaries) more than
    // one trajectory associated to the track will be created, but
    // the event manager will merge them at some point.
    G4VTrajectory *trj = new Trajectory(Track);

    // Set the trajectory in the tracking manager
    fpTrackingManager->SetStoreTrajectory(true);
    fpTrackingManager->SetTrajectory(trj);
}

void GarfieldTrackingAction::PostUserTrackingAction(const G4Track *Track)
{
    // Do nothing if the track is an optical photon or an S2 Photon
    if (Track->GetDefinition() == G4OpticalPhoton::Definition() ||
            Track->GetDefinition() == S2Photon::Definition())
        return;

    Trajectory *trj = (Trajectory *)TrajectoryMap::Get(Track->GetTrackID());

    // Do nothing if the track has no associated trajectory in the map
    if (!trj) return;

    // Record final time and position of the track
    trj->SetFinalPosition(Track->GetPosition());
    trj->SetFinalTime(Track->GetGlobalTime());
    trj->SetTrackLength(Track->GetTrackLength());
    trj->SetFinalVolume(Track->GetVolume()->GetName());
    trj->SetFinalMomentum(Track->GetMomentum());

    // Record last process of the track
    G4String proc_name = Track->GetStep()->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    trj->SetFinalProcess(proc_name);
}


