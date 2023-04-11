//
// Created by ilker on 7/19/22.
//

#include "CRABOpticalTrackingAction.h"
// ----------------------------------------------------------------------------
// nexus | OpticalTrackingAction.cc
//
// This class saves the trajectories of optical photons, in addition to the
// particles saved by the default tracking action. Its purpose is to store
// optical photon information in the output file.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "OpticalTrackingAction.h"

#include "Trajectory.h"
#include "TrajectoryMap.h"
#include "FactoryBase.h"

#include <G4Track.hh>
#include <G4TrackingManager.hh>
#include <G4Trajectory.hh>
#include <G4OpticalPhoton.hh>
#include "IonizationElectron.h"
#include <G4GenericMessenger.hh>
#include "G4OpticalPhoton.hh"

using namespace nexus;

REGISTER_CLASS(CRABOpticalTrackingAction, G4UserTrackingAction)

CRABOpticalTrackingAction::CRABOpticalTrackingAction(): G4UserTrackingAction(),Particle("optical")
{
    msg_ = new G4GenericMessenger(this, "/TrackingAction/CRABOpticalTrackingAction/");

    G4GenericMessenger::Command& thresh_cmd =
            msg_->DeclareProperty("Particle_Type", Particle,
                                  "Select which secondary particle needs to be saved to file (opticalphoton,ionizationelectron(-ie),or All");
}



CRABOpticalTrackingAction::~CRABOpticalTrackingAction()
{
}



void CRABOpticalTrackingAction::PreUserTrackingAction(const G4Track* track)
{
    // Create a new trajectory associated to the track.
    // N.B. If the processesing of a track is interrupted to be resumed
    // later on (to process, for instance, its secondaries) more than
    // one trajectory associated to the track will be created, but
    // the event manager will merge them at some point.
    G4VTrajectory* trj = new Trajectory(track);


    if (Particle=="ie") //Skip Optical Photons and Keep Ionization Electrons
    {
        if(track->GetDefinition() == IonizationElectron::Definition())
            fpTrackingManager->SetStoreTrajectory(true);

        else if(track->GetDefinition() == G4OpticalPhoton::Definition())
            fpTrackingManager->SetStoreTrajectory(false);
    }
    else if(Particle=="optical"){ //Skip Ionization Electrons and Keep Optical Photons
        if(track->GetDefinition() == IonizationElectron::Definition())
            fpTrackingManager->SetStoreTrajectory(false);
        else if(track->GetDefinition() == G4OpticalPhoton::Definition())
            fpTrackingManager->SetStoreTrajectory(true);



    }else // Keep ALL Be aware file with get large with this setting
        fpTrackingManager->SetStoreTrajectory(true); // Set the trajectory in the tracking manager


    fpTrackingManager->SetTrajectory(trj);

}



void CRABOpticalTrackingAction::PostUserTrackingAction(const G4Track* track)
{
    Trajectory* trj = (Trajectory*) TrajectoryMap::Get(track->GetTrackID());

    // Do nothing if the track has no associated trajectory in the map
    if (!trj) return;


    if (track->GetDefinition() == G4OpticalPhoton::Definition() && Particle=="ie") //Skip Optical Photons and Keep Ionization Electrons
            return;
    else if(track->GetDefinition() == IonizationElectron::Definition() && Particle=="optical")  //Skip Ionization Electrons and Keep Optical Photons
            return;

    // Record final time and position of the track
    trj->SetFinalPosition(track->GetPosition());
    trj->SetFinalTime(track->GetGlobalTime());
    trj->SetTrackLength(track->GetTrackLength());
    trj->SetFinalMomentum(track->GetMomentum());

    // In case of optical photons
    if (track->GetDefinition() == G4OpticalPhoton::Definition()) {
        // If optical-photon has no NextVolume (escaping from the world)
        // Assign current volume as the decay one
        if (track->GetNextVolume()) trj->SetFinalVolume(track->GetNextVolume()->GetName());
        else                        trj->SetFinalVolume(track->GetVolume()->GetName());
    }
        // Final Volume of non optical photons
    else trj->SetFinalVolume(track->GetVolume()->GetName());

    // Record last process of the track
    G4String final_process = track->GetStep()->GetPostStepPoint()
            ->GetProcessDefinedStep()->GetProcessName();
    trj->SetFinalProcess(final_process);
}
