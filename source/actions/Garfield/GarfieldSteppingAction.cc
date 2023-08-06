#include "GarfieldSteppingAction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VTouchable.hh"
#include "G4Step.hh"
#include "DetectorConstruction.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "G4StepPoint.hh"
#include "G4SDManager.hh"
#include "G4Run.hh"
#include "G4EventManager.hh"
#include "GasBoxSD.hh"
#include "S2Photon.hh"


SteppingAction::SteppingAction(EventAction *eva) : fEventAction(eva),  ev_shift(0) {

  msg_ = new G4GenericMessenger(this, "/Action/SteppingAction/",
    "Control commands of the stepping action.");

  msg_->DeclareProperty("event_shift", ev_shift, "Set the event ID Shift number");
  
}



void SteppingAction::UserSteppingAction(const G4Step *aStep)

{

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4int  event = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
  const G4ThreeVector pos(aStep->GetPreStepPoint()->GetPosition());
  const G4ThreeVector tpos(aStep->GetPostStepPoint()->GetPosition());
  const G4ThreeVector sdir(aStep->GetPreStepPoint()->GetMomentumDirection());
  const G4ThreeVector tdir(aStep->GetPostStepPoint()->GetMomentumDirection());

  
  std::string startp("null");
  std::string endp("null");

  const G4StepPoint* endPoint = aStep->GetPostStepPoint();

  const G4VProcess* tprocess   = endPoint->GetProcessDefinedStep();
  const G4VProcess* sprocess   = aStep->GetPreStepPoint()->GetProcessDefinedStep();

  tprocess = aStep->GetPostStepPoint()->GetProcessDefinedStep();

  G4Track* track = aStep->GetTrack();
  G4double tID = track->GetTrackID();

  if (trackID != track->GetTrackID()){
    trackID =  track->GetTrackID();
    reflected = false;
  }

  const G4ParticleDefinition* particle = track->GetParticleDefinition();
  G4int pID       = particle->GetPDGEncoding();
  G4double time   = aStep->GetPreStepPoint()->GetGlobalTime();

  G4OpBoundaryProcess* boundary = 0;
  // if (!boundary &&  particle->GetParticleName() == "S2Photon") {
  if (!boundary){
      
      G4ProcessVector* pv = particle->GetProcessManager()->GetProcessList();
      for (size_t i=0; i<pv->size(); i++) {
          
          if ((*pv)[i]->GetProcessName() == "OpBoundary") {
              boundary = (G4OpBoundaryProcess*) (*pv)[i];
              
              // This is a reflection on the SS
              if(aStep->GetPreStepPoint()->GetMaterial()->GetName() == "GXe" && aStep->GetPostStepPoint()->GetMaterial()->GetName() == "Steel" && boundary->GetStatus() == SpikeReflection){
                  reflected = true;
               }

              // This gives the material that we just reflected off
              if(aStep->GetPreStepPoint()->GetMaterial()->GetName() == "Steel" && aStep->GetPostStepPoint()->GetMaterial()->GetName() == "GXe" && boundary->GetStatus() == StepTooSmall){
                Material_Store = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName();
              }

              // This gives the material that we just reflected off
              if(aStep->GetPostStepPoint()->GetMaterial()->GetName() == "MgF2" && boundary->GetStatus() == FresnelReflection){
                reflected = true;
                Material_Store = aStep->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName();
              }

              // Total internal reflection from MgF2
              if(aStep->GetPreStepPoint()->GetMaterial()->GetName() == "MgF2" && boundary->GetStatus() == TotalInternalReflection){
                reflected = true;
                Material_Store = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName();
              }

              break;
          }
      }
  }

  G4int id(0);

  
  if (sprocess)
      startp = sprocess->GetProcessName();
  if (tprocess)
      endp = tprocess->GetProcessName();

  
  G4LogicalVolume* lVolume = aStep->GetPreStepPoint()->GetTouchableHandle()
                             ->GetVolume()->GetLogicalVolume();

   
  G4TouchableHandle touch = endPoint->GetTouchableHandle();
  G4VPhysicalVolume* eVolume = touch->GetVolume();
  G4String eVname("null");

  if (pID==11 && track->GetKineticEnergy()/keV>0.100 && (lVolume->GetName().find("GAS")!=std::string::npos)) // don't count the thermale's, just G4 e's
    fEventAction->EDepPrim(aStep->GetTotalEnergyDeposit());


  // Here we select if we got an S1 photon or S2 photon
  // S1  = 1, S2 = 2
  G4int PhotonType = 2;
  if(track->GetParticleDefinition()==S2Photon::OpticalPhoton()){
    PhotonType = 2;
  }
  else if(track->GetParticleDefinition()==G4OpticalPhoton::OpticalPhoton()){
    PhotonType = 1;
  }
 
 

  
  if (eVolume){
    
    eVname = eVolume->GetName();

    // Camera
    if (lVolume->GetName().find("camLogical")!=std::string::npos){
      track->SetTrackStatus(fStopAndKill);
      analysisManager->FillNtupleDColumn(id,0, event+ev_shift);
      analysisManager->FillNtupleDColumn(id,1, pID);
      analysisManager->FillNtupleDColumn(id,2, time/ns);
      analysisManager->FillNtupleDColumn(id,3, pos[0]/mm);
      analysisManager->FillNtupleDColumn(id,4, pos[1]/mm);
      analysisManager->FillNtupleDColumn(id,5, pos[2]/mm);
      
      // Reflected photon
      if (reflected){
        analysisManager->FillNtupleIColumn(id,6, 1);
        analysisManager->FillNtupleSColumn(id,7, Material_Store);
      }
      else {
        analysisManager->FillNtupleIColumn(id,6, 0);
        analysisManager->FillNtupleSColumn(id,7, "None");
      }

      analysisManager->FillNtupleIColumn(id,8, PhotonType);
    
      analysisManager->AddNtupleRow(id);

      // if (reflected) std::cout << "Parent ID from reflected photon Detected: " << track->GetTrackID() << "  Material:  " << Material_Store << std::endl;
      // else  std::cout << "Photon arrived but was not reflected: " << track->GetTrackID() << std::endl;

    }
    

    // Lens
    id = 4;
    if (lVolume->GetName().find("Lens")!=std::string::npos)  {

        // analysisManager->FillNtupleDColumn(id,0, event+ev_shift);
        // analysisManager->FillNtupleDColumn(id,1, pID);
        // analysisManager->FillNtupleDColumn(id,2, time/ns);
        // analysisManager->FillNtupleDColumn(id,3, pos[0]/mm);
        // analysisManager->FillNtupleDColumn(id,4, pos[1]/mm);
        // analysisManager->FillNtupleDColumn(id,5, pos[2]/mm);
        // analysisManager->AddNtupleRow(id);

	  }

    // PMT location
    id = 5;
    if (lVolume->GetName().find("S1_WINDOW")!=std::string::npos){
      analysisManager->FillNtupleDColumn(id,0, event+ev_shift);
      analysisManager->FillNtupleDColumn(id,1, pID);
      analysisManager->FillNtupleDColumn(id,2, time/ns);
      analysisManager->FillNtupleDColumn(id,3, pos[0]/mm);
      analysisManager->FillNtupleDColumn(id,4, pos[1]/mm);
      analysisManager->FillNtupleDColumn(id,5, pos[2]/mm);	  
      
      // Reflected photon
      if (reflected){
        analysisManager->FillNtupleIColumn(id,6, 1);
        analysisManager->FillNtupleSColumn(id,7, Material_Store);
      }
      else {
        analysisManager->FillNtupleIColumn(id,6, 0);
        analysisManager->FillNtupleSColumn(id,7, "None");
      }

      analysisManager->FillNtupleIColumn(id,8, PhotonType);

      analysisManager->AddNtupleRow(id);
    
    }

  }
  
  // if particle == thermale, opticalphoton and parent == primary and stepID==1, or trackID<=2
  // count the NEST e-s/photons into a class variable from the primary particle. Retrieve at EndEvent().
}



  
