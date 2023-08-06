#ifndef TrackingAction_h
#define TrackingAction_h 1

#include "G4UserTrackingAction.hh"
#include "G4Types.hh"
#include "G4String.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4Event.hh"
#include "G4PrimaryVertex.hh"
#include "Analysis.hh"

#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


class TrackingAction : public G4UserTrackingAction {
 public:
  TrackingAction();
  ~TrackingAction(){};

  virtual void PreUserTrackingAction(const G4Track *);
  //  virtual void PostUserTrackingAction(const G4Track *);
 
 private:
  G4int fPPID;
  G4double fPKE;
  G4double fFPKE;
 	
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
