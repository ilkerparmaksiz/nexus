#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4Types.hh"
#include "G4String.hh"
#include "G4OpBoundaryProcess.hh"
#include "GarfieldEventAction.hh"
#include "Analysis.hh"
#include "G4GenericMessenger.hh"

#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


class SteppingAction : public G4UserSteppingAction {
 public:
  SteppingAction(EventAction *eva);
  ~SteppingAction(){};

  void UserSteppingAction(const G4Step *);
 
 private:
  EventAction* fEventAction;

  G4GenericMessenger* msg_;

  G4int ev_shift;

  G4bool reflected;
  G4int trackID = -99;
  G4String Material_Store;
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
