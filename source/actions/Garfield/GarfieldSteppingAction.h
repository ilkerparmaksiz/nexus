#ifndef GarfieldSteppingAction_h
#define GarfieldSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4Types.hh"
#include "G4String.hh"
#include "G4OpBoundaryProcess.hh"
#include "GarfieldEventAction.h"
#include "G4GenericMessenger.hh"
#include "GarfieldEventAction.h"
#include "Analysis.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace nexus {
    class GarfieldSteppingAction : public G4UserSteppingAction {
     public:
        GarfieldSteppingAction(GarfieldEventAction *eva);
        ~GarfieldSteppingAction();

      void UserSteppingAction(const G4Step *);

     private:
        GarfieldEventAction* fEventAction;
        G4GenericMessenger* msg_;

        G4int ev_shift;

        G4bool reflected;
        G4int trackID = -99;
        G4String Material_Store;

    };
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
