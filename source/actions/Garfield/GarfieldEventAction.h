#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"
#include <vector>

class G4VPhysicalVolume;
class SteppingAction;
class G4Event;
class G4GenericMessenger;

namespace nexus {
    class GarfieldEventAction : public G4UserEventAction {
     public:
        GarfieldEventAction();
      ~GarfieldEventAction();

     public:
      void BeginOfEventAction(const G4Event *);
      void EndOfEventAction(const G4Event *);
      void EDepPrim(const G4double&);

    private:
        G4GenericMessenger* msg_;
        G4int nevt_, nupdate_;
        G4double energy_min_;
        G4double energy_max_;
        G4double fEDepPrim;



    };

}
#endif
