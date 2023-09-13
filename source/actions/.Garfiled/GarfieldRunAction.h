#ifndef RunAction_hh
#define RunAction_hh 1



#include "TROOT.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TFile.h"

// using namespace std;         //needed for using standard libraries

#include <G4UserRunAction.hh>

// Run action class, carries out tasks at the begin and end of each run.
// The concept of a run incorporates a fixed geometry, fixed beam conditions,
// simulation of number of primaries.
// Begins with /run/beamOn command and finishes with tracking of last secondary
// to zero energy:
class GarfieldRunAction : public G4UserRunAction {
 public:
  // Run action class needs a pointer of the detector construction class in
  // order to get details of the readout geometry.
  // Accepts pointer to detector construction class:
  GarfieldRunAction();
  ~GarfieldRunAction();

  void BeginOfRunAction(const G4Run *);
  void EndOfRunAction(const G4Run *);


 private:

};
#endif
