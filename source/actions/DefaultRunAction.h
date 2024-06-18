// ----------------------------------------------------------------------------
// nexus | DefaultRunAction.h
//
// This is the default run action of the NEXT simulations.
// A message at the beginning and at the end of the simulation is printed.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#ifndef DEFAULT_RUN_ACTION_H
#define DEFAULT_RUN_ACTION_H

#include <G4UserRunAction.hh>
#include "DefaultRunAction.h"
#include "chrono"
#include "vector"

namespace nexus {

  class DefaultRunAction: public G4UserRunAction
  {
  public:
    /// Constructor
    DefaultRunAction();
    /// Destructor
    ~DefaultRunAction();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);
    std::chrono::time_point<std::chrono::high_resolution_clock> startTime;
    std::chrono::time_point<std::chrono::high_resolution_clock> endTime;
    std::vector<double>* EventCompletionTime;
    std::vector<double>* NumberofPhotonsPerEvent;
    double Runtime;

  };

}

#endif
