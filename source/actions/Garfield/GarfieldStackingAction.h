// ----------------------------------------------------------------------------
// nexus | GarfieldStackingAction.h
//
// This class is an example of how to implement a stacking action, if needed.
// At the moment, it is not used in the NEXT simulations.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#ifndef DEFAULT_STACKING_ACTION_H
#define DEFAULT_STACKING_ACTION_H
#include "NEST/G4/NESTStackingAction.hh"
#include <G4UserStackingAction.hh>


namespace nexus {

  // General-purpose user stacking action

  class GarfieldStackingAction: public G4UserStackingAction
  {
  public:
    /// Constructor
    GarfieldStackingAction();
    /// Destructor
    ~GarfieldStackingAction();

  };

} // end namespace nexus

#endif
