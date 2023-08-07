// ----------------------------------------------------------------------------
// nexus | GarfieldStackingAction.cc
//
// This class is an example of how to implement a stacking action, if needed.
// At the moment, it is not used in the NEXT simulations.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "GarfieldStackingAction.h"
#include "FactoryBase.h"


using namespace nexus;

REGISTER_CLASS(GarfieldStackingAction, NESTStackingAction)

GarfieldStackingAction::GarfieldStackingAction(): NESTStackingAction()
{
}
GarfieldStackingAction::~GarfieldStackingAction()  {}


