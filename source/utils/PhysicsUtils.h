// ----------------------------------------------------------------------------
// nexus | Physics Utils.h
//
// Commonly used functions by physics classes
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include <Randomize.hh>
#include <G4PhysicsTable.hh>
#include <G4PhysicsOrderedFreeVector.hh>
#include <G4ThreeVector.hh>


#ifndef PHYS_U_H
#define PHYS_U_H

using namespace CLHEP;

namespace nexus {

    void BuildThePhysicsTable( G4PhysicsTable* &theFastIntegralTable);
    void ComputeCumulativeDistribution(const G4PhysicsOrderedFreeVector&, G4PhysicsOrderedFreeVector&);

    // Get the photon polarization vector
    void GetPhotonPol(G4ThreeVector &momentum, G4ThreeVector &polarization);


}

#endif
