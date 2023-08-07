// ----------------------------------------------------------------------------
// nexus | IonizationElectron.h
//
// Definition of the ionization electron.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#ifndef S2_PHOTON_H
#define S2_PHOTON_H
#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include <G4ParticleDefinition.hh>



namespace nexus{
  class S2Photon: public G4ParticleDefinition
  {
  public:
    /// Returns a pointer to the only instance
    static S2Photon* Definition();

    static S2Photon* OpticalPhotonDefinition();
    static S2Photon* OpticalPhoton();

    /// Destructor
    ~S2Photon(){};

  private:
    /// Default constructor is hidden.
    /// No instance of this class can be created.
    S2Photon(){};

  private:
      static S2Photon* theInstance;
  };

}
#endif
