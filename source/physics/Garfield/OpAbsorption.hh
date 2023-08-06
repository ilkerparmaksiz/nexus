//
// Created by ilker on 3/24/23.
//

#ifndef CRAB_OPABSORPTION_HH
#define CRAB_OPABSORPTION_HH

#include "G4OpAbsorption.hh"
#include "S2Photon.hh"
namespace nexus{
    class OpAbsorption : public G4OpAbsorption{
    public:
         using G4OpAbsorption::G4OpAbsorption;
         virtual G4bool IsApplicable(const G4ParticleDefinition &aParticleType) override ;
    };


    #endif //CRAB_OPABSORPTION_HH

    inline G4bool OpAbsorption::IsApplicable(const G4ParticleDefinition &aParticleType) {
        if(&aParticleType==G4OpticalPhoton::OpticalPhoton() || &aParticleType==S2Photon::OpticalPhoton()){
            return true;
        }
        return false;

    }
}