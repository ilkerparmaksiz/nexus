//
// Created by ilker on 3/24/23.
//

#ifndef CRAB_OPWLS_HH
#define CRAB_OPWLS_HH
#include "G4OpWLS.hh"
#include "S2Photon.hh"
namespace nexus{
    class OpWLS:public G4OpWLS{
    public:
       using G4OpWLS::G4OpWLS;
       virtual G4bool IsApplicable(const G4ParticleDefinition &aParticleType) override;
    };

    #endif //CRAB_OPWLS_HH

    inline G4bool OpWLS::IsApplicable(const G4ParticleDefinition &aParticleType) {
        if(&aParticleType==G4OpticalPhoton::OpticalPhoton() || &aParticleType==S2Photon::OpticalPhoton()){
            return true;
        }
        return false;

    }
}