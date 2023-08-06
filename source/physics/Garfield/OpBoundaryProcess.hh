//
// Created by ilker on 3/24/23.
//

#ifndef CRAB_OPBOUNDARYPROCESS_HH
#define CRAB_OPBOUNDARYPROCESS_HH

#include "G4OpBoundaryProcess.hh"
#include "S2Photon.hh"
namespace nexus {
    class OpBoundaryProcess : public G4OpBoundaryProcess{
        using G4OpBoundaryProcess::G4OpBoundaryProcess;
        virtual G4bool IsApplicable(const G4ParticleDefinition& aParticleType) override;
        virtual void SetInvokeSD(G4bool) override;
    };

    #endif //CRAB_OPBOUNDARYPROCESS_HH

    inline G4bool OpBoundaryProcess::IsApplicable(const G4ParticleDefinition &aParticleType) {
        if(&aParticleType==G4OpticalPhoton::OpticalPhoton() || &aParticleType==S2Photon::OpticalPhoton()){
            return true;
        }
        return false;

    }

    inline void OpBoundaryProcess::SetInvokeSD(G4bool){
        // G4OpBoundaryProcess::SetInvokeSD(flag);
        return;
    }
}