/*
 * DegradModel.h
 *
 *  Created on: Apr 9, 2014
 *      Author: dpfeiffe
 */

#ifndef DEGRADMODEL_H_
#define DEGRADMODEL_H_

#include "GarfieldHelper.h"

#include <G4VFastSimulationModel.hh>
#include <G4PhysicsTable.hh>
#include <G4PhysicsOrderedFreeVector.hh>

namespace nexus{
    class DegradModel : public G4VFastSimulationModel {
     public:
        //-------------------------
        // Constructor, destructor
        //-------------------------
        DegradModel(G4String, G4Region*, GarfieldHelper GH);
        G4ThreeVector GetTrackEndPoint();
        G4double GetTrackEndTime();
        ~DegradModel();


        virtual G4bool IsApplicable(const G4ParticleDefinition&);
        virtual G4bool ModelTrigger(const G4FastTrack&);
        virtual void DoIt(const G4FastTrack&, G4FastStep&);
        inline G4bool FindParticleName(G4String s){if(s=="e-") return true; return false;};
        void SetPrimaryKE(G4double KE) {fPrimKE = KE;};

        private:
        void GetElectronsFromDegrad(G4FastStep& fastStep,G4ThreeVector degradPos,G4double degradTime);
        void SetTrackEndPoint(G4ThreeVector pos, G4double time);

        G4double fPrimKE; // Primary kinetic energy of the particle

        GarfieldHelper GH_;

        G4double      end_time;
        G4ThreeVector track_end_pos;

        // Scintillation timing 
        G4double slow_comp_; // ns
        G4double slow_prob_; // %
        G4double fast_comp_; // ns
        G4double fast_prob_; // %

        G4PhysicsTable* theFastIntegralTable_;
        G4PhysicsOrderedFreeVector* spectrum_integral;

    };
}
#endif /* DegradModel_H_ */
