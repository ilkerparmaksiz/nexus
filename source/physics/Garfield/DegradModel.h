/*
 * DegradModel.h
 *
 *  Created on: Apr 9, 2014
 *      Author: dpfeiffe
 */

#ifndef DEGRADMODEL_H_
#define DEGRADMODEL_H_

#include "G4ThreeVector.hh"
#include "G4VFastSimulationModel.hh"
#include "GasModelParameters.h"
#include "GasBoxSD.h"

class G4VPhysicalVolume;

namespace nexus{
    class GasBoxSD;
    class DetectorConstruction;
    class DegradModel : public G4VFastSimulationModel {
     public:
        //-------------------------
        // Constructor, destructor
        //-------------------------
        DegradModel(GasModelParameters*, G4String, G4Region*,DetectorConstruction*,GasBoxSD*);
        ~DegradModel();


        virtual G4bool IsApplicable(const G4ParticleDefinition&);
        virtual G4bool ModelTrigger(const G4FastTrack&);
        virtual void DoIt(const G4FastTrack&, G4FastStep&);
        inline G4bool FindParticleName(G4String s){if(s=="e-") return true; return false;};
        inline void Reset(){processOccured=false;};
        void SetPrimaryKE(G4double KE) {fPrimPhotonKE = KE;};

        private:
        void GetElectronsFromDegrad(G4FastStep& fastStep,G4ThreeVector degradPos,G4double degradTime);


        G4double thermalE;
        G4double fPrimPhotonKE;
        DetectorConstruction* detCon;
        GasBoxSD* fGasBoxSD;
        G4bool processOccured;

        char* crab_path; // Path to the root directory


    };
}
#endif /* DegradModel_H_ */
