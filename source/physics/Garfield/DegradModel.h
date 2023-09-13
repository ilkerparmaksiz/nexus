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
#include "IonizationSD.h"


class G4VPhysicalVolume;
class G4GenericMessenger;

namespace nexus{
    class IonizationSD;
    class GeometryBase;
    class DegradModel : public G4VFastSimulationModel {
     public:
        //-------------------------
        // Constructor, destructor
        //-------------------------
        DegradModel(G4String, G4Region*,IonizationSD*);
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
        IonizationSD* fIonizationSD;
        G4bool processOccured;

        G4double GasPressure_;

        char* crab_path; // Path to the root directory

        // Messenger for the definition of control commands
        G4GenericMessenger* msg_;

    };
}
#endif /* DegradModel_H_ */
