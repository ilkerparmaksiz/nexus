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
        G4ThreeVector GetTrackEndPoint(G4int trk_id);
        G4double GetTrackEndTime(G4int trk_id);
        G4int GetCurrentTrackIndex(G4int trk_id); // Get the index of the track in the vector
        G4double GetAvgIoniEnergy(G4int trk_id);
        G4int GetTotIonizations(G4int trk_id);
        G4double GetTrackLength(G4int trk_id);
        G4int GetBremID(G4int trk_id, G4int hit_id);
        ~DegradModel();


        virtual G4bool IsApplicable(const G4ParticleDefinition&);
        virtual G4bool ModelTrigger(const G4FastTrack&);
        virtual void DoIt(const G4FastTrack&, G4FastStep&);
        void Reset();
        void AddTrack(G4int trk_id);

        private:
        void GetElectronsFromDegrad(G4FastStep& fastStep,G4ThreeVector degradPos,G4double degradTime, G4int trk_id);
        void SetTrackEndPoint(G4ThreeVector pos, G4double time, G4int trk_index);
        void SetNioni(G4int Ne, G4int trk_index); // Set the number of ionization particles
        void AddTrackLength(G4int trk_id); // Get track length
        void AddBremID(G4int trk_index, G4int brem_id); // Add the brem id to a vector
        G4double GetScintTime(); // Get timing delay from scintillaiton

        G4double fPrimKE; // Primary kinetic energy of the particle

        GarfieldHelper GH_;

        G4int event_id_;
        std::vector<G4double>      end_times_;
        std::vector<G4int>         N_ioni_;
        std::vector<G4double>      ke_vec_;
        std::vector<G4double>      trk_len_vec_;
        std::vector<G4ThreeVector> track_end_pos_;
        std::vector<std::vector<G4int>> brem_id_vec_;
        G4bool degrad_status_; // Checks if degrad has been run

        std::vector<G4int> track_ids_;

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
