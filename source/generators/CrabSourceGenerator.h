//
// Created by ilker on 9/7/21.
// implemented  from(Austin) https://github.com/Q-Pix/Q_PIX_GEANT4/blob/UTAH/src/Supernova.cpp and NextIonGenerator
//

#ifndef NEXUS_CRABSOURCEGENERATOR_H
#define NEXUS_CRABSOURCEGENERATOR_H


#include <G4VPrimaryGenerator.hh>

class G4Event;
class G4GenericMessenger;
class G4ParticleDefinition;


namespace nexus{

    class GeometryBase;

    class CrabSourceGenerator: public G4VPrimaryGenerator
    {
    public:
        // Constructor
        CrabSourceGenerator();
        // Destructor
        ~CrabSourceGenerator();

        // This method is invoked at the beginning of the event,
        // setting a primary vertex that contains the chosen ion
        void GeneratePrimaryVertex(G4Event*);
        void EventsWithWindow(G4Event *,G4double decay_time);
    private:
        G4ParticleDefinition* IonDefinition();


    private:
        G4int atomic_number_, mass_number_;
        G4double energy_level_;
        G4bool decay_at_time_zero_;
        G4double decay_time_;
        G4double event_window_;
        G4int N_Decays_per_s_;
        G4String region_;
        G4GenericMessenger* msg_;
        const GeometryBase* geom_;
        G4int NDecays;

    };

} // end namespace nexus




#endif //NEXUS_CRABSOURCEGENERATOR_H
