//
// Created by ilker on 9/2/21.
//

#ifndef NEXUS_GRAB_H
#define NEXUS_GRAB_H
#include "GeometryBase.h"

class G4GenericMessenger;

namespace nexus {

    class CRAB: public GeometryBase
    {
        public:
            /// Constructor
            CRAB();
            /// Destructor
            ~CRAB();

            /// Return vertex within region <region> of the chamber

            virtual void Construct();
            virtual G4ThreeVector GenerateVertex(const G4String& region) const;

        private:
            /// Messenger for the definition of control commands
            G4GenericMessenger* msg_;
            G4double Lab_size;
            G4double chamber_diam   ;
            G4double chamber_length ;
            G4double chamber_thickn ;
            G4double SourceEn_diam   ;
            G4double SourceEn_length ;
            G4double SourceEn_thickn ;
            G4double SourceEn_holedia ;
            G4double gas_pressure_;
            G4ThreeVector vtx_;
            G4ThreeVector vertex;
            void ConstructLab();
            void PlaceVolumes();
            void AssignVisuals();
            void PrintParam();

    };

} // end namespace nexus


#endif //NEXUS_GRAB_H
