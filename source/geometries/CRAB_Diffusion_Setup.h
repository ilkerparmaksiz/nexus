//
// Created by ilker on 9/2/21.
//

#ifndef NEXUS_CRAB_Diffusion_Setup_H
#define NEXUS_CRAB_Diffusion_Setup_H

#include "GeometryBase.h"
#include "PmtR7378A.h"
#include "CylinderPointSampler2020.h"
#include "CRAB_Detector.h"
#include "CRAB_EMCCD.h"
#include "CRABImageIntensifier.h"
class G4GenericMessenger;
namespace nexus {

    class CylinderPointSampler2020;
    class CRAB_Diffusion_Setup: public GeometryBase
    {
        public:

            /// Constructor
            CRAB_Diffusion_Setup();

            /// Destructor
            ~CRAB_Diffusion_Setup();

            /// Return vertex within region <region> of the chamber
            virtual void Construct();
            //virtual G4ThreeVector GenerateVertex(const G4String& region) const;

        private:

            /// Messenger for the definition of control commands
            G4GenericMessenger* msg_;

            // Importing Geometries
            CRAB_Detector * det;
            PmtR7378A * AnodePMT;
            CRAB_EMCCD * theCamera;
            CRABImageIntensifier * II;

            void PlaceVolumes();
            void AssignVisuals();
            void PrintParam();


    };


} // end namespace nexus


#endif //NEXUS_CRAB_Diffusion_Setup_H
