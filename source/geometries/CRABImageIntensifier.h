//
// Created by ilker on 1/14/23.
//

#ifndef NEXUS_CRABIMAGEINTENSIFIER_H
#define NEXUS_CRABIMAGEINTENSIFIER_H

#include "GeometryBase.h"
#include "G4VPhysicalVolume.hh"
class G4GenericMessenger;
namespace  nexus {


    class CRABImageIntensifier: public GeometryBase{

    public:
        // Construct
        CRABImageIntensifier();


        // Destruct Geometry
        ~CRABImageIntensifier();


        /// Define Volumes of the geometry
        virtual void Construct();
        // Generate a vertex within a given region of geometry
       virtual G4ThreeVector GenerateVertex(const G4String &region) const;


       //
       void SetLabLogical(G4LogicalVolume *);
        G4LogicalVolume* GetLabLogical();
        void SetLabPhysical(G4VPhysicalVolume *);
        G4VPhysicalVolume* GetLabPhysical();
        G4double Offsetz_;


    private:
        G4GenericMessenger * msg_;
        long PhotonGain;
        bool fvisual;
        G4VPhysicalVolume * Lab_Physical;
        G4LogicalVolume * Lab_Logical;

    };

    inline void CRABImageIntensifier::SetLabLogical(G4LogicalVolume * lv) {
        Lab_Logical=lv;
    }
    inline void CRABImageIntensifier::SetLabPhysical(G4VPhysicalVolume * phys) {

        Lab_Physical=phys;
    }
    inline G4LogicalVolume* CRABImageIntensifier::GetLabLogical() {
        return Lab_Logical;
    }
    inline G4VPhysicalVolume* CRABImageIntensifier::GetLabPhysical() {
        return Lab_Physical;
    }

} // end namespace nexus

#endif //NEXUS_CRABIMAGEINTENSIFIER_H
