//
// Created by ilker on 1/14/23.
//

#ifndef NEXUS_CRAB_EMCCD_H
#define NEXUS_CRAB_EMCCD_H
#include "GeometryBase.h"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolumeStore.hh"

class G4GenericMessenger;
namespace nexus {
    class CRAB_EMCCD : public GeometryBase{
    public:
        // Construct EMCCD_Class
        CRAB_EMCCD();

        // Destruct
        ~CRAB_EMCCD();

        // Define the volumes of the geometry
        virtual void Construct ();

        // Generate a vertex within a given region of the geometry
        virtual G4ThreeVector GenerateVertex(const G4String &region) const;

        void SetLabLogical(G4LogicalVolume *);
        G4LogicalVolume* GetLabLogical();
        void SetLabPhysical(G4VPhysicalVolume *);
        G4VPhysicalVolume* GetLabPhysical();
        G4double Offsetz_;
    private:
        G4GenericMessenger *msg_;
        bool fvisual;
        G4VPhysicalVolume * Lab_Physical;
        G4LogicalVolume * Lab_Logical;



    };

    inline void CRAB_EMCCD::SetLabLogical(G4LogicalVolume * lv) {
        Lab_Logical=lv;
    }
    inline void CRAB_EMCCD::SetLabPhysical(G4VPhysicalVolume * phys) {

        Lab_Physical=phys;
    }
    inline G4LogicalVolume* CRAB_EMCCD::GetLabLogical() {
        return Lab_Logical;
    }
    inline G4VPhysicalVolume* CRAB_EMCCD::GetLabPhysical() {
        return Lab_Physical;
    }
}
#endif //NEXUS_CRAB_EMCCD_H
