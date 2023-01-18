//
// Created by ilker on 1/14/23.
//

#include "CRAB_EMCCD.h"
#include "FactoryBase.h"
#include "G4Tubs.hh"
#include "G4GenericMessenger.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "MaterialsList.h"
#include "G4NistManager.hh"
#include "G4MaterialTable.hh"
#include "G4PVPlacement.hh"
#include "Visibilities.h"
namespace nexus {

    using namespace CLHEP;

    // Register the classs file
    REGISTER_CLASS(CRAB_EMCCD,GeometryBase)
    //Construct
    CRAB_EMCCD::CRAB_EMCCD() :GeometryBase(),msg_(nullptr),fvisual(true){
        msg_=new G4GenericMessenger(this,"/Geometry/CRAB_EMMCD","This is the geometry of the EMCCD");

    }

    // Destruct
    CRAB_EMCCD::~CRAB_EMCCD() noexcept {}

    void CRAB_EMCCD::Construct() {
        // C91000-23B Camera Pixel Size
        float pixelSize=16*um;
        float effectiveXandY=8.192*mm;

        // LAB
        G4String lab_name = "Lab";
        G4Box * lab_solid_volume = new G4Box(lab_name,0.4/2*m,0.4/2*m,0.4/2*m);
        G4LogicalVolume * lab_logic_volume= new G4LogicalVolume(lab_solid_volume,G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"),lab_name) ;
        auto labPhysical= new G4PVPlacement(0,G4ThreeVector(),lab_logic_volume,lab_logic_volume->GetName(),0,false,0, false);


        // Solids
        G4Box * CamBody_solid=new G4Box("CamBody",(110/2)*mm,(110/2)*mm,(156/2)*mm);

        G4Tubs * CamFront_solid=new G4Tubs("CamFront",0,(76/2)*mm,(25/2)*mm,0,twopi);

        G4Box * PixelLocation_solid= new G4Box("PixelLocation",effectiveXandY/2,effectiveXandY/2,1*mm);

        G4Box * Pixel_solid= new G4Box("Pixel",pixelSize/2,pixelSize/2,pixelSize/4);

        G4Material *Alum=  G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
        G4Material *Si=G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");

        // Logical Volume
        G4LogicalVolume * CamBody_logic=new G4LogicalVolume(CamBody_solid,Alum,"CamBody_Logic");

        G4LogicalVolume * CamFront_logic=new G4LogicalVolume(CamFront_solid,Alum,"CamFront_Logic");

        G4LogicalVolume * EffectiveArea_logic=new G4LogicalVolume(PixelLocation_solid,Alum,"EffectiveArea_Logic");

        G4LogicalVolume * Pixel_Logic=new G4LogicalVolume(Pixel_solid,Si,"Pixel");


        // Place Solids

        new G4PVPlacement(0,G4ThreeVector(0,0,0),CamBody_solid->GetName(),CamBody_logic,labPhysical,false,0,false);
        new G4PVPlacement(0,G4ThreeVector(0,0,(-156/2-25/2)*mm),CamFront_logic->GetName(),CamFront_logic,labPhysical,false,0,false);
        new G4PVPlacement(0,G4ThreeVector(0,0,(-156/2-25/2-1)*mm),EffectiveArea_logic->GetName(),EffectiveArea_logic,labPhysical,false,0,false);
        new G4PVPlacement(0,G4ThreeVector(0,0,(-156/2-25/2-1-pixelSize/4)*mm),Pixel_Logic->GetName(),Pixel_Logic,labPhysical,true,0,false);

        // Visuals
        if (fvisual){
            G4VisAttributes * Camfront = new G4VisAttributes(nexus::DarkRed());
            G4VisAttributes * CamBody = new G4VisAttributes(nexus::DarkGrey());
            G4VisAttributes * EffectiveArea = new G4VisAttributes(nexus::Yellow());
            G4VisAttributes * Pixels = new G4VisAttributes(nexus::DarkRedAlpha());
            Camfront->SetForceSolid(true);
            CamBody->SetForceSolid(true);
            EffectiveArea->SetForceSolid(true);
            Pixels->SetForceSolid(true);

            CamBody_logic->SetVisAttributes(CamBody);

            EffectiveArea_logic->SetVisAttributes(EffectiveArea);
            Pixel_Logic->SetVisAttributes(Pixels);
            CamFront_logic->SetVisAttributes(Camfront);
        }

        this->SetLogicalVolume(lab_logic_volume);




    }

    G4ThreeVector CRAB_EMCCD::GenerateVertex(const G4String &region) const {
        return G4ThreeVector (0,0,0);
    }
}