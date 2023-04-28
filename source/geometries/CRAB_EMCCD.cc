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
    CRAB_EMCCD::CRAB_EMCCD() :GeometryBase(),msg_(nullptr),fvisual(true),Lab_Logical(nullptr),Lab_Physical(nullptr),Offsetz_(0){
        msg_=new G4GenericMessenger(this,"/Geometry/CRAB_EMMCD/","This is the geometry of the EMCCD");
    }

    // Destruct
    CRAB_EMCCD::~CRAB_EMCCD() noexcept {}

    void CRAB_EMCCD::Construct() {
        // C91000-23B Camera Pixel Size
        float pixelSize=16*um;
        //float pixelSize=5*mm;
        float effectiveXandY=8.192*mm;

        // LAB

        if(Lab_Logical== nullptr){
            G4String lab_name = "Lab";
            G4Box * lab_solid_volume = new G4Box(lab_name,0.4/2*m,0.4/2*m,0.4/2*m);
            Lab_Logical= new G4LogicalVolume(lab_solid_volume,G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"),lab_name) ;
        }
        if(Lab_Physical== nullptr)  Lab_Physical= new G4PVPlacement(0,G4ThreeVector(),Lab_Logical,Lab_Logical->GetName(),0,false,0, false);

        // Solids
        G4Box * CamBody_solid=new G4Box("CamBody",(110/2)*mm,(110/2)*mm,(156/2)*mm);

        G4Tubs * CamFront_solid=new G4Tubs("CamFront",(25.4/2)*mm,(76/2)*mm,(25/2)*mm,0,twopi);

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
        new G4PVPlacement(0,G4ThreeVector(0,0,Offsetz_),CamBody_solid->GetName(),CamBody_logic,Lab_Physical,false,0,false);
        new G4PVPlacement(0,G4ThreeVector(0,0,Offsetz_+(-156/2-25/2)*mm),CamFront_logic->GetName(),CamFront_logic,Lab_Physical,false,0,false);
        new G4PVPlacement(0,G4ThreeVector(0,0,Offsetz_+(-156/2-1)*mm),EffectiveArea_logic->GetName(),EffectiveArea_logic,Lab_Physical,false,0,false);
        //new G4PVPlacement(0,G4ThreeVector(0,0,(-156/2-1-pixelSize/(4))*mm),Pixel_Logic->GetName(),Pixel_Logic,labPhysical,true,0,false);

        G4int PixelSize=512;
        G4int XOffset=256;
        G4int YOffset=256;

        /*G4float zposition=(-156/2-1-pixelSize/(4)-25/2)*mm;
        G4int PixelCount = 0;
        for (int px=0;px<PixelSize;px++){
            for (int py=0;py<PixelSize;py++){
                G4String PixelName=Pixel_Logic->GetName()+"_"+std::to_string(PixelCount);
                new G4PVPlacement(0,G4ThreeVector((px-XOffset)*pixelSize,(py-YOffset)*pixelSize,zposition),PixelName,Pixel_Logic,labPhysical,true,PixelCount,false);
                PixelCount++;
            }
        }*/
        // Visuals
        if (fvisual){
            G4VisAttributes * Camfront = new G4VisAttributes(nexus::DarkRed());
            G4VisAttributes * CamBody = new G4VisAttributes(nexus::DarkGrey());
            G4VisAttributes * EffectiveArea = new G4VisAttributes(nexus::Yellow());
            G4VisAttributes * Pixels = new G4VisAttributes(nexus::Red());
            Camfront->SetForceSolid(true);
            CamBody->SetForceSolid(true);
            EffectiveArea->SetForceSolid(true);
            Pixels->SetForceSolid(true);
            CamBody_logic->SetVisAttributes(CamBody);
            EffectiveArea_logic->SetVisAttributes(EffectiveArea);
            Pixel_Logic->SetVisAttributes(Pixels);
            CamFront_logic->SetVisAttributes(Camfront);
        }
        this->SetLogicalVolume(Lab_Logical);



    }

    G4ThreeVector CRAB_EMCCD::GenerateVertex(const G4String &region) const {
        return G4ThreeVector (0,0,0);
    }
}