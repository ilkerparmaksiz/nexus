//
// Created by ilker on 1/14/23.
//

#include "CRABImageIntensifier.h"
#include "FactoryBase.h"

#include <G4Box.hh>
#include <G4Tubs.hh>
#include <MaterialsList.h>
#include "G4GenericMessenger.hh"
#include "G4Tubs.hh"
#include "MaterialsList.h"
#include "SensorSD.h"
#include "OpticalMaterialProperties.h"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"
#include "Visibilities.h"
#include "G4VPhysicalVolume.hh"
#include <G4PVPlacement.hh>

namespace nexus {
    using namespace CLHEP;

    REGISTER_CLASS(CRABImageIntensifier, GeometryBase)

    CRABImageIntensifier::CRABImageIntensifier() : GeometryBase(), msg_(nullptr), fvisual(true),PhotonGain(100),Lab_Logical(nullptr),Lab_Physical(nullptr),Offsetz_(0){
        msg_=new G4GenericMessenger(this,"/Geometry/CRABImageIntesifier/","Control command sof geometry of CRAB TPC");
        msg_->DeclareProperty("gain",PhotonGain,"This is the amount of the photons that are produced");
    }

    CRABImageIntensifier::~CRABImageIntensifier(){
    }
    void CRABImageIntensifier::Construct() {
        double II_z=(25.4/2)*mm;
        double II_r=(25.4/2)*mm;
        double windowthickness=(3/2)*mm;
        G4Tubs * front_solid = new G4Tubs ("II_Front_Window",0,II_r,windowthickness,0,twopi);
        G4Tubs * II_body_solid=new G4Tubs ("II_body",0,II_r,II_z,0,twopi);
        G4Tubs * II_PhotoCathode_solid=new G4Tubs ("II_body",0,II_r,windowthickness,0,twopi);
        G4Tubs * back_solid = new G4Tubs ("II_Back_Window",0,II_r,windowthickness,0,twopi);

        // Change this later on
        G4Material * sapphire_mat =materials::Sapphire();
        G4Material * aluminum=G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");

        if(Lab_Logical== nullptr){
            G4String lab_name="LAB";
            G4Box * lab_solid_volume = new G4Box(lab_name,0.4/2*m,0.4/2*m,0.4/2*m);
            Lab_Logical = new G4LogicalVolume(lab_solid_volume,G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"),lab_name) ;
        }

        if(Lab_Physical== nullptr) Lab_Physical= new G4PVPlacement(0,G4ThreeVector(),Lab_Logical,Lab_Logical->GetName(),0,false,0, false);
        G4LogicalVolume * front_logical= new G4LogicalVolume (front_solid,sapphire_mat,"front_logical");
        G4LogicalVolume * back_logical= new G4LogicalVolume (back_solid,sapphire_mat,"front_logical");
        G4LogicalVolume * II_body_logical= new G4LogicalVolume (II_body_solid,aluminum,"II_body_logical");
        G4LogicalVolume * II_PhotoCathode= new G4LogicalVolume (II_PhotoCathode_solid,aluminum,"II_PhotoCathode_logical");


        // Setting the logical volume
        //this->SetLogicalVolume(II_body_logical);

        // Body Placement
        new G4PVPlacement(0,G4ThreeVector(0,0,Offsetz_),"IIbody",II_body_logical,Lab_Physical, false,0,true);
        new G4PVPlacement(0,G4ThreeVector(0,0,Offsetz_-II_z-windowthickness/2),"IICathode",II_PhotoCathode,Lab_Physical, false,0,true);
        new G4PVPlacement(0,G4ThreeVector(0,0,Offsetz_-II_z-windowthickness),"IIWindowFront",front_logical,Lab_Physical, false,0,true);
        new G4PVPlacement(0,G4ThreeVector(0,0,Offsetz_+II_z+windowthickness/2),"IIWindowBack",back_logical,Lab_Physical, false,0,true);



        // Sensitive Detector
        SensorSD * IIdet= new SensorSD("II_PhotoCathode");
        IIdet->SetDetectorVolumeDepth(2);
        IIdet->SetTimeBinning(100*nanosecond);
        G4SDManager::GetSDMpointer()->AddNewDetector(IIdet);
        II_PhotoCathode->SetSensitiveDetector(IIdet);

        // Will look in to producing optical

        // Visuals
        if (fvisual){
            G4VisAttributes * front = new G4VisAttributes(nexus::YellowAlpha());
            G4VisAttributes * back = new G4VisAttributes(nexus::LightGreenAlpha());
            G4VisAttributes * body = new G4VisAttributes(nexus::DarkGreyAlpha());
            G4VisAttributes * IIcath = new G4VisAttributes(nexus::BlueAlpha());
            front->SetForceSolid(true);
            body->SetForceSolid(true);
            back->SetForceSolid(true);
            IIcath->SetForceSolid(true);

             front_logical->SetVisAttributes(front);
             back_logical->SetVisAttributes(back);
             II_body_logical->SetVisAttributes(body);
             II_PhotoCathode->SetVisAttributes(IIcath);
        }

        this->SetLogicalVolume(Lab_Logical);

    }

    G4ThreeVector CRABImageIntensifier::GenerateVertex(const G4String &region) const {
        return G4ThreeVector (0,0,0);
    }

} // end of the name space