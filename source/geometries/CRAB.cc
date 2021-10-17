//
// Created by ilker on 9/2/21.
//

#include "GRAB.h"

#include "MaterialsList.h"
#include "OpticalMaterialProperties.h"
#include "UniformElectricDriftField.h"
#include "IonizationSD.h"
#include "FactoryBase.h"
#include <G4Box.hh>
#include <G4SubtractionSolid.hh>
#include <G4GenericMessenger.hh>
#include <G4Tubs.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4SDManager.hh>
#include <G4NistManager.hh>
#include <G4VisAttributes.hh>
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"

#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Units/PhysicalConstants.h>


namespace nexus{
    using namespace CLHEP;
    REGISTER_CLASS(CRAB, GeometryBase)


             CRAB::CRAB():
             GeometryBase(),
             msg_(nullptr),
             Lab_size(1. *m),
             chamber_diam   (15. * cm),
             chamber_length (25. * cm),
             chamber_thickn (1. * mm),
             SourceEn_diam   (1. * cm),
             SourceEn_length (3. * cm),
             SourceEn_thickn (2. * mm),
             SourceEn_holedia (2. * mm),
             gas_pressure_(10. * bar),
             vtx_(0,0,0)


    {
        msg_ = new G4GenericMessenger(this, "/Geometry/CRAB/","Control commands of geometry of CRAB TPC");
        G4GenericMessenger::Command&  Pressure_cmd =msg_->DeclarePropertyWithUnit("gas_pressure","bar",gas_pressure_,"Pressure of Gas");
        Pressure_cmd.SetParameterName("XeNonPressure", false);


        G4GenericMessenger::Command&  chamber_diam_cmd =msg_->DeclarePropertyWithUnit("chamber_diam","cm",chamber_diam,"ChamberDiam");
        chamber_diam_cmd.SetParameterName("chamberdiam", false);

        G4GenericMessenger::Command&  chamber_length_cmd =msg_->DeclarePropertyWithUnit("chamber_length","cm",chamber_length,"Chamberlength");
        chamber_length_cmd.SetParameterName("chamberlength", false);

        G4GenericMessenger::Command&  chamber_thickn_cmd =msg_->DeclarePropertyWithUnit("chamber_thickn","mm",chamber_thickn,"Chamberthickn");
        chamber_thickn_cmd .SetParameterName("chamberthickn", false);


        G4GenericMessenger::Command&  source_position_cmd =msg_->DeclarePropertyWithUnit("SourcePosition","cm",vtx_,"vtx");
        source_position_cmd.SetParameterName("vtx", false);

        G4GenericMessenger::Command&  SourceEn_diam_cmd =msg_->DeclarePropertyWithUnit("SourceEn_diam","cm",SourceEn_diam,"SourceEnDiam");
        chamber_diam_cmd.SetParameterName("SourceEndiam", false);

        G4GenericMessenger::Command&  SourceEn_length_cmd =msg_->DeclarePropertyWithUnit("SourceEn_length","cm",SourceEn_length,"SourceEnlength");
        SourceEn_length_cmd.SetParameterName("SourceEnlength", false);

        G4GenericMessenger::Command&  SourceEn_holedi_cmd =msg_->DeclarePropertyWithUnit("SourceEn_holedi","cm",SourceEn_holedia,"SourceEnholedi");
        SourceEn_length_cmd.SetParameterName("SourceEnholedi", false);

    }

    CRAB::~CRAB()
    {
        //delete msg_;
    }

    void CRAB::Construct(){


        //Constructing Lab Space
        G4String lab_name="LAB";
        G4Box * lab_solid_volume = new G4Box(lab_name,Lab_size/2,Lab_size/2,Lab_size/2);
        G4LogicalVolume * lab_logic_volume= new G4LogicalVolume(lab_solid_volume,G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"),lab_name) ;

        //lab_logic_volume->SetVisAttributes(G4VisAttributes::Invisible);


        //Creating the Steel Cylinder that we use
        G4Tubs* chamber_solid =new G4Tubs("CHAMBER", 0., (chamber_diam/2. + chamber_thickn),(chamber_length/2. + chamber_thickn), 0.,twopi);
        G4LogicalVolume* chamber_logic =new G4LogicalVolume(chamber_solid,materials::Steel(), "CHAMBER"); //
        // Placing the gas in the chamber

        G4Tubs* gas_solid =new G4Tubs("GAS", 0., chamber_diam/2., chamber_length/2., 0., twopi);
        // Radioactive Source Encloser
        //Source
        G4Tubs* SourceHolChamber_solid =new G4Tubs("SourceHolChamber", SourceEn_diam/2, (SourceEn_diam/2. + SourceEn_thickn),(SourceEn_length/2. + SourceEn_thickn), 0.,twopi);
        G4Tubs* SourceHolChamberBlock_solid =new G4Tubs("SourceHolChBlock",0,(SourceEn_diam/2 + SourceEn_thickn),( SourceEn_thickn), 0.,twopi);


        //G4VSolid *SourceHolderGas_solid= new G4SubtractionSolid("SourceHolderGas",Source_Chm_solid,SourceHolder_solid);


        G4Material* gxe = materials::GXe(gas_pressure_);
        gxe->SetMaterialPropertiesTable(opticalprops::GXe(gas_pressure_, 68));
        G4LogicalVolume* gas_logic = new G4LogicalVolume(gas_solid, gxe, "GAS");

        //G4LogicalVolume* SourceHolderGas_logic = new G4LogicalVolume(SourceHolderGas_solid, gxe, "SourceHolderGAS_logic");

        G4LogicalVolume* SourceHolChamber_logic = new G4LogicalVolume(SourceHolChamber_solid,materials::Steel(), "SourceHolChamber_logic");
        G4LogicalVolume* SourceHolChamberBlock_logic = new G4LogicalVolume(SourceHolChamberBlock_solid,materials::Steel(), "SourceHolChBlock_logic");


        // Place the Volumes
        new G4PVPlacement(0,G4ThreeVector(),lab_logic_volume,lab_logic_volume->GetName(),0,false,0, true);
        new G4PVPlacement(0,G4ThreeVector(0.,0.,0.) ,chamber_logic, chamber_solid->GetName(), lab_logic_volume, false, 0,true);
        new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), gas_logic, gas_solid->GetName(),chamber_logic, false, 0, true);
        new G4PVPlacement(0, G4ThreeVector(-5*cm,0.,0.), SourceHolChamber_logic, SourceHolChamber_solid->GetName(),gas_logic, false, 0, true);
        new G4PVPlacement(0, G4ThreeVector(-5*cm,0,-SourceEn_length - SourceEn_thickn), SourceHolChamberBlock_logic, SourceHolChamberBlock_solid->GetName(),gas_logic, false, 0, true);




        // Define this volume as an ionization sensitive detector
        IonizationSD* sensdet = new IonizationSD("/CRAB/GAS");
        gas_logic->SetSensitiveDetector(sensdet);
        G4SDManager::GetSDMpointer()->AddNewDetector(sensdet);

        AssignVisuals();
        this->SetLogicalVolume(lab_logic_volume);
        //this->SetLogicalVolume(chamber_logic);

    }


    void CRAB::AssignVisuals() {
        // Chamber
        G4LogicalVolumeStore* lvStore = G4LogicalVolumeStore::GetInstance();
        G4LogicalVolume* Chamber = lvStore->GetVolume("CHAMBER");
        G4LogicalVolume* Lab = lvStore->GetVolume("LAB");
        G4LogicalVolume* SourceHolder = lvStore->GetVolume("SourceHolChamber_logic");
        G4LogicalVolume* SourceHolderBlock = lvStore->GetVolume("SourceHolChBlock_logic");
        G4VisAttributes *SourceHolderVa=new G4VisAttributes(G4Colour(2,2,2));
        G4VisAttributes *ChamberVa=new G4VisAttributes(G4Colour(1,1,1));
        G4VisAttributes *LabVa=new G4VisAttributes(G4Colour(2,2,2));

        LabVa->SetForceWireframe(true);
        SourceHolderVa->SetForceWireframe(true);
        ChamberVa->SetForceSolid(true);
        SourceHolder->SetVisAttributes(SourceHolderVa);
        SourceHolderBlock->SetVisAttributes(SourceHolderVa);
        Lab->SetVisAttributes(LabVa);
        Chamber->SetVisAttributes(G4VisAttributes::Invisible);

    }

    void CRAB::PrintParam() {

    }

    G4ThreeVector CRAB::GenerateVertex(const G4String& region) const
    {

        if(!(region=="LAB" || region=="GAS" || region=="ACTIVE" )){

            G4Exception("[CRAB]", "GenerateVertex()", FatalException,
                        "Unknown vertex generation region.");
        }
        return vtx_;
    }


}
