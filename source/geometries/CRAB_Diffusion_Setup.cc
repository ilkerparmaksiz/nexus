//
// Created by ilker on 9/2/21.
//

#include "CRAB_Diffusion_Setup.h"
#include "Visibilities.h"
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
#include "XenonProperties.h"
#include "G4UnitsTable.hh"
#include "HexagonMeshTools.h"

#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Units/PhysicalConstants.h>
#include <G4UnionSolid.hh>
#include <G4UserLimits.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4LogicalBorderSurface.hh>
#include "../physics/UltraFresnelLens.hh"
#include "G4Polyhedra.hh"
#include "G4ExtrudedSolid.hh"
#include "G4Orb.hh"
#include "G4IntersectionSolid.hh"
#include "CRAB_EMCCD.h"
namespace nexus{
    using namespace CLHEP;
    REGISTER_CLASS(CRAB_Diffusion_Setup, GeometryBase)
        CRAB_Diffusion_Setup::CRAB_Diffusion_Setup():
             GeometryBase(),
             msg_(nullptr)
    {
        msg_ = new G4GenericMessenger(this, "/Geometry/CRAB_Diffusion_Setup/","Control commands of geometry of CRAB TPC");
        det=new CRAB_Detector();
        AnodePMT = new PmtR7378A();
        theCamera = new CRAB_EMCCD();
        II = new CRABImageIntensifier();
    }

    CRAB_Diffusion_Setup::~CRAB_Diffusion_Setup()
    {
        delete msg_;
    }

    void CRAB_Diffusion_Setup::Construct() {
        G4String lab_name="LAB";
        G4Box * lab_solid_volume = new G4Box(lab_name,(1*m)/2,(1*m)/2,(1*m)/2);
        G4LogicalVolume * lab_logic_volume= new G4LogicalVolume(lab_solid_volume,G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"),lab_name) ;

        //Materials
        G4Material *MgF2=materials::MgF2();
        G4Material *Steel=materials::Steel();
        G4Material *vacuum=G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");

        //Set This Lab to be Logical Volumes for all the geometries

        // Lab Logical
        det->SetLabLogical(lab_logic_volume);
        II->SetLabLogical(lab_logic_volume);
        theCamera->SetLabLogical(lab_logic_volume);


        //Placements
        II->Offsetz_=30*cm;
        theCamera->Offsetz_=40*cm;


        // call the construct function for all the geometries
        theCamera->Construct(); // EMCCD
        II->Construct(); // II
        det->Construct(); // CRAB_Detector
        AnodePMT->SetPMTName("Anode"); // Set the name of the PMT
        AnodePMT->Construct(); // Anode PMT

        // Lens

        // Image Intensifier

        // CRAB Detector

        //Component Placement


        //G4VPhysicalVolume *Camera_Physical=new G4PVPlacement(0,G4ThreeVector(0,0,40*cm),theCamera->GetLogicalVolume(),theCamera->GetLogicalVolume()->GetName(),lab_logic_volume,false,0,false);
        //G4VPhysicalVolume *II_Physical=new G4PVPlacement(0,G4ThreeVector(0,0,30*cm),II->GetLogicalVolume(),II->GetLogicalVolume()->GetName(),lab_logic_volume,false,0,false);

        //Internal Components

               // Needle 1

               // Needle 2

               // Needle 3


        // Anode PMT

        G4ThreeVector Anode_PMT_Pos;
        Anode_PMT_Pos=G4ThreeVector(0,0,-26.72*cm);
        // Place the PMT to the position
        G4double AnodePMTTubeDiam=2.54*cm;
        G4double offset=1.65*cm;
        G4double AnodePMT_Tube_Block_Thickness=0.2*cm;
        G4double AnodePMT_offset=0.2*cm;
        G4double AnodePMT_Tube_Length=(AnodePMT->Length()+0.5*cm)/2 + offset -AnodePMT_offset-0.05*cm ;

        G4Tubs * AnodePMT_Tube_solid=new G4Tubs("AnodeAnodePMT_TUBE",(AnodePMTTubeDiam/2)+0.5*cm,(AnodePMTTubeDiam/2)+0.7*cm,AnodePMT_Tube_Length,0,twopi);
        G4LogicalVolume * AnodePMT_Tube_Logic=new G4LogicalVolume(AnodePMT_Tube_solid,materials::Steel(),AnodePMT_Tube_solid->GetName());

        G4Tubs * AnodePMT_Block_solid=new G4Tubs("AnodeAnodePMT_TUBE_BLOCK",0,(AnodePMTTubeDiam/2+0.5*cm),AnodePMT_Tube_Block_Thickness,0,twopi);
        G4LogicalVolume * AnodePMT_Block_Logic=new G4LogicalVolume(AnodePMT_Block_solid,materials::Steel(),AnodePMT_Block_solid->GetName());

        // Vacuum Tube For Anode AnodePMT
        G4Tubs * InsideTheAnodePMT_Tube_solid=new G4Tubs("AnodeAnodePMT_TUBE_VACUUM",0,(AnodePMTTubeDiam/2+0.5*cm),AnodePMT_Tube_Length,0,twopi);
        G4LogicalVolume * InsideTheAnodePMT_Tube_Logic=new G4LogicalVolume(InsideTheAnodePMT_Tube_solid,vacuum,InsideTheAnodePMT_Tube_solid->GetName());

        G4Tubs * AnodePMT_Block_solid1=new G4Tubs("AnodeAnodePMT_TUBE_BLOCK",0,(AnodePMTTubeDiam/2+0.5*cm),AnodePMT_Tube_Block_Thickness,0,twopi);

        // Placing PMT Related Geometries
        G4VPhysicalVolume *AnodePMT_Tube_Phys=new G4PVPlacement(0,Anode_PMT_Pos-G4ThreeVector(0,0,(2*AnodePMT_offset)-offset),AnodePMT_Tube_Logic,AnodePMT_Tube_Logic->GetName(),lab_logic_volume,false,0,false);
        G4VPhysicalVolume *AnodePMT_TubeBlock_Phys=new G4PVPlacement(0,Anode_PMT_Pos-G4ThreeVector(0,0,(2*AnodePMT_offset)-offset+AnodePMT_Tube_Length+AnodePMT_Tube_Block_Thickness/2),AnodePMT_Block_Logic,AnodePMT_Tube_Logic->GetName(),lab_logic_volume,false,0,false);
        new G4PVPlacement(0,Anode_PMT_Pos-G4ThreeVector(0,0,(2*AnodePMT_offset)-offset),InsideTheAnodePMT_Tube_Logic,InsideTheAnodePMT_Tube_Logic->GetName(),lab_logic_volume,false,0,false);
        G4VPhysicalVolume *AnodePMT_Phys=new G4PVPlacement(0,G4ThreeVector (0,0,AnodePMT_Tube_Length/2+AnodePMT->Length()/2-5*cm),AnodePMT->GetLogicalVolume(),"Anode_PMT",InsideTheAnodePMT_Tube_Logic,0,false);

        //Visuals for PMT TUBE
        G4VisAttributes PmttubeVacuumVis=nexus::DarkGreyAlpha();
        PmttubeVacuumVis.SetForceCloud(true);
        InsideTheAnodePMT_Tube_Logic->SetVisAttributes(PmttubeVacuumVis);

        // Select the logical volume
        this->SetLogicalVolume(lab_logic_volume);

    }
}
