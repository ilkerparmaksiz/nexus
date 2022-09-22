//
// Created by ilker on 9/2/21.
//

#include "CRAB.h"
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

#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Units/PhysicalConstants.h>
#include <G4UnionSolid.hh>
#include <G4UserLimits.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4LogicalBorderSurface.hh>



namespace nexus{
    using namespace CLHEP;
    REGISTER_CLASS(CRAB, GeometryBase)


             CRAB::CRAB():
             GeometryBase(),
             msg_(nullptr),
             Lab_size(1. *m),
             chamber_diam   (15. * cm),
             chamber_length (25. * cm),
             chamber_thickn (2. * mm),
             SourceEn_offset (5.2 *cm),
             SourceEn_diam   (1. * cm),
             SourceEn_length (1 * cm),
             SourceEn_thickn (2. * mm),
             SourceEn_holedia (2. * mm),
             gas_pressure_(10. * bar),
             vtx_(0,0,0),
            sc_yield_(25510./MeV),
            e_lifetime_(1000. * ms),
             pmt_hole_length_ (18.434 * cm),
             MgF2_window_thickness_ (6. * mm),
             MgF2_window_diam_ (10 * mm),
             wndw_ring_stand_out_ (1.5 * mm), //how much the ring around sapph windows stands out of them
             pedot_coating_thickness_ (200. * nanometer), // copied from NEW
             optical_pad_thickness_ (1. * mm), // copied from NEW
             pmt_base_diam_ (47. * mm),
             pmt_base_thickness_ (5. * mm),
             efield_(true),
             HideSourceHolder_(true),
             max_step_size_(1.*mm),
             ElGap_(7*mm),
             ELyield_(970/cm),
             PMT1_Pos_(2.32*cm),
             PMT3_Pos_(3.52*cm)

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

        G4GenericMessenger::Command&  SourceEn_offset_cmd =msg_->DeclarePropertyWithUnit("SourceEn_offset","cm",SourceEn_offset,"SourceEnDiam");
        SourceEn_offset_cmd.SetParameterName("SourceEnoffset", false);

        G4GenericMessenger::Command&  SourceEn_diam_cmd =msg_->DeclarePropertyWithUnit("SourceEn_diam","cm",SourceEn_diam,"SourceEnDiam");
        SourceEn_diam_cmd.SetParameterName("SourceEndiam", false);

        G4GenericMessenger::Command&  SourceEn_length_cmd =msg_->DeclarePropertyWithUnit("SourceEn_length","cm",SourceEn_length,"SourceEnlength");
        SourceEn_length_cmd.SetParameterName("SourceEnlength", false);

        G4GenericMessenger::Command&  SourceEn_holedi_cmd =msg_->DeclarePropertyWithUnit("SourceEn_holedi","cm",SourceEn_holedia,"SourceEnholedi");
        SourceEn_holedi_cmd.SetParameterName("SourceEnholedi", false);

        G4GenericMessenger::Command&  Active_diam_cmd =msg_->DeclarePropertyWithUnit("Active_diam","cm",Active_diam,"ActiveDiam");
        Active_diam_cmd.SetParameterName("Activediam", false);

        G4GenericMessenger::Command&  Active_length_cmd =msg_->DeclarePropertyWithUnit("Active_length","cm",Active_length,"Activelength");
        Active_length_cmd.SetParameterName("Activelength", false);

        G4GenericMessenger::Command&  eliftime_cmd =msg_->DeclarePropertyWithUnit("ElecLifTime","ms",e_lifetime_,"Electron LifeTime");
        eliftime_cmd.SetParameterName("ElecLifTime", false);
        new G4UnitDefinition("1/MeV","1/MeV", "1/Energy", 1/MeV);

        G4GenericMessenger::Command&  scinYield_cmd =msg_->DeclareProperty("scinYield",sc_yield_,"Scintilation Yield Photons/MeV");
        scinYield_cmd.SetParameterName("scinYield", false);
        scinYield_cmd.SetUnitCategory("1/Energy");

        G4GenericMessenger::Command&  ELGap_cmd =msg_->DeclarePropertyWithUnit("ELGap","mm",ElGap_,"Gap which electroluminesence happens");
        Active_length_cmd.SetParameterName("ELGap", false);

        new G4UnitDefinition("1/cm","1/cm", "1/cm", 1/cm);
        G4GenericMessenger::Command&  ELYield_cmd =msg_->DeclareProperty("ELYield",ELyield_,"EL Yield photons/cm");
        ELYield_cmd.SetParameterName("ELYield", false);
        ELYield_cmd.SetUnitCategory("1/cm");

        G4GenericMessenger::Command&  PMT1_Pos_cmd =msg_->DeclarePropertyWithUnit("PMT1_Pos","cm",PMT1_Pos_,"PMT1 Pos");
        PMT1_Pos_cmd.SetParameterName("PMT1_Pos", false);
        G4GenericMessenger::Command&  PMT3_Pos_cmd =msg_->DeclarePropertyWithUnit("PMT3_Pos","cm",PMT3_Pos_,"PMT3 Pos");
        PMT3_Pos_cmd.SetParameterName("PMT3_Pos", false);

        pmt1_=new PmtR7378A();
        pmt2_=new PmtR7378A();

    }

    CRAB::~CRAB()
    {

        delete msg_;
        delete pmt1_;
        delete pmt2_;
    }

    void CRAB::Construct(){


        //Materials
        G4Material* gxe = materials::GXe(gas_pressure_,68);
        G4Material *MgF2=materials::MgF2();
        G4Material *Steel=materials::Steel();
        G4Material *vacuum=G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");



        // Optical Properties Assigned here
        MgF2->SetMaterialPropertiesTable(opticalprops::MgF2());
        vacuum->SetMaterialPropertiesTable(opticalprops::Vacuum());
        gxe->SetMaterialPropertiesTable(opticalprops::GXe(gas_pressure_, 68,sc_yield_,e_lifetime_));
        //Steel->SetMaterialPropertiesTable(opticalprops::STEEL());
        //Constructing Lab Space
        G4String lab_name="LAB";
        G4Box * lab_solid_volume = new G4Box(lab_name,Lab_size/2,Lab_size/2,Lab_size/2);
        G4LogicalVolume * lab_logic_volume= new G4LogicalVolume(lab_solid_volume,G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"),lab_name) ;


        //lab_logic_volume->SetVisAttributes(G4VisAttributes::Invisible);

        //Creating the Steel Cylinder that we need

        /// First Creating the Ends of the Cylinder with Proper Holes
        G4Tubs* chamber_flange_solid =new G4Tubs("CHAMBER_FLANGE", MgF2_window_diam_/2, (chamber_diam/2. + chamber_thickn),( chamber_thickn), 0.,twopi);
        G4LogicalVolume* chamber_flange_logic =new G4LogicalVolume(chamber_flange_solid,materials::Steel(), "CHAMBER_FLANGE");

        // Now Creating The Chamber with Without Ends
        G4Tubs* chamber_solid =new G4Tubs("CHAMBER", chamber_diam/2., (chamber_diam/2. + chamber_thickn),(chamber_length/2. + chamber_thickn), 0.,twopi);
        G4LogicalVolume* chamber_logic =new G4LogicalVolume(chamber_solid,materials::Steel(), "CHAMBER"); //



        /// Sapphire window ///
        G4Tubs* MgF2_window_solid = new G4Tubs("MgF2_WINDOW", 0., MgF2_window_diam_/2.,
                                                   (MgF2_window_thickness_ )/2., 0., twopi);
        G4LogicalVolume* MgF2_window_logic= new G4LogicalVolume(MgF2_window_solid, MgF2, "MgF2_WINDOW");


        //////////////////////////////////////////

        G4double FielCageGap=(160.3+29.55)*mm;
        // Placing the gas in the chamber
        G4Tubs* gas_solid =new G4Tubs("GAS", 0., chamber_diam/2., chamber_length/2., 0., twopi);
        G4LogicalVolume* gas_logic = new G4LogicalVolume(gas_solid, gxe, "GAS");

        //Defining the Detection Area
        //G4Tubs* Active_solid =new G4Tubs("ACTIVE", 0., Active_diam/2., Active_length/2., 0., twopi);
        //G4LogicalVolume* Active_logic = new G4LogicalVolume(Active_solid, gxe, "ACTIVE");


        // EL Region
        G4Tubs* EL_solid =new G4Tubs("EL_GAP", 0., Active_diam/2.,ElGap_/2 , 0., twopi);
        G4LogicalVolume* EL_logic = new G4LogicalVolume(EL_solid, gxe, "EL_GAP");
        //EL_logic->SetUserLimits(new G4UserLimits(1*mm));

        //FieldCage
        G4Tubs* FieldCage_Solid =new G4Tubs("FIELDCAGE", 0., Active_diam/2.,FielCageGap/2 , 0., twopi);
        G4LogicalVolume* FieldCage_Logic = new G4LogicalVolume(FieldCage_Solid, gxe, "FIELDCAGE");
        //FieldCage_Logic->SetUserLimits(new G4UserLimits(1*mm));


        // Radioactive Source Encloser
        //Source
        G4Tubs* SourceHolChamber_solid =new G4Tubs("SourceHolChamber", SourceEn_holedia/2, (SourceEn_diam/2. + SourceEn_thickn),(SourceEn_length/2. + SourceEn_thickn),0,twopi);
        G4LogicalVolume* SourceHolChamber_logic = new G4LogicalVolume(SourceHolChamber_solid,materials::Steel(), "SourceHolChamber_logic");

        G4Tubs* SourceHolChamberBlock_solid =new G4Tubs("SourceHolChBlock",0,(SourceEn_holedia/2),( SourceEn_thickn/2), 0.,twopi);
        G4LogicalVolume* SourceHolChamberBlock_logic = new G4LogicalVolume(SourceHolChamberBlock_solid,materials::Steel(), "SourceHolChBlock_logic");


        //G4LogicalVolume* SourceHolderGas_logic = new G4LogicalVolume(SourceHolChamber_solid, gxe, "SourceHolderGAS_logic");


        //Adding the PMTs in here
        pmt1_->SetPMTName("S1");
        pmt2_->SetPMTName("S2");
        pmt1_->Construct();
        pmt2_->Construct();


        // Adding Logical Volumes for PMTs
        G4LogicalVolume * pmt1_logic=pmt1_->GetLogicalVolume();
        G4LogicalVolume * pmt2_logic=pmt2_->GetLogicalVolume();

        // Particle Source Holder
        //Rotation Matrix

        G4RotationMatrix* rotateHolder = new G4RotationMatrix();
        rotateHolder->rotateY(90.*deg);

        if(!HideSourceHolder_){
            new G4PVPlacement(rotateHolder, G4ThreeVector(-SourceEn_offset,0,0), SourceHolChamber_logic, SourceHolChamber_solid->GetName(),gas_logic, false, 0, false);
            new G4PVPlacement(rotateHolder, G4ThreeVector(-SourceEn_offset-SourceEn_length/2,0,0), SourceHolChamberBlock_logic, SourceHolChamberBlock_solid->GetName(),gas_logic, false, 0, false);
        }


        // PMT1 and PMT3

        // PMTs
        G4double PMT_offset=0.2*cm;
        G4double PMT_pos=(chamber_length/2)+chamber_thickn+(pmt1_->Length()/2)+MgF2_window_thickness_+PMT_offset;
        G4RotationMatrix* pmt1rotate = new G4RotationMatrix();
        pmt1rotate->rotateY(180.*deg);

        //// PMT Covering Tube ///
        G4double offset=1.65*cm;
        G4double PMT_Tube_Length1=MgF2_window_thickness_+(pmt1_->Length()+0.5*cm)/2 + offset -PMT_offset-0.05*cm ;
        G4double PMT_Tube_Length0=(17*cm+pmt1_->Length())/2 - 3.87*cm-PMT_offset;
        G4double PMT_Tube_Block_Thickness=0.2*cm;
        G4double LongPMTTubeOffset=7.5*cm-3.9*cm;
        G4double PMTTubeDiam=2.54*cm;

        // Tube Away from EL
        G4Tubs * PMT_Tube_solid0=new G4Tubs("PMT_TUBE0",(PMTTubeDiam/2)+0.5*cm,(PMTTubeDiam/2)+0.7*cm,PMT_Tube_Length0,0,twopi);
        G4LogicalVolume * PMT_Tube_Logic0=new G4LogicalVolume(PMT_Tube_solid0,materials::Steel(),PMT_Tube_solid0->GetName());

        G4Tubs * PMT_Block_solid0=new G4Tubs("PMT_TUBE_BLOCK0",0,(PMTTubeDiam/2+0.5*cm),PMT_Tube_Block_Thickness,0,twopi);
        G4LogicalVolume * PMT_Block_Logic0=new G4LogicalVolume(PMT_Block_solid0,materials::Steel(),PMT_Block_solid0->GetName());
        // Vacuum for PMT TUBE0
        /// add 2.54*cm to config file
        G4Tubs * InsideThePMT_Tube_solid0=new G4Tubs("PMT_TUBE_VACUUM0",0,(PMTTubeDiam/2+0.4*cm),PMT_Tube_Length0,0,twopi);
        G4LogicalVolume * InsideThePMT_Tube_Logic0=new G4LogicalVolume(InsideThePMT_Tube_solid0,vacuum,InsideThePMT_Tube_solid0->GetName());

        // Tube Close to EL
        G4Tubs * PMT_Tube_solid1=new G4Tubs("PMT_TUBE1",(PMTTubeDiam/2)+0.5*cm,(PMTTubeDiam/2)+0.7*cm,PMT_Tube_Length1,0,twopi);
        G4LogicalVolume * PMT_Tube_Logic1=new G4LogicalVolume(PMT_Tube_solid1,materials::Steel(),PMT_Tube_solid1->GetName());
        G4Tubs * PMT_Block_solid1=new G4Tubs("PMT_TUBE_BLOCK1",0,(PMTTubeDiam/2+0.5*cm),PMT_Tube_Block_Thickness,0,twopi);
        G4LogicalVolume * PMT_Block_Logic=new G4LogicalVolume(PMT_Block_solid1,materials::Steel(),PMT_Block_solid1->GetName());

        // Vacuum for PMT TUBE1

        G4Tubs * InsideThePMT_Tube_solid1=new G4Tubs("PMT_TUBE_VACUUM1",0,(PMTTubeDiam/2+0.4*cm),PMT_Tube_Length1,0,twopi);
        G4LogicalVolume * InsideThePMT_Tube_Logic1=new G4LogicalVolume(InsideThePMT_Tube_solid1,vacuum,InsideThePMT_Tube_solid1->GetName());


        // Place the Volumes

        //LAB

        new G4PVPlacement(0,G4ThreeVector(),lab_logic_volume,lab_logic_volume->GetName(),0,false,0, false);

        //Flanges on the Chamber
        G4VPhysicalVolume *Left_Flange_phys = new G4PVPlacement(0,G4ThreeVector(0,0,(chamber_length/2)),chamber_flange_logic,chamber_flange_solid->GetName(),lab_logic_volume,true,0,false);
        G4VPhysicalVolume *Right_Flange_phys = new G4PVPlacement(0,G4ThreeVector(0,0,-(chamber_length/2)),chamber_flange_logic,chamber_flange_solid->GetName(),lab_logic_volume,true,1,false);

        //Chamber
        G4VPhysicalVolume * chamber_phys=  new G4PVPlacement(0,G4ThreeVector(0.,0.,0) ,chamber_logic, chamber_solid->GetName(), lab_logic_volume, false, 0,false);

        // Xenon Gas in Active Area and Non-Active Area
        G4VPhysicalVolume * gas_phys= new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), gas_logic, gas_solid->GetName(),lab_logic_volume, false, 0, false);
        //new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), Active_logic, Active_solid->GetName(),gas_logic, false, 0, false);


        // FieldCage
        G4double FieldCagePos=chamber_length/2-((129)*mm)-FielCageGap/2-ElGap_/2;
        G4double EL_pos=chamber_length/2-FielCageGap/2-((326)*mm)-ElGap_/2;

        new G4PVPlacement(0,G4ThreeVector(0,0,FieldCagePos/2),FieldCage_Logic,FieldCage_Logic->GetName(),gas_logic, 0,0,false);

        //EL_Gap
        new G4PVPlacement(0,G4ThreeVector(0.,0.,EL_pos/2),EL_logic,EL_solid->GetName(),gas_logic, 0,0, false);


        //  MgF2 Windows
        //G4double window_posz = chamber_length/2+chamber_thickn/2;
        G4double window_posz = chamber_length/2;
        new G4PVPlacement(0, G4ThreeVector(0., 0., window_posz), MgF2_window_logic,"MgF2_WINDOW1", lab_logic_volume,false, 0, false);
        new G4PVPlacement(pmt1rotate, G4ThreeVector(0., 0., -window_posz), MgF2_window_logic,"MgF2_WINDOW2", lab_logic_volume, false, 1, false);

        //PMT Tubes
        G4VPhysicalVolume *PMT_Tube_Phys0=new G4PVPlacement(0,G4ThreeVector(0,0,PMT_pos+LongPMTTubeOffset),PMT_Tube_Logic0,PMT_Tube_Logic0->GetName(),lab_logic_volume,false,0,false);
        G4VPhysicalVolume *PMT_Tube_Phys1=new G4PVPlacement(pmt1rotate,G4ThreeVector(0,0,-(PMT_pos-PMT_offset)-offset),PMT_Tube_Logic1,PMT_Tube_Logic1->GetName(),lab_logic_volume,false,0,false);
        //PMT Tube Vacuum

        G4VPhysicalVolume *PMT_Tube_Vacuum_Phys0=new G4PVPlacement(0,G4ThreeVector(0,0,PMT_pos+LongPMTTubeOffset),InsideThePMT_Tube_Logic0,"PMT_TUBE_VACUUM0",lab_logic_volume,false,0,false);
        G4VPhysicalVolume *PMT_Tube_Vacuum_Phys1=new G4PVPlacement(pmt1rotate,G4ThreeVector(0,0,-(PMT_pos-PMT_offset)-offset),InsideThePMT_Tube_Logic1,"PMT_TUBE_VACUUM1",lab_logic_volume,false,0,false);

        // PMT Tube Block
        new G4PVPlacement(0,G4ThreeVector(0,0,PMT_pos-PMT_offset+PMT_Tube_Length0-PMT_Tube_Block_Thickness/2+LongPMTTubeOffset),PMT_Block_Logic0,PMT_Block_Logic->GetName(),lab_logic_volume,false,0,false);
        new G4PVPlacement(pmt1rotate,G4ThreeVector(0,0,-(PMT_pos-PMT_offset+PMT_Tube_Length1-PMT_Tube_Block_Thickness/2)-offset),PMT_Block_Logic,PMT_Block_Logic->GetName(),lab_logic_volume,false,1,false);

        // PMTs
        new G4PVPlacement(pmt1rotate,G4ThreeVector (0,0,((PMT3_Pos_)-pmt1_->Length()/2-PMT_Tube_Length1/2-MgF2_window_thickness_/2)),pmt1_logic,pmt1_->GetPMTName(),InsideThePMT_Tube_Logic0,true,0,false);
        new G4PVPlacement(pmt1rotate,G4ThreeVector (0,0,(PMT1_Pos_-pmt1_->Length()/2-MgF2_window_thickness_/2)),pmt2_logic,pmt2_->GetPMTName(),InsideThePMT_Tube_Logic1,true,0,false);



        // Define this volume as an ionization sensitive detector
        FieldCage_Logic->SetUserLimits(new G4UserLimits(1*mm));
        IonizationSD* sensdet = new IonizationSD("/CRAB/FIELDCAGE");
        //Active_logic->SetSensitiveDetector(sensdet);
        FieldCage_Logic->SetSensitiveDetector(sensdet);
        G4SDManager::GetSDMpointer()->AddNewDetector(sensdet);


        // Electrical Field
        if(efield_){

            UniformElectricDriftField* field = new UniformElectricDriftField();

            field->SetCathodePosition(FieldCagePos/2+FielCageGap/2);
            field->SetAnodePosition(EL_pos/2+ElGap_/2);
            //field->SetAnodePosition(EL_pos/2);
            field->SetDriftVelocity(.90*mm/microsecond);
            field->SetTransverseDiffusion(.92*mm/sqrt(cm));
            field->SetLongitudinalDiffusion(.36*mm/sqrt(cm));
            G4Region* drift_region = new G4Region("DRIFT");

            drift_region->SetUserInformation(field);
            drift_region->AddRootLogicalVolume(FieldCage_Logic);
            // For CRAB Assuming we have 10 bar gas and Efield is 19,298.20 V/cm
            /// THIS NEEDS TO BE CHANGED

            UniformElectricDriftField * EfieldForEL=new UniformElectricDriftField();
            //EfieldForEL->SetCathodePosition(EL_pos/2+EL_Gap/2);
            EfieldForEL->SetCathodePosition(EL_pos/2+ElGap_/2);
            EfieldForEL->SetAnodePosition(EL_pos/2-ElGap_/2);
            EfieldForEL->SetDriftVelocity(4.61*mm/microsecond);
            EfieldForEL->SetTransverseDiffusion(0.24*mm/sqrt(cm));
            EfieldForEL->SetLongitudinalDiffusion(0.17*mm/sqrt(cm));
            // ELRegion->SetLightYield(xgp.ELLightYield(24.8571*kilovolt/cm));//value for E that gives Y=1160 photons per ie- in normal conditions
            //EfieldForEL->SetLightYield(XenonELLightYield(20*kilovolt/cm, gas_pressure_));
            EfieldForEL->SetLightYield(ELyield_);
            G4Region* el_region = new G4Region("EL_REGION");
            el_region->SetUserInformation(EfieldForEL);
            el_region->AddRootLogicalVolume(EL_logic);


        }
        AssignVisuals();

        /// OpticalSurface
        G4OpticalSurface * OpSteelSurf=new G4OpticalSurface("SteelSurface",unified,polished,dielectric_metal);
        OpSteelSurf->SetMaterialPropertiesTable(opticalprops::STEEL());
        new G4LogicalBorderSurface("SteelSurface_Chamber",gas_phys,chamber_phys,OpSteelSurf);
        new G4LogicalBorderSurface("SteelSurface_LeftFlange",gas_phys,Left_Flange_phys,OpSteelSurf);
        new G4LogicalBorderSurface("SteelSurface_RightFlange",gas_phys,Right_Flange_phys,OpSteelSurf);
        //new G4LogicalBorderSurface("SteelSurface_PMT3_Enclosing",PMT_Tube_Vacuum_Phys0,PMT_Tube_Phys0,OpSteelSurf);
        //new G4LogicalBorderSurface("SteelSurface_PMT1_Enclosing",PMT_Tube_Vacuum_Phys1,PMT_Tube_Phys1,OpSteelSurf);

        //G4Optic
        this->SetLogicalVolume(lab_logic_volume);
        //this->SetLogicalVolume(chamber_logic);

    }


    void CRAB::AssignVisuals() {
        // Chamber
        G4LogicalVolumeStore* lvStore = G4LogicalVolumeStore::GetInstance();




        // Lab
        G4LogicalVolume* Lab = lvStore->GetVolume("LAB");
        G4VisAttributes *LabVa=new G4VisAttributes(G4Colour(2,2,2));
        LabVa->SetForceWireframe(false);

        //Chamber
        G4LogicalVolume* Chamber = lvStore->GetVolume("CHAMBER");
        G4VisAttributes *ChamberVa=new G4VisAttributes(G4Colour(1,1,1));
        ChamberVa->SetForceSolid(true);
        Chamber->SetVisAttributes(G4VisAttributes::GetInvisible());



        //GAS
        G4LogicalVolume* Gas = lvStore->GetVolume("GAS");
        G4VisAttributes *GasVa=new G4VisAttributes(G4Colour(2,2,2));
        GasVa->SetForceCloud(true);
        Gas->SetVisAttributes(GasVa);

        //Source Enclosure Related
        G4LogicalVolume* SourceHolder = lvStore->GetVolume("SourceHolChamber_logic");
        G4LogicalVolume* SourceHolderBlock = lvStore->GetVolume("SourceHolChBlock_logic");
        G4VisAttributes *SourceHolderVa=new G4VisAttributes(G4Colour(2,2,2));
        SourceHolderVa->SetForceSolid(true);

        // Flange
        G4LogicalVolume* flangeLog = lvStore->GetVolume("CHAMBER_FLANGE");
        G4VisAttributes flangeVis=nexus::DarkGreyAlpha();
        flangeVis.SetForceSolid(true);
        flangeLog->SetVisAttributes(ChamberVa);


        //PMT TUBE AND PMT BLOCK
        G4LogicalVolume * PmttubeLog0=lvStore->GetVolume("PMT_TUBE0");
        PmttubeLog0->SetVisAttributes(G4VisAttributes::GetInvisible());
        G4LogicalVolume * PmttubeBlockLog0=lvStore->GetVolume("PMT_TUBE_BLOCK0");
        G4LogicalVolume * PmttubeLog1=lvStore->GetVolume("PMT_TUBE1");
        PmttubeLog1->SetVisAttributes(G4VisAttributes::GetInvisible());
        G4LogicalVolume * PmttubeBlockLog1=lvStore->GetVolume("PMT_TUBE_BLOCK1");
        PmttubeBlockLog0->SetVisAttributes(ChamberVa);
        PmttubeBlockLog1->SetVisAttributes(ChamberVa);
        G4LogicalVolume * PmttubeVacuumLog1=lvStore->GetVolume("PMT_TUBE_VACUUM0");
        G4LogicalVolume * PmttubeVacuumLog2=lvStore->GetVolume("PMT_TUBE_VACUUM1");
        G4VisAttributes PmttubeVacuumVis=nexus::Yellow();
        PmttubeVacuumVis.SetForceCloud(true);
        PmttubeVacuumLog1->SetVisAttributes(PmttubeVacuumVis);
        PmttubeVacuumLog2->SetVisAttributes(PmttubeVacuumVis);


        //MgF2Window
        G4LogicalVolume* MgF2WindowLog = lvStore->GetVolume("MgF2_WINDOW");
        G4VisAttributes  MgF2WindowVis=nexus::LillaAlpha();
        MgF2WindowVis.SetForceSolid(true);
        MgF2WindowLog->SetVisAttributes(MgF2WindowVis);


        // EL-Region
        G4LogicalVolume * ELLogic=lvStore->GetVolume("EL_GAP");
        G4VisAttributes ELVis=nexus::BlueAlpha();
        ELVis.SetForceCloud(true);
        ELLogic->SetVisAttributes(ELVis);

        // FieldCage
        G4LogicalVolume * FieldCage=lvStore->GetVolume("FIELDCAGE");
        G4VisAttributes FielCageVis=nexus::Red();
        FielCageVis.SetForceCloud(true);
        FieldCage->SetVisAttributes(FielCageVis);




        SourceHolder->SetVisAttributes(SourceHolderVa);
        SourceHolderBlock->SetVisAttributes(SourceHolderVa);
        Lab->SetVisAttributes(G4VisAttributes::GetInvisible());

    }


    void CRAB::PrintParam() {

    }

    G4ThreeVector CRAB::GenerateVertex(const G4String& region) const
    {

        if(!(region=="LAB" || region=="GAS" || region=="ACTIVE" || region=="SourceHolChamber" || region=="FIELDCAGE")){

            G4Exception("[CRAB]", "GenerateVertex()", FatalException,
                        "Unknown vertex generation region.");
        }
        return vtx_;
    }


}
