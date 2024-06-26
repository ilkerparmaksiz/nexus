#include "CRAB0.h"
#include "G4PVParameterised.hh"
#include "G4PVReplica.hh"
#include "G4RotationMatrix.hh"
#include "G4UnitsTable.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4OpticalSurface.hh"
#include "G4Trd.hh"
#include "G4RegionStore.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4Cons.hh"
#include "G4IntersectionSolid.hh"
#include "G4Trd.hh"
#include "config.h"

#ifdef With_GarField
#include "DegradModel.h"
#include "GarfieldVUVPhotonModel.h"
#include "GarfieldHelper.h"
#else
#include "UniformElectricDriftField.h"
#endif
#include "G4SDManager.hh"
#include "G4Polyhedra.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4ExtrudedSolid.hh"
#include "G4MultiUnion.hh"
#include "Visibilities.h"
#include "HexagonMeshTools.h"
#include "FactoryBase.h"
#include "IonizationSD.h"
#include "MaterialsList.h"
#include "SensorSD.h"
#include "G4MultiUnion.hh"
#include "G4SubtractionSolid.hh"
#include "G4BooleanSolid.hh"
#include <G4GenericMessenger.hh>


namespace nexus {
    
    REGISTER_CLASS(CRAB0,GeometryBase)
    
    CRAB0::CRAB0() :
            checkOverlaps(0),
            temperature(300 * kelvin), // temperature
            Lab_size(1 * m),
            chamber_diam(16.4 * cm),
            chamber_length(43.18 * cm), // Config files vary
            chamber_thickn(7 * mm),
            SourceEn_offset(5.7 * cm),
            SourceEn_diam(1.0 * cm),
            SourceEn_length(1 * cm),
            SourceEn_thickn(2. * mm),
            SourceEn_holedia(5. * mm),
            gas_pressure_(10 * bar),
            vtx_(-1.6 * cm , 0, -5 * cm),
            Active_diam(8.6 * cm),
            sc_yield_(25510. / MeV),
            e_lifetime_(3000. * ms),
            MgF2_window_thickness_(6. * mm),
            Anode_window_diam_(16.22*mm),
            Cathode_window_diam_(16.58*mm),
            HideSourceHolder_(false),
            max_step_size_(1. * mm),
            ElGap_(7 * mm),
            ELyield_(925 / cm),
            PMT1_Pos_(2.08 * cm),
            PMT3_Pos_(3.52 * cm),
            HideCollimator_(true),
            specific_vertex_{},
            fOffset(-0.8*cm),
            GasFile_("data/Xenon_10Bar.gas"),
            useCAD_(false)
           {

            // Messenger
            msg_ = new G4GenericMessenger(this, "/Geometry/CRAB0/","Control commands of geometry of CRAB0.");

            G4GenericMessenger::Command& DetChamberR = msg_->DeclareProperty("DetChamberD", chamber_diam, "Set the detector chamber diameter");
            DetChamberR.SetUnitCategory("Length");
            DetChamberR.SetParameterName("DetChamberD", false);
            DetChamberR.SetRange("DetChamberD>0.");

            G4GenericMessenger::Command& DetChamberL = msg_->DeclareProperty("DetChamberL", chamber_length, "Set the detector chamber length");
            DetChamberL.SetUnitCategory("Length");
            DetChamberL.SetParameterName("DetChamberL", false);
            DetChamberL.SetRange("DetChamberL>0.");

            G4GenericMessenger::Command& DetActiveD = msg_->DeclareProperty("DetActiveD", Active_diam, "Set the detector active diameter");
            DetActiveD.SetUnitCategory("Length");
            DetActiveD.SetParameterName("DetActiveD", false);
            DetActiveD.SetRange("DetActiveD>0.");

            G4GenericMessenger::Command& DetActiveL = msg_->DeclareProperty("DetActiveL", FielCageGap, "Set the detector active length");
            DetActiveL.SetUnitCategory("Length");
            DetActiveL.SetParameterName("DetActiveL", false);
            DetActiveL.SetRange("DetActiveL>0.");

            G4GenericMessenger::Command& GasPressure = msg_->DeclarePropertyWithUnit("GasPressure","bar",gas_pressure_, "Set the gas pressure");
            GasPressure.SetUnitCategory("Pressure");
            GasPressure.SetParameterName("GasPressure", false);
            GasPressure.SetRange("GasPressure>0.");

            G4GenericMessenger::Command& fieldDrift = msg_->DeclareProperty("fieldDrift", fieldDrift_, "Set the drift field [V/cm]");
            fieldDrift.SetParameterName("fieldDrift", false);
            fieldDrift.SetRange("fieldDrift>0.");

            G4GenericMessenger::Command& fieldEL = msg_->DeclareProperty("fieldEL", fieldEL_, "Set the EL field [V/cm]");
            fieldEL.SetParameterName("fieldEL", false);
            fieldEL.SetRange("fieldEL>0.");

            msg_->DeclarePropertyWithUnit("specific_vertex", "mm",  specific_vertex_, "Set generation vertex.");

            msg_->DeclareProperty("useCAD", useCAD_, "Use CAD geometry or G4 bools");

            msg_->DeclareProperty("GasFile", GasFile_, "Use CAD geometry or G4 bools");


        Sampler=std::make_shared<SampleFromSurface>(SampleFromSurface("Needles"));

        }

    CRAB0::~CRAB0() {
    }

    void CRAB0::Construct() {

        //  ------------------------ Materials --------------------------------
        gxe = materials::GXe(gas_pressure_, 68);
        MgF2 = materials::MgF2();
        Steel = materials::Steel();
        PEEK = materials::PEEK();
        vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
        teflon = G4NistManager::Instance()->FindOrBuildMaterial("G4_TEFLON");
        G4Material *Air=G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");



        //  ------------------- Optical Properties ----------------------------
        MgF2->SetMaterialPropertiesTable(opticalprops::MgF2());
        Air->SetMaterialPropertiesTable(opticalprops::Vacuum());
        gxe->SetMaterialPropertiesTable(opticalprops::GXe(gas_pressure_, 68, sc_yield_, e_lifetime_));
        Steel->SetMaterialPropertiesTable(opticalprops::STEEL());
        vacuum->SetMaterialPropertiesTable(opticalprops::Vacuum());
        //teflon->SetMaterialPropertiesTable(opticalprops::STEEL());
        //PEEK->SetMaterialPropertiesTable(opticalprops::STEEL());




        //  ----------------------- Lab Space ---------------------------------
        G4String lab_name = "LAB";
        G4Box *lab_solid_volume = new G4Box(lab_name, Lab_size / 2, Lab_size / 2, Lab_size / 2);
        G4LogicalVolume *lab_logic_volume = new G4LogicalVolume(lab_solid_volume,Air, lab_name);



        // Xenon Gas
        G4Tubs *gas_solid = new G4Tubs("GAS", 0., chamber_diam / 2.+ chamber_thickn, chamber_length / 2. + chamber_thickn, 0., twopi);



        // Define a rotation matrix to orient all detector pieces along y direction
        G4RotationMatrix *rotateX = new G4RotationMatrix();
        rotateX->rotateX(90. * deg);

        // Define a rotation matrix for the meshes so they point in the right direction
        G4RotationMatrix *rotateMesh = new G4RotationMatrix();
        rotateMesh->rotateZ(30. * deg);


        //  ------------------------ Vessel -----------------------------------
        // Steel end caps with holes


        G4Tubs* chamber_flange_solid_Anode = new G4Tubs("CHAMBER_FLANGE_ANODE", Anode_window_diam_/2, (chamber_diam/2. + chamber_thickn), chamber_thickn/2.0, 0., twopi);
        G4LogicalVolume* chamber_flange_logic_Anode =new G4LogicalVolume(chamber_flange_solid_Anode,materials::Steel(), "CHAMBER_FLANGE_ANODE");

        G4Tubs* chamber_flange_solid_Cathode = new G4Tubs("CHAMBER_FLANGE_CATHODE", Cathode_window_diam_/2, (chamber_diam/2. + chamber_thickn), chamber_thickn/2.0, 0., twopi);
        G4LogicalVolume* chamber_flange_logic_Cathode =new G4LogicalVolume(chamber_flange_solid_Cathode,materials::Steel(), "CHAMBER_FLANGE_CATHODE");


        // Chamber barrel
        G4Tubs *chamber_solid = new G4Tubs("CHAMBER", chamber_diam / 2., (chamber_diam / 2. + chamber_thickn),
                                           (chamber_length / 2), 0., twopi);
        G4LogicalVolume *chamber_logic = new G4LogicalVolume(chamber_solid,Steel, "CHAMBER");


        gas_logic = new G4LogicalVolume(gas_solid, gxe, "GAS");


        // Xenon Gas in Active Area and Non-Active Area
        G4VPhysicalVolume *gas_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), gas_logic, gas_solid->GetName(),
                                                        lab_logic_volume, false, 0, false);
        // --- Placement ---

        // Flanges on the Chamber, place in the gas logic so we include the aperature region
        G4VPhysicalVolume *Left_Flange_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, chamber_length / 2 +
                                                                                       chamber_thickn / 2.0),
                                                                chamber_flange_logic_Cathode, chamber_flange_solid_Cathode->GetName(),
                                                                gas_logic, false, 0, checkOverlaps);

        G4VPhysicalVolume *Right_Flange_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, -chamber_length / 2 -
                                                                                        chamber_thickn / 2.0),
                                                                 chamber_flange_logic_Anode, chamber_flange_solid_Anode->GetName(),
                                                                 gas_logic, false, 0, checkOverlaps);


        G4VPhysicalVolume *chamber_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., 0), chamber_logic,chamber_solid->GetName(), gas_logic, false, 0,checkOverlaps);


        // OpticalSurface


        // Add optical surface
        G4OpticalSurface* OpSteelSurf = new G4OpticalSurface("OPSURF");
        OpSteelSurf->SetMaterialPropertiesTable(opticalprops::STEEL());
        OpSteelSurf->SetType(dielectric_metal);
        OpSteelSurf->SetModel(unified);
        OpSteelSurf->SetFinish(polished);
        // gas_mesh_opsur->SetSigmaAlpha(0.0);
        // --- Optical ---
        new G4LogicalBorderSurface("SteelSurface_Chamber",gas_phys,chamber_phys, OpSteelSurf);
        new G4LogicalBorderSurface("SteelSurface_Flange_Right", gas_phys,Right_Flange_phys, OpSteelSurf);
        new G4LogicalBorderSurface("SteelSurface_Flange_Left", gas_phys,Left_Flange_phys, OpSteelSurf);

        //  --------------------- Window/lens ---------------------------------
        // MgF2 window
        G4Tubs *MgF2_window_solid = new G4Tubs("MgF2_WINDOW", 0., Anode_window_diam_ / 2.,
                                               (MgF2_window_thickness_) / 2., 0., twopi);
        G4LogicalVolume *MgF2_window_logic = new G4LogicalVolume(MgF2_window_solid, MgF2, "MgF2_WINDOW");
        // Lens
        const G4double lensRcurve(2.83 * cm); // radius of curvature of MgF2 Lens
        const G4ThreeVector posLensTubeIntersect(0., 0., -lensRcurve);

        // Lens is made from the intersection of a sphere and a cylinder
        G4double maxLensLength = 4 * mm;
        G4Tubs *sLensTube = new G4Tubs("sLensSphereTube", 0, Cathode_window_diam_ / 2, maxLensLength, 0.,
                                       twopi); // 4 mm is the max lens length
        G4Orb *sLensOrb = new G4Orb("sLensSphere", lensRcurve);
        G4IntersectionSolid *sLens = new G4IntersectionSolid("sLens", sLensTube, sLensOrb, 0, posLensTubeIntersect);

        // Lens logical
        G4LogicalVolume *lensLogical = new G4LogicalVolume(sLens, MgF2, "Lens");

        // --- Placement ---
        G4double window_posz = chamber_length / 2 + chamber_thickn;
        // G4VPhysicalVolume* lensPhysical = new G4PVPlacement(0, G4ThreeVector(0., 0., window_posz), MgF2_window_logic,"MgF2_WINDOW1", lab_logic_volume,false, 0, checkOverlaps);



        G4VPhysicalVolume *lensPhysical = new G4PVPlacement(0, G4ThreeVector(0., 0., window_posz + maxLensLength / 2.0),lensLogical, "MgF2_LENS_CATHODE", gas_logic, false, 0, checkOverlaps);
        G4VPhysicalVolume *MgF2WindowPhysical= new G4PVPlacement(0, G4ThreeVector(0., 0., -window_posz), MgF2_window_logic, "MgF2_WINDOW_ANODE", gas_logic, false,1, checkOverlaps);


        //  --------------------------- EL ------------------------------------
        G4double PeekRodExtend=0.6*cm;
        FielCageGap = 21.26 * cm-PeekRodExtend;

        // EL Region
        G4Tubs *EL_solid = new G4Tubs("EL_GAP", 0., Active_diam / 2., ElGap_ / 2, 0., twopi);
        G4LogicalVolume *EL_logic = new G4LogicalVolume(EL_solid, gxe, "EL_GAP");

        G4double EL_ID         = 7.2 * cm;
        G4double EL_OD         = 12.0 * cm;
        G4double EL_thick      = 1.3 * cm;
        G4double EL_mesh_thick = 0.1 * mm;
        G4double EL_hex_size   = 2.5 * mm;

        G4Tubs *EL_ring_solid = new G4Tubs("EL_Ring", EL_ID / 2., EL_OD / 2.0, EL_thick / 2, 0., twopi);
        G4LogicalVolume *EL_ring_logic = new G4LogicalVolume(EL_ring_solid, Steel, "EL_Ring");
        G4LogicalVolume *Cathode_ring_logic = new G4LogicalVolume(EL_ring_solid, Steel, "Cathode_Ring");





        // Define a hexagonal prism
///// Off for Opticks /////
        // EL Mesh
        /*
        // Dist from centre of hex to hex vertex, excluding the land width (circumradius)
        G4double hex_circumR = EL_hex_size / std::sqrt(3);

        // Number of hexagons needed -- need to use fixed amount, too many and nexus will crash
        G4int nHole = 16;

        // Define the Stainless steel mesh cylinder to subtract hex pattern from
        G4Tubs *Mesh_Disk = new G4Tubs("Mesh_Disk", 0., EL_OD / 2.0, EL_mesh_thick / 2., 0.,
                                       twopi); // Use OD so mesh stays within the logical

        HexagonMeshTools::HexagonMeshTools *HexCreator; // Hexagonal Mesh Tool
        G4ExtrudedSolid *HexPrism = HexCreator->CreateHexagon(EL_mesh_thick, hex_circumR);
        G4LogicalVolume *ELP_Disk_logic     = new G4LogicalVolume(Mesh_Disk, Steel, "ELP_Mesh_Logic");
        G4LogicalVolume *ELPP_Disk_logic    = new G4LogicalVolume(Mesh_Disk, Steel, "ELPP_Mesh_Logic");
        G4LogicalVolume *Cathode_Disk_logic = new G4LogicalVolume(Mesh_Disk, Steel, "Cathode_Mesh_Logic");
        G4LogicalVolume *EL_Hex_logic       = new G4LogicalVolume(HexPrism,  gxe,   "Mesh_Hex");
        */
///// Off for Opticks End /////


        //  ----------------------- Field Cage --------------------------------

        G4double FieldCagePos = -0 * cm+PeekRodExtend;
        G4double EL_pos = -10.98 * cm+PeekRodExtend;

        // FieldCage -- needs to be updated to rings and PEEK rods
        G4Tubs *FieldCage_Solid = new G4Tubs("ACTIVE", 0., Active_diam / 2., FielCageGap / 2, 0., twopi);
        G4LogicalVolume *FieldCage_Logic = new G4LogicalVolume(FieldCage_Solid, gxe, "ACTIVE");
#ifndef With_GarField
        new G4PVPlacement(0,G4ThreeVector(0,0,FieldCagePos/2),FieldCage_Logic,FieldCage_Logic->GetName(),gas_logic, 0,0,checkOverlaps);
#endif
        // Field Rings
        G4double FR_ID = 8.6 * cm; // Field Ring Inner Diameter
        G4double FR_OD = 11.5 * cm; // Field Ring Outer Diameter
        G4double FR_thick = 3.5 * mm; // Field Ring thickness

        G4Tubs *FR_Solid = new G4Tubs("FR", FR_ID / 2., FR_OD / 2., FR_thick / 2.0, 0., twopi);
        G4LogicalVolume *FR_logic = new G4LogicalVolume(FR_Solid, Steel, "FR");

        // PEEK Rods
        G4double PEEK_Rod_OD = 1.6 * cm;    // PEEK Rods Outer Diameter
        G4double PEEK_Rod_thick = 1.27 * cm; // PEEK Rods thickness

        G4Tubs *PEEK_Rod_Solid = new G4Tubs("PEEK_Rod", 0., PEEK_Rod_OD / 2., PEEK_Rod_thick / 2.0, 0., twopi);
        G4LogicalVolume *PEEK_logic = new G4LogicalVolume(PEEK_Rod_Solid, PEEK, "PEEK_Rod");

        G4Tubs *PEEK_Rod_Solid_cathode = new G4Tubs("PEEK_Rod", 0., PEEK_Rod_OD / 2., 1 * cm / 2.0, 0., twopi);
        G4LogicalVolume *PEEK_logic_cathode = new G4LogicalVolume(PEEK_Rod_Solid_cathode, PEEK, "PEEK_Rod_C");

        G4Tubs *PEEK_Rod_Solid_buffer = new G4Tubs("PEEK_Rod", 0., PEEK_Rod_OD / 2., 1.1 * cm / 2.0, 0., twopi);
        G4LogicalVolume *PEEK_logic_buffer = new G4LogicalVolume(PEEK_Rod_Solid_buffer, PEEK, "PEEK_Rod_B");

        G4Tubs *PEEK_Rod_Solid_buffer_end = new G4Tubs("PEEK_Rod", 0., PEEK_Rod_OD / 2., 3.37 * cm / 2.0+PeekRodExtend/2, 0., twopi);
        G4LogicalVolume *PEEK_logic_buffer_end = new G4LogicalVolume(PEEK_Rod_Solid_buffer_end, PEEK, "PEEK_Rod_BE");



        // --- Placement ---
       std::vector<G4VPhysicalVolume *> FieldRingsPhysical;
        // Field Cage Rings
        for (int i=1; i<5;i++){
            FieldRingsPhysical.push_back( new G4PVPlacement(0, G4ThreeVector(0, 0, -FR_thick / 2.0 +i * (FR_thick + PEEK_Rod_thick)),FR_logic, FR_logic->GetName(), gas_logic, 1, i, checkOverlaps));
            FieldRingsPhysical.push_back( new G4PVPlacement(0, G4ThreeVector(0, 0, -FR_thick / 2.0 -i * (FR_thick + PEEK_Rod_thick)),FR_logic, FR_logic->GetName(), gas_logic, 1, i+5, checkOverlaps));


        }


        FieldRingsPhysical.push_back( new G4PVPlacement(0, G4ThreeVector(0, 0, -FR_thick / 2.0 ),FR_logic, FR_logic->GetName(), gas_logic, 1, 5, checkOverlaps));
        FieldRingsPhysical.push_back( new G4PVPlacement(0, G4ThreeVector(0, 0, -FR_thick / 2.0 +5 * (FR_thick + PEEK_Rod_thick)),FR_logic, FR_logic->GetName(), gas_logic, 1, 11, checkOverlaps));

        // // EL Field Rings
        FieldRingsPhysical.push_back(new G4PVPlacement(0, G4ThreeVector(0, 0, PeekRodExtend-FR_thick / 2.0 -
                                                                              4 * (FR_thick + PEEK_Rod_thick) -
                                                                              2.5 * cm - 1.3 * cm - 0.7 * cm -
                                                                              1.3 * cm - 2 * cm - FR_thick -
                                                                              2 * (FR_thick + PEEK_Rod_thick)),
                                                       FR_logic, FR_logic->GetName(), gas_logic, 1, 12, checkOverlaps));
        FieldRingsPhysical.push_back(new G4PVPlacement(0, G4ThreeVector(0, 0, PeekRodExtend-FR_thick / 2.0 -
                                                                              4 * (FR_thick + PEEK_Rod_thick) -
                                                                              2.5 * cm - 1.3 * cm - 0.7 * cm -
                                                                              1.3 * cm - 2 * cm - FR_thick -
                                                                              1 * (FR_thick + PEEK_Rod_thick)),
                                                       FR_logic, FR_logic->GetName(), gas_logic, 1, 13, checkOverlaps));
        FieldRingsPhysical.push_back(new G4PVPlacement(0, G4ThreeVector(0, 0, PeekRodExtend-FR_thick / 2.0 -
                                                                              4 * (FR_thick + PEEK_Rod_thick) -
                                                                              2.5 * cm - 1.3 * cm - 0.7 * cm -
                                                                              1.3 * cm - 2 * cm - FR_thick), FR_logic,
                                                       FR_logic->GetName(), gas_logic, 1, 14, checkOverlaps));

        // Reflections
        for (int i=0;i<FieldRingsPhysical.size();i++)
            new G4LogicalBorderSurface("SteelSurface_Rings_"+std::to_string(i),gas_phys,FieldRingsPhysical.at(i), OpSteelSurf);



        // PEEK Rods

        G4double x_rot_3 = 5 * std::sin(120 * deg) * cm;
        G4double y_rot_3 = -5 * std::cos(120 * deg) * cm;
        G4double x_rot_2 = 5 * std::sin(-120 * deg) * cm;
        G4double y_rot_2 = -5 * std::cos(-120 * deg) * cm;
        std::vector<G4double> x_rot_v = {0, x_rot_2, x_rot_3};
        std::vector<G4double> y_rot_v = {-5 * cm, y_rot_2, y_rot_3};
        G4double EL_PEEK_ROD_SHIFT = - 4 * (FR_thick + PEEK_Rod_thick) - 2.5 * cm - EL_thick - ElGap_ - EL_thick - 2 * cm - FR_thick;

        // Loop rotations
        for (G4int j = 0; j <=2; j++) {
            // The PEEK rods at the end
            new G4PVPlacement(0, G4ThreeVector(x_rot_v[j], y_rot_v[j], 1 * cm / 2.0 + 5 * (FR_thick + PEEK_Rod_thick)), PEEK_logic_cathode, PEEK_logic->GetName(), gas_logic, 1, j, checkOverlaps);

            // The other FC PEEK rods
            for (G4int i = 4; i >= -4; i--) {
                new G4PVPlacement(0, G4ThreeVector(x_rot_v[j], y_rot_v[j], PEEK_Rod_thick / 2.0 + i * (FR_thick + PEEK_Rod_thick)), PEEK_logic, PEEK_logic->GetName(), gas_logic, 1, j+4, checkOverlaps);
            }

            // PEEK EL
            new G4PVPlacement(0, G4ThreeVector(x_rot_v[j], y_rot_v[j], 1.1 * cm / 2.0       + EL_PEEK_ROD_SHIFT+PeekRodExtend), PEEK_logic_buffer, PEEK_logic->GetName(), gas_logic, 1, 0, checkOverlaps);
            new G4PVPlacement(0, G4ThreeVector(x_rot_v[j], y_rot_v[j], PEEK_Rod_thick / 2.0 + EL_PEEK_ROD_SHIFT - 1 * (FR_thick + PEEK_Rod_thick)+PeekRodExtend), PEEK_logic, PEEK_logic->GetName(), gas_logic, 1, j+1, checkOverlaps);
            new G4PVPlacement(0, G4ThreeVector(x_rot_v[j], y_rot_v[j], PEEK_Rod_thick / 2.0 + EL_PEEK_ROD_SHIFT - 2 * (FR_thick + PEEK_Rod_thick)+PeekRodExtend), PEEK_logic, PEEK_logic->GetName(), gas_logic, 1, j+2, checkOverlaps);
            new G4PVPlacement(0, G4ThreeVector(x_rot_v[j], y_rot_v[j], 3.37 * cm / 2.0      + EL_PEEK_ROD_SHIFT - 2 * (FR_thick + PEEK_Rod_thick) - FR_thick - 3.37 * cm+PeekRodExtend/2), PEEK_logic_buffer_end, PEEK_logic->GetName(), gas_logic, 1, j+3, checkOverlaps);
        }

#ifndef With_GarField
        // EL_Gap
         new G4PVPlacement(0, G4ThreeVector(0., 0., EL_pos), EL_logic, EL_solid->GetName(), gas_logic, 0, 0, checkOverlaps);
#endif
        G4VPhysicalVolume *EL_Ring_Plus = new G4PVPlacement(0, G4ThreeVector(0., 0., PeekRodExtend+EL_thick / 2.0 - FR_thick -
                                                                                     4 * (FR_thick + PEEK_Rod_thick) -
                                                                                     2.5 * cm - EL_thick),
                                                            EL_ring_logic, "EL_Ring_Plus", gas_logic, 0, 0, checkOverlaps);
        new G4LogicalBorderSurface("SteelSurfaceELPRing", gas_phys,EL_Ring_Plus, OpSteelSurf);

///// Off for Opticks /////
        // Place the Mesh bits
        /*4VPhysicalVolume *EL_Mesh_Plus_plus = new G4PVPlacement(rotateMesh, G4ThreeVector(0., 0.,PeekRodExtend+
                                                                                           EL_thick / 2.0 - FR_thick -
                                                                                           4 *
                                                                                           (FR_thick + PEEK_Rod_thick) -
                                                                                           2.5 * cm - EL_thick -
                                                                                           EL_thick / 2.0),
                                                                 ELP_Disk_logic, ELP_Disk_logic->GetName(), gas_logic,
                                                                 0, 0, checkOverlaps);
        HexCreator->PlaceHexagons(nHole, EL_hex_size, EL_mesh_thick, ELP_Disk_logic, EL_Hex_logic);
         */
///// Off for Opticks End /////
        G4VPhysicalVolume *EL_Ring_Plus_plus = new G4PVPlacement(0, G4ThreeVector(0., 0., PeekRodExtend+EL_thick / 2.0 - FR_thick -
                                                                                          4 *
                                                                                          (FR_thick + PEEK_Rod_thick) -
                                                                                          2.5 * cm - EL_thick - ElGap_ -
                                                                                          EL_thick), EL_ring_logic,
                                                                 "EL_Ring_Plus_plus", gas_logic, 0, 0, checkOverlaps);
        new G4LogicalBorderSurface("SteelSurfaceELRing", gas_phys,EL_Ring_Plus_plus, OpSteelSurf);

///// Off for Opticks /////
        // Place the Mesh bits
        /*G4VPhysicalVolume *EL_Mesh_Plus = new G4PVPlacement(0, G4ThreeVector(0., 0., PeekRodExtend+EL_thick / 2.0 - FR_thick -4 * (FR_thick + PEEK_Rod_thick) -2.5 * cm - EL_thick - ElGap_ -EL_thick + EL_thick / 2.0),ELPP_Disk_logic, ELPP_Disk_logic->GetName(), gas_logic, 0,0, checkOverlaps);
        HexCreator->PlaceHexagons(nHole, EL_hex_size, EL_mesh_thick, ELPP_Disk_logic, EL_Hex_logic);
         */
///// Off for Opticks End /////


        // Cathode
        G4VPhysicalVolume * Cathode = new G4PVPlacement(0, G4ThreeVector(0., 0., EL_thick / 2.0 + 1 * cm +
                                                                                 5 * (FR_thick + PEEK_Rod_thick)),
                                                        Cathode_ring_logic, "CATHODE", gas_logic, 0, 0, checkOverlaps);
        new G4LogicalBorderSurface("SteelSurfaceFR", gas_phys,Cathode, OpSteelSurf);
///// Off for Opticks /////
        // Place the Mesh bits
        /*G4VPhysicalVolume *Cathode_EL_Mesh = new G4PVPlacement(rotateMesh, G4ThreeVector(0., 0.,
                                                                                         EL_thick / 2.0 + 1 * cm + 5 *
                                                                                                                   (FR_thick +
                                                                                                                    PEEK_Rod_thick) -
                                                                                         EL_thick / 2.0),
                                                               Cathode_Disk_logic, Cathode_Disk_logic->GetName(),
                                                               gas_logic, 0, 0, checkOverlaps);
        HexCreator->PlaceHexagons(nHole, EL_hex_size, EL_mesh_thick, Cathode_Disk_logic, EL_Hex_logic);
        new G4LogicalSkinSurface("GAS_ELPMESH_OPSURF", ELP_Disk_logic, OpSteelSurf);
        new G4LogicalSkinSurface("GAS_ELPPMESH_OPSURF", ELPP_Disk_logic, OpSteelSurf);
        new G4LogicalSkinSurface("GAS_CATHODEMESH_OPSURF",Cathode_Disk_logic, OpSteelSurf);
         */
///// Off for Opticks End /////


        // --- Optical ---
        //new G4LogicalSkinSurface("SteelSurfaceFR", FR_logic, OpSteelSurf);
        //new G4LogicalSkinSurface("SteelSurfaceELRing1", EL_ring_logic, OpSteelSurf);
        //new G4LogicalSkinSurface("SteelSurfaceCathodeRing", Cathode_ring_logic, OpSteelSurf);



        //  ----------------------- Needle Source -----------------G4LogicalSkinSurfaceG4LogicalSkinSurface------------
        // Source





        // --- Placement ---

        // Source Holder
        G4VPhysicalVolume *Needle_Phys;
        if (!HideSourceHolder_) {
            /// Needle Source
            G4double NeedleyepRMin = 0;
            G4double NeedleyepRMax = (0.42) * mm;
            G4double NeedleyepDz = (2 / 2) * mm;
            G4double NeedleHalfLength = (2.56 / 2) * cm;
            G4double NeedleTailDiam = (0.6 / 2) * mm;
            G4double NeedleOffset = 1 * mm;

            G4Tubs *NeedleEye = new G4Tubs("NeedleEye", NeedleyepRMin, NeedleyepRMax, NeedleyepDz, 0., twopi);
            G4Tubs *NeedleTail = new G4Tubs("NeedleTail", NeedleyepRMin, NeedleTailDiam, NeedleHalfLength, 0., twopi);


            // Combining them to create the Needle
            G4VSolid *Needle = new G4UnionSolid("Needle", NeedleEye, NeedleTail, 0, G4ThreeVector(0, 0, NeedleHalfLength));


            G4LogicalVolume *Needle_Logic = new G4LogicalVolume(Needle, Steel, "Needle");
            // Needle Solid
            G4RotationMatrix *NeedleRotate = new G4RotationMatrix();
            NeedleRotate->rotateY(90. * deg);

            G4ThreeVector NeedlePos = {vtx_[0]- NeedleOffset, vtx_[1] , vtx_[2] - FieldCagePos / 2};
            G4ThreeVector CollPosition = {NeedlePos[0]- 5 * mm, NeedlePos[1] , NeedlePos[2]};
#ifndef With_GarField
            Needle_Phys = new G4PVPlacement(NeedleRotate, NeedlePos, Needle_Logic, Needle->GetName(), FieldCage_Logic, true,0, checkOverlaps);
#else
            Needle_Phys = new G4PVPlacement(NeedleRotate, NeedlePos, Needle_Logic, Needle->GetName(), gas_logic, true,0, checkOverlaps);
#endif
            new G4LogicalBorderSurface("SteelSurface_Needle", gas_phys, Needle_Phys, OpSteelSurf);
            // Collimator
            if (!HideCollimator_) {
                G4Tubs *SourceHolChamber_solid = new G4Tubs("SourceHolChamber", SourceEn_holedia / 2,
                                                            (SourceEn_diam / 2. + SourceEn_thickn),
                                                            (SourceEn_length / 2. + SourceEn_thickn), 0, twopi);
                G4LogicalVolume *SourceHolChamber_logic = new G4LogicalVolume(SourceHolChamber_solid, Steel,
                                                                              "SourceHolChamber_logic");

                G4Tubs *SourceHolChamberBlock_solid = new G4Tubs("SourceHolChBlock", 0, (SourceEn_holedia / 2),
                                                                 (SourceEn_thickn / 2), 0., twopi);
                G4LogicalVolume *SourceHolChamberBlock_logic = new G4LogicalVolume(SourceHolChamberBlock_solid,
                                                                                   Steel,
                                                                                   "SourceHolChBlock_logic");
                G4Tubs *Collimator = new G4Tubs("Collimator", (5.74 / 2) * mm, (11.5 / 2) * mm, 0.5 * cm, 0., twopi);
                G4Tubs *CollimatorBlock = new G4Tubs("CollimatorBlock", NeedleTailDiam, (5.74 / 2) * mm, 0.2 * cm, 0., twopi);
                G4VSolid *CollimatorWithBlock = new G4UnionSolid("CollimatorWithBlock", Collimator, CollimatorBlock, 0,
                                                                 G4ThreeVector(0, 0, 0.5 * cm - 0.2 * cm));
                G4LogicalVolume *Coll_Logic = new G4LogicalVolume(CollimatorWithBlock, PEEK,"CollimatorWithBlock");

                new G4PVPlacement(NeedleRotate, CollPosition, Coll_Logic, CollimatorWithBlock->GetName(), gas_logic,
                                  true, 0, checkOverlaps);
            }

        }



        //  --------------------------- PMT -----------------------------------
        //Adding the PMTs in here
        pmt1_ = new PmtR7378A();
        pmt1_->Construct();


        // Adding Logical Volumes for PMTs
        G4LogicalVolume *pmt1_logic = pmt1_->GetLogicalVolume();


        // // PMTs
        G4double PMT_offset = 0.2 * cm;
        G4double PMT_pos =
                (chamber_length / 2) + chamber_thickn + (pmt1_->Length() / 2) + MgF2_window_thickness_ + PMT_offset;
        G4RotationMatrix *pmt1rotate = new G4RotationMatrix();
        pmt1rotate->rotateX(180. * deg);


        //// PMT Covering Tube ///
        G4double offset = 1.65 * cm;
        G4double PMT_Tube_Length1 =MgF2_window_thickness_ + (pmt1_->Length() + 0.5 * cm) / 2 + offset - PMT_offset - 0.05 * cm;
        G4double PMT_Tube_Length0 = (17 * cm + pmt1_->Length()) / 2 - 3.87 * cm - PMT_offset;
        G4double PMT_Tube_Block_Thickness = 0.2 * cm;
        G4double LongPMTTubeOffset = 7.5 * cm - 3.9 * cm;
        G4double PMTTubeDiam = 2.54 * cm;

        // Tube Away from EL
        G4Tubs *PMT_Tube_solid0 = new G4Tubs("PMT_TUBE0", (PMTTubeDiam / 2) + 0.5 * cm, (PMTTubeDiam / 2) + 0.7 * cm,
                                             PMT_Tube_Length0, 0, twopi);
        G4LogicalVolume *PMT_Tube_Logic0 = new G4LogicalVolume(PMT_Tube_solid0, Steel,
                                                               PMT_Tube_solid0->GetName());

        G4Tubs *PMT_Block_solid0 = new G4Tubs("PMT_TUBE_BLOCK0", 0, (PMTTubeDiam / 2 + 0.5 * cm),
                                              PMT_Tube_Block_Thickness, 0, twopi);
        G4LogicalVolume *PMT_Block_Logic0 = new G4LogicalVolume(PMT_Block_solid0, Steel,
                                                                PMT_Block_solid0->GetName());

        // Vacuum for PMT TUBE0
        /// add 2.54*cm to config file
        G4Tubs *InsideThePMT_Tube_solid0 = new G4Tubs("PMT_TUBE_VACUUM0", 0, (PMTTubeDiam / 2 + 0.5 * cm),
                                                      PMT_Tube_Length0, 0, twopi);
        G4LogicalVolume *InsideThePMT_Tube_Logic0 = new G4LogicalVolume(InsideThePMT_Tube_solid0, vacuum,
                                                                        InsideThePMT_Tube_solid0->GetName());

        // Tube Close to EL
        G4Tubs *PMT_Tube_solid1 = new G4Tubs("PMT_TUBE1", (PMTTubeDiam / 2) + 0.5 * cm, (PMTTubeDiam / 2) + 0.7 * cm,
                                             PMT_Tube_Length1, 0, twopi);
        G4LogicalVolume *PMT_Tube_Logic1 = new G4LogicalVolume(PMT_Tube_solid1, Steel,
                                                               PMT_Tube_solid1->GetName());
        G4Tubs *PMT_Block_solid1 = new G4Tubs("PMT_TUBE_BLOCK1", 0, (PMTTubeDiam / 2 + 0.5 * cm),
                                              PMT_Tube_Block_Thickness, 0, twopi);
        G4LogicalVolume *PMT_Block_Logic = new G4LogicalVolume(PMT_Block_solid1, Steel,
                                                               PMT_Block_solid1->GetName());

        // Vacuum for PMT TUBE1

        G4Tubs *InsideThePMT_Tube_solid1 = new G4Tubs("PMT_TUBE_VACUUM1", 0, (PMTTubeDiam / 2 + 0.5 * cm),
                                                      PMT_Tube_Length1, 0, twopi);
        G4LogicalVolume *InsideThePMT_Tube_Logic1 = new G4LogicalVolume(InsideThePMT_Tube_solid1, vacuum,
                                                                        InsideThePMT_Tube_solid1->GetName());


        // --- Placement ---

        // PMT Tubes
        G4VPhysicalVolume *PMT_Tube_Phys0 = new G4PVPlacement(0, G4ThreeVector(0, 0, PMT_pos + LongPMTTubeOffset),
                                                              PMT_Tube_Logic0, PMT_Tube_Logic0->GetName(),
                                                              lab_logic_volume, false, 0, checkOverlaps);
        G4VPhysicalVolume *PMT_Tube_Phys1 = new G4PVPlacement(0, G4ThreeVector(0, 0, -(PMT_pos - PMT_offset) - offset),
                                                              PMT_Tube_Logic1, PMT_Tube_Logic1->GetName(),
                                                              lab_logic_volume, false, 0, checkOverlaps);

        // PMT Tube Vacuum
        G4VPhysicalVolume *PMT_Tube_Vacuum_Phys0 = new G4PVPlacement(0,
                                                                     G4ThreeVector(0, 0, PMT_pos + LongPMTTubeOffset+300*um),
                                                                     InsideThePMT_Tube_Logic0, "PMT_TUBE_VACUUM_CATH",
                                                                     lab_logic_volume, false, 0, checkOverlaps);

        G4VPhysicalVolume *PMT_Tube_Vacuum_Phys1 = new G4PVPlacement(0, G4ThreeVector(0, 0,
                                                                                      -(PMT_pos - PMT_offset) - offset),
                                                                     InsideThePMT_Tube_Logic1, "PMT_TUBE_VACUUM_ANODE",
                                                                     lab_logic_volume, false, 0, checkOverlaps);

        // PMT Tube Block
        new G4PVPlacement(0, G4ThreeVector(0, 0,
                                           PMT_pos - PMT_offset + PMT_Tube_Length0 - PMT_Tube_Block_Thickness / 2 +
                                           LongPMTTubeOffset), PMT_Block_Logic0, PMT_Block_Logic0->GetName(),
                          lab_logic_volume, false, 0, checkOverlaps);
        new G4PVPlacement(0, G4ThreeVector(0, 0,
                                           -(PMT_pos - PMT_offset + PMT_Tube_Length1 - PMT_Tube_Block_Thickness / 2) -
                                           offset), PMT_Block_Logic, PMT_Block_Logic->GetName(), lab_logic_volume,
                          false, 1, checkOverlaps);



        // --- Optical ---
        new G4LogicalBorderSurface("SteelSurface_Camera_Enclosing",PMT_Tube_Vacuum_Phys0 ,PMT_Tube_Phys0, OpSteelSurf);
        new G4LogicalBorderSurface("SteelSurface_PMT_Enclosing",PMT_Tube_Vacuum_Phys1,PMT_Tube_Phys1, OpSteelSurf);


        //  ------------------------ Camera -----------------------------------

        // Window
        G4double camHalfLength = 0.5 * mm;
        G4double camRadius = 12.7 * mm;
        G4VSolid *camSolid = new G4Tubs("camWindow", 0., camRadius, camHalfLength, 0., twopi);
        G4LogicalVolume *camLogical = new G4LogicalVolume(camSolid, MgF2, "camLogical");



        // --- Placement ---
        G4double ImageDist = 7.945 * cm; // Got from trial and error
        G4VPhysicalVolume *camPhysical = new G4PVPlacement(0, G4ThreeVector(0, 0, (chamber_length / 2 + chamber_thickn +ImageDist) - PMT_pos -LongPMTTubeOffset), camLogical,"camWindow", InsideThePMT_Tube_Logic0, false, 0, checkOverlaps);

        // PMT Placement
        G4VPhysicalVolume *pmt1_phys =new G4PVPlacement(0, G4ThreeVector(0, 0., (PMT1_Pos_ - pmt1_->Length() / 2 - MgF2_window_thickness_ / 2)),pmt1_logic, "S1", InsideThePMT_Tube_Logic1, false, 0, checkOverlaps);



        //  ------------------------ EL Brackets ------------------------------

        G4Box *bracket_body = new G4Box("Box_body", 26 * mm / 2.0, 26 * mm / 2.0, 51 * mm / 2.0);
        G4Box *bracket_subtract = new G4Box("Box_subtract", 27 * mm / 2.0, 26 * mm / 2.0, 33 * mm / 2.0);

        G4SubtractionSolid *bracket_solid = new G4SubtractionSolid("Bracket", bracket_body, bracket_subtract, 0,
                                                                   G4ThreeVector(0, 3 * mm, 0));
        G4LogicalVolume *bracket_logical = new G4LogicalVolume(bracket_solid, teflon, "bracketLogical");


        // Define a rotation matrix to orient all detector pieces along y direction
        G4RotationMatrix *rotateZ_120 = new G4RotationMatrix();
        rotateZ_120->rotateZ(120. * deg);
        G4RotationMatrix *rotateZ_m120 = new G4RotationMatrix();
        rotateZ_m120->rotateZ(-120. * deg);

        x_rot_3 = 5.7 * std::sin(120 * deg) * cm;
        y_rot_3 = -5.7 * std::cos(120 * deg) * cm;
        x_rot_2 = 5.7 * std::sin(-120 * deg) * cm;
        y_rot_2 = -5.7 * std::cos(-120 * deg) * cm;

        // --- Placement ---

        G4VPhysicalVolume *bracketPhysical1 = new G4PVPlacement(0, G4ThreeVector(0, (-5.7) * cm, EL_pos),
                                                                bracket_logical, "bracketPhysical", gas_logic, true, 0,
                                                                checkOverlaps);
        G4VPhysicalVolume *bracketPhysical2 = new G4PVPlacement(rotateZ_120, G4ThreeVector(x_rot_2, y_rot_2, EL_pos),
                                                                bracket_logical, "bracketPhysical", gas_logic, true, 1,
                                                                checkOverlaps);
        G4VPhysicalVolume *bracketPhysical3 = new G4PVPlacement(rotateZ_m120, G4ThreeVector(x_rot_3, y_rot_3, EL_pos),
                                                                bracket_logical, "bracketPhysical", gas_logic, true, 2,checkOverlaps);

#ifndef With_GarField
        // Electrical Field
      if(true){
            FielCageGap=(160.3+29.55)*mm;
            FieldCagePos=chamber_length/2-((129)*mm)-FielCageGap/2-ElGap_/2;
            EL_pos=chamber_length/2-FielCageGap/2-((326)*mm)-ElGap_/2;

            UniformElectricDriftField* field = new UniformElectricDriftField();

            field->SetCathodePosition(FieldCagePos/2+FielCageGap/2);
            field->SetAnodePosition(EL_pos/2+ElGap_/2);
            //field->SetAnodePosition(EL_pos/2);
            field->SetDriftVelocity(.90*mm/microsecond);
            field->SetTransverseDiffusion(.92*mm/sqrt(cm));
            field->SetLongitudinalDiffusion(.36*mm/sqrt(cm));
            /*if(!HideSourceHolder_){
                if(!HideCollimator_) field->SetStepLimit(3*mm,0.3*mm);
                else field->SetStepLimit(1*mm,0.2*mm);
            }*/

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
            //EfieldForEL->SetELGap(ElGap_*cm);
	        EfieldForEL->SetLightYield(ELyield_);
            G4Region* el_region = new G4Region("EL_GAP");
            el_region->SetUserInformation(EfieldForEL);
            el_region->AddRootLogicalVolume(EL_logic);

        }
#endif

/* Not needed
        // Add the camera as a sensitive detector
        SensorSD* camerasd = new SensorSD("/Camera/Cm");
        camerasd->SetDetectorVolumeDepth(1);
        camerasd->SetTimeBinning(100*ns);
        camerasd->SetDetectorNamingOrder(1);
        camerasd->SetMotherVolumeDepth(2);
        G4SDManager::GetSDMpointer()->AddNewDetector(camerasd);
        camLogical->SetSensitiveDetector(camerasd);
*/

        // Camera
        G4OpticalSurface *opXenon_Glass = new G4OpticalSurface("CamSurfaceBorder");
        opXenon_Glass->SetMaterialPropertiesTable(opticalprops::PerfectDetector());
        opXenon_Glass->SetModel(unified);                  // SetModel
        opXenon_Glass->SetType(dielectric_dielectric);   // SetType
        opXenon_Glass->SetFinish(polished);                 // SetFinish
        //new G4LogicalBorderSurface("CamSurfaceBorder",PMT_Tube_Vacuum_Phys0,camPhysical,opXenon_Glass);
        new G4LogicalSkinSurface("CamSurfaceBorder",camLogical,opXenon_Glass);




        // ____________________________________________________________________
        // ================= Detector Properties  =============================

        G4SDManager *SDManager = G4SDManager::GetSDMpointer();
        IonizationSD* ionisd = new IonizationSD("/CRAB0/GAS");
        SDManager->SetVerboseLevel(1);
        SDManager->AddNewDetector(ionisd);
#ifdef  With_GarField
        gas_logic->SetSensitiveDetector(ionisd);
#else
        FieldCage_Logic->SetSensitiveDetector(ionisd);

#endif
        gas_logic->SetSensitiveDetector(ionisd);
        // Construct a G4Region, connected to the logical volume in which you want to use the G4FastSimulationModel
        G4Region *regionGas = new G4Region("GasRegion");
        regionGas->AddRootLogicalVolume(gas_logic);
#ifdef With_GarField
        GarfieldHelper GH(chamber_diam/2.0/cm, chamber_length/cm, Active_diam/2.0/cm , FielCageGap/cm, gas_pressure_, ElGap_, fieldDrift_, fieldEL_);
        GH.SetGasFile(GasFile_);
#endif
        // Visuals
        AssignVisuals();

        //These commands generate the four gas models and connect it to the GasRegion
#ifdef With_GarField
        G4Region *region = G4RegionStore::GetInstance()->GetRegion("GasRegion");
        new DegradModel("DegradModel",region, GH);
        new GarfieldVUVPhotonModel("GarfieldVUVPhotonModel",region, GH, ionisd);
#endif
        this->SetLogicalVolume(lab_logic_volume);
        this->SetSpan(Lab_size*1.02);



    }

    void CRAB0::AssignVisuals() {
        // Chamber
        G4LogicalVolumeStore *lvStore = G4LogicalVolumeStore::GetInstance();


        // Lab
        G4LogicalVolume *Lab = lvStore->GetVolume("LAB");
        G4VisAttributes *LabVa = new G4VisAttributes(G4Colour(2, 2, 2));
        LabVa->SetForceWireframe(false);
        //Chamber
        G4LogicalVolume *Chamber = lvStore->GetVolume("CHAMBER");
        G4VisAttributes *ChamberVa = new G4VisAttributes(G4Colour(1, 1, 1));
        ChamberVa->SetForceSolid(true);
        Chamber->SetVisAttributes(G4VisAttributes::GetInvisible());


        //GAS
        G4LogicalVolume *Gas = lvStore->GetVolume("GAS");
        G4VisAttributes *GasVa = new G4VisAttributes(nexus::YellowAlpha());
        GasVa->SetForceCloud(true);
        Gas->SetVisAttributes(GasVa);

        /*
        //Source Enclosure Related
        G4LogicalVolume *SourceHolder = lvStore->GetVolume("SourceHolChamber_logic");
        G4LogicalVolume *Collimator = lvStore->GetVolume("CollimatorWithBlock");
        G4VisAttributes *CollimatorVa = new G4VisAttributes(nexus::YellowAlpha());
        CollimatorVa->SetForceSolid(true);
        G4LogicalVolume *SourceHolderBlock = lvStore->GetVolume("SourceHolChBlock_logic");
        G4VisAttributes *SourceHolderVa = new G4VisAttributes(G4Colour(2, 2, 2));
        SourceHolderVa->SetForceSolid(true);
        Collimator->SetVisAttributes(CollimatorVa);
        */
        G4LogicalVolume *Needle = lvStore->GetVolume("Needle");
        Needle->SetVisAttributes(ChamberVa);



        // Flange

        // Anode Flange
        G4LogicalVolume* flangeLog_anode = lvStore->GetVolume("CHAMBER_FLANGE_ANODE");
        G4VisAttributes flangeVis=nexus::DarkGreyAlpha();
        flangeVis.SetForceSolid(true);
        flangeLog_anode->SetVisAttributes(ChamberVa);

        G4LogicalVolume* flangeLog_cathode = lvStore->GetVolume("CHAMBER_FLANGE_CATHODE");
        flangeVis.SetForceSolid(true);
        flangeLog_cathode->SetVisAttributes(ChamberVa);

        // Field Rings
        G4LogicalVolume *FRLog = lvStore->GetVolume("FR");
        G4VisAttributes FReVis = nexus::CopperBrownAlpha();
        FReVis.SetForceSolid(true);
        FRLog->SetVisAttributes(FReVis);

        // EL Rings
        G4LogicalVolume *EL_RingLog = lvStore->GetVolume("EL_Ring");
        G4VisAttributes EL_RingVis = nexus::DarkGreyAlpha();
        EL_RingVis.SetForceSolid(true);
        EL_RingLog->SetVisAttributes(EL_RingVis);


        // Cathode Rings
        G4LogicalVolume *CATHODE = lvStore->GetVolume("Cathode_Ring");
        G4VisAttributes CATHODEVis = nexus::DarkGreyAlpha();
        CATHODEVis.SetForceSolid(true);
        CATHODE->SetVisAttributes(CATHODEVis);

        /*
        // Brackets
        G4LogicalVolume *BracketLog = lvStore->GetVolume("bracketLogical");
        G4VisAttributes BracketVis = nexus::DirtyWhiteAlpha();
        BracketVis.SetForceSolid(true);
        BracketLog->SetVisAttributes(BracketVis);
        */

        // PEEK
        G4LogicalVolume *PEEKLog = lvStore->GetVolume("PEEK_Rod");
        G4VisAttributes PEEKVis = nexus::YellowAlpha();
        PEEKVis.SetForceSolid(true);
        PEEKLog->SetVisAttributes(PEEKVis);

        PEEKLog = lvStore->GetVolume("PEEK_Rod_C");
        PEEKLog->SetVisAttributes(PEEKVis);

        PEEKLog = lvStore->GetVolume("PEEK_Rod_B");
        PEEKLog->SetVisAttributes(PEEKVis);

        PEEKLog = lvStore->GetVolume("PEEK_Rod_BE");
        PEEKLog->SetVisAttributes(PEEKVis);


        //PMT TUBE AND PMT BLOCK
        G4LogicalVolume *PmttubeLog0 = lvStore->GetVolume("PMT_TUBE0");
        PmttubeLog0->SetVisAttributes(G4VisAttributes::GetInvisible());
        G4LogicalVolume *PmttubeBlockLog0 = lvStore->GetVolume("PMT_TUBE_BLOCK0");
        G4LogicalVolume *PmttubeLog1 = lvStore->GetVolume("PMT_TUBE1");
        PmttubeLog1->SetVisAttributes(G4VisAttributes::GetInvisible());
        G4LogicalVolume *PmttubeBlockLog1 = lvStore->GetVolume("PMT_TUBE_BLOCK1");
        PmttubeBlockLog0->SetVisAttributes(ChamberVa);
        PmttubeBlockLog1->SetVisAttributes(ChamberVa);
        G4LogicalVolume *PmttubeVacuumLog1 = lvStore->GetVolume("PMT_TUBE_VACUUM0");
        G4LogicalVolume *PmttubeVacuumLog2 = lvStore->GetVolume("PMT_TUBE_VACUUM1");
        G4VisAttributes PmttubeVacuumVis = nexus::DarkGreen();
        PmttubeVacuumVis.SetForceCloud(true);
        PmttubeVacuumLog1->SetVisAttributes(PmttubeVacuumVis);
        PmttubeVacuumLog2->SetVisAttributes(PmttubeVacuumVis);


        //MgF2Window
        G4LogicalVolume *lensLogical = lvStore->GetVolume("Lens");
        G4VisAttributes MgF2LensVis = nexus::DarkGreen();
        MgF2LensVis.SetForceSolid(true);
        lensLogical->SetVisAttributes(MgF2LensVis);

        G4LogicalVolume *MgF2WindowLog = lvStore->GetVolume("MgF2_WINDOW");
        G4VisAttributes MgF2WindowVis = nexus::DarkGreen();
        MgF2WindowVis.SetForceSolid(true);
        MgF2WindowLog->SetVisAttributes(MgF2WindowVis);

        /* // Lens
        G4LogicalVolume * LensLog=lvStore->GetVolume("LensMotherPV");
        G4VisAttributes  LensVis=Red();
        LensVis.SetForceSolid(true);
        LensLog->SetVisAttributes(LensVis);
        */

        // Camera
        G4LogicalVolume *CAMLog = lvStore->GetVolume("camLogical");
        G4VisAttributes CAMVis = nexus::DarkRedAlpha();
        CAMVis.SetForceSolid(true);
        CAMLog->SetVisAttributes(CAMVis);
#ifndef With_GarField
        // EL-Region
        G4LogicalVolume *ELLogic = lvStore->GetVolume("EL_GAP");
        G4VisAttributes ELVis = nexus::BlueAlpha();
        ELVis.SetForceCloud(true);
        //ELVis.SetForceSolid(true);
        ELLogic->SetVisAttributes(ELVis);
        // FieldCage
        G4LogicalVolume * FieldCage=lvStore->GetVolume("ACTIVE");
        G4VisAttributes FielCageVis=nexus::Red();
        FielCageVis.SetForceCloud(true);
        FieldCage->SetVisAttributes(FielCageVis);
#endif
        // Cathode Grid
        // G4VisAttributes mesh_col = nexus::DarkGrey();
        // mesh_col.SetForceSolid(true);
        // G4LogicalVolume* Meshlog = lvStore->GetVolume("ELP_Mesh_Logic");
        // Meshlog->SetVisAttributes(G4VisAttributes::GetInvisible());
        // Meshlog = lvStore->GetVolume("ELPP_Mesh_Logic");
        // Meshlog->SetVisAttributes(G4VisAttributes::GetInvisible());
        // Meshlog = lvStore->GetVolume("Cathode_Mesh_Logic");
        // Meshlog->SetVisAttributes(G4VisAttributes::GetInvisible());

        // mesh_col = nexus::Empty();
        // mesh_col.SetForceSolid(true);
        // Meshlog = lvStore->GetVolume("Mesh_Hex");
        // Meshlog->SetVisAttributes(G4VisAttributes::GetInvisible());

        /*// ACTIVE
        G4LogicalVolume *ACTIVE = lvStore->GetVolume("ACTIVE");
        G4VisAttributes ACTIVEVis = nexus::Red();
        ACTIVEVis.SetForceCloud(true);
        ACTIVE->SetVisAttributes(G4VisAttributes::GetInvisible());
        */

        Lab->SetVisAttributes(G4VisAttributes::GetInvisible());

    }



    G4ThreeVector CRAB0::GenerateVertex(const G4String& region) const{
        G4ThreeVector vertex(0., 0., 0.);

        if (region == "CENTER") {
            vertex = G4ThreeVector(0., 0., 0.);
        }

        else if (region == "AD_HOC") {
            // AD_HOC does not need to be shifted because it is passed by the user
            vertex = specific_vertex_;
            return vertex;
        }

        else if (region == "NEEDLE") {
            // AD_HOC does not need to be shifted because it is passed by the user
            vertex = specific_vertex_+vtx_;
            return vertex;
        }


        else {
            G4Exception("[CRAB0]", "GenerateVertex()", FatalException,
            "Unknown vertex generation region!");
        }

        return vertex;
    }


}
