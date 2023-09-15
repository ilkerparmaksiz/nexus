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
#include "DegradModel.h"
#include "GarfieldVUVPhotonModel.h"
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
#include "GarfieldHelper.h"


namespace nexus {
    
    REGISTER_CLASS(CRAB0,GeometryBase)
    
    CRAB0::CRAB0() :
            checkOverlaps(1),
            temperature(300 * kelvin), // temperature
            Lab_size(3. * m),
            chamber_diam(16.4 * cm),
            chamber_length(43.18 * cm), // Config files vary
            chamber_thickn(7. * mm),
            SourceEn_offset(5.7 * cm),
            SourceEn_diam(1.0 * cm),
            SourceEn_length(1 * cm),
            SourceEn_thickn(2. * mm),
            SourceEn_holedia(5. * mm),
            gas_pressure_(10 * bar),
            vtx_(0 * cm, -1.6 * cm, -5.25 * cm),
            Active_diam(8.6 * cm),
            sc_yield_(25510. / MeV),
            e_lifetime_(1000. * ms),
            MgF2_window_thickness_(6. * mm),
            MgF2_window_diam_(16.5 * mm),
            HideSourceHolder_(false),
            max_step_size_(1. * mm),
            ElGap_(7 * mm),
            ELyield_(970 / cm),
            PMT1_Pos_(2.32 * cm),
            PMT3_Pos_(3.52 * cm),
            HideCollimator_(true),
            specific_vertex_{} {

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


        }

    CRAB0::~CRAB0() {
    }

    void CRAB0::Construct() {

        //  ------------------------ Materials --------------------------------
        G4Material *gxe    = materials::GXe(gas_pressure_, 68);
        G4Material *MgF2   = materials::MgF2();
        G4Material *Steel  = materials::Steel();
        G4Material *PEEK   = materials::PEEK();
        G4Material *vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
        G4Material *teflon = G4NistManager::Instance()->FindOrBuildMaterial("G4_TEFLON");

        //  ------------------- Optical Properties ----------------------------
        MgF2  ->SetMaterialPropertiesTable(opticalprops::MgF2());
        vacuum->SetMaterialPropertiesTable(opticalprops::Vacuum());
        gxe   ->SetMaterialPropertiesTable(opticalprops::GXe(gas_pressure_, 68, sc_yield_, e_lifetime_));

        // OpticalSurface
        G4OpticalSurface *OpSteelSurf = new G4OpticalSurface("SteelSurface", unified, polished, dielectric_metal);
        OpSteelSurf->SetMaterialPropertiesTable(opticalprops::STEEL());

        //  ----------------------- Lab Space ---------------------------------
        G4String lab_name = "LAB";
        G4Box *lab_solid_volume = new G4Box(lab_name, Lab_size / 2, Lab_size / 2, Lab_size / 2);
        G4LogicalVolume *lab_logic_volume = new G4LogicalVolume(lab_solid_volume,
                                                                G4NistManager::Instance()->FindOrBuildMaterial(
                                                                        "G4_AIR"), lab_name);

        this->SetLogicalVolume(lab_logic_volume);


        // Xenon Gas
        G4Tubs *gas_solid = new G4Tubs("GAS", 0., chamber_diam / 2., chamber_length / 2. + chamber_thickn, 0., twopi);
        gas_logic = new G4LogicalVolume(gas_solid, gxe, "GAS");


        // Define a rotation matrix to orient all detector pieces along y direction
        G4RotationMatrix *rotateX = new G4RotationMatrix();
        rotateX->rotateX(90. * deg);

        // Define a rotation matrix for the meshes so they point in the right direction
        G4RotationMatrix *rotateMesh = new G4RotationMatrix();
        rotateMesh->rotateZ(30. * deg);


        // --- Placement ---
        auto labPhysical = new G4PVPlacement(0, G4ThreeVector(), lab_logic_volume, lab_logic_volume->GetName(), 0,
                                             false, 0, false);


        // Xenon Gas in Active Area and Non-Active Area
        G4VPhysicalVolume *gas_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), gas_logic, gas_solid->GetName(),
                                                        lab_logic_volume, false, 0, false);

        
        // ____________________________________________________________________
        // ==================== Detector Geometry =============================
        //  ------------------------ Vessel -----------------------------------
    
        // Steel end caps with holes
        G4Tubs *chamber_flange_solid = new G4Tubs("CHAMBER_FLANGE", MgF2_window_diam_ / 2,
                                                  (chamber_diam / 2. + chamber_thickn), chamber_thickn / 2.0, 0.,
                                                  twopi);
        G4LogicalVolume *chamber_flange_logic = new G4LogicalVolume(chamber_flange_solid, materials::Steel(),
                                                                    "CHAMBER_FLANGE");

        // Chamber barrel
        G4Tubs *chamber_solid = new G4Tubs("CHAMBER", chamber_diam / 2., (chamber_diam / 2. + chamber_thickn),
                                           (chamber_length / 2), 0., twopi);
        G4LogicalVolume *chamber_logic = new G4LogicalVolume(chamber_solid, materials::Steel(), "CHAMBER");

        // --- Placement ---

        // Flanges on the Chamber, place in the gas logic so we include the aperature region
        G4VPhysicalVolume *Left_Flange_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, chamber_length / 2 +
                                                                                       chamber_thickn / 2.0),
                                                                chamber_flange_logic, chamber_flange_solid->GetName(),
                                                                gas_logic, true, 0, false);
        G4VPhysicalVolume *Right_Flange_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, -chamber_length / 2 -
                                                                                        chamber_thickn / 2.0),
                                                                 chamber_flange_logic, chamber_flange_solid->GetName(),
                                                                 gas_logic, true, 1, false);

        G4VPhysicalVolume *chamber_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., 0), chamber_logic,
                                                            chamber_solid->GetName(), lab_logic_volume, false, 0,
                                                            false);

        // --- Optical ---
        new G4LogicalBorderSurface("SteelSurface_Chamber", gas_phys, chamber_phys, OpSteelSurf);
        new G4LogicalBorderSurface("SteelSurface_LeftFlange", gas_phys, Left_Flange_phys, OpSteelSurf);
        new G4LogicalBorderSurface("SteelSurface_RightFlange", gas_phys, Right_Flange_phys, OpSteelSurf);


        //  --------------------- Window/lens ---------------------------------
        // MgF2 window
        G4Tubs *MgF2_window_solid = new G4Tubs("MgF2_WINDOW", 0., MgF2_window_diam_ / 2.,
                                               (MgF2_window_thickness_) / 2., 0., twopi);
        G4LogicalVolume *MgF2_window_logic = new G4LogicalVolume(MgF2_window_solid, MgF2, "MgF2_WINDOW");

        // Lens
        const G4double lensRcurve(2.83 * cm); // radius of curvature of MgF2 Lens
        const G4ThreeVector posLensTubeIntersect(0., 0., -lensRcurve);

        // Lens is made from the intersection of a sphere and a cylinder
        G4double maxLensLength = 4 * mm;
        G4Tubs *sLensTube = new G4Tubs("sLensSphereTube", 0, MgF2_window_diam_ / 2, maxLensLength, 0.,
                                       twopi); // 4 mm is the max lens length
        G4Orb *sLensOrb = new G4Orb("sLensSphere", lensRcurve);
        G4IntersectionSolid *sLens = new G4IntersectionSolid("sLens", sLensTube, sLensOrb, 0, posLensTubeIntersect);

        // Lens logical
        G4LogicalVolume *lensLogical = new G4LogicalVolume(sLens, MgF2, "Lens");

        // --- Placement ---
        G4double window_posz = chamber_length / 2 + chamber_thickn;
        // G4VPhysicalVolume* lensPhysical = new G4PVPlacement(0, G4ThreeVector(0., 0., window_posz), MgF2_window_logic,"MgF2_WINDOW1", lab_logic_volume,false, 0, false);

        G4VPhysicalVolume *lensPhysical = new G4PVPlacement(0, G4ThreeVector(0., 0., window_posz + maxLensLength / 2.0),
                                                            lensLogical, "MgF2_WINDOW1", gas_logic, false, 0, false);
        new G4PVPlacement(0, G4ThreeVector(0., 0., -window_posz), MgF2_window_logic, "MgF2_WINDOW2", gas_logic, false,
                          1, false);

        // --- Optical ---
        G4OpticalSurface *opXenon_Glass2 = new G4OpticalSurface("XenonLensSurface");
        opXenon_Glass2->SetModel(glisur);                  // SetModel
        opXenon_Glass2->SetType(dielectric_dielectric);   // SetType
        opXenon_Glass2->SetFinish(polished);                 // SetFinish
        opXenon_Glass2->SetPolish(0.0);
        new G4LogicalBorderSurface("XenonLensSurface", gas_phys, lensPhysical, opXenon_Glass2);

        //  --------------------------- EL ------------------------------------

        FielCageGap = 21.26 * cm;
        
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


        // EL Mesh
        // Dist from centre of hex to hex vertex, excluding the land width (circumradius)
        G4double hex_circumR = EL_hex_size / std::sqrt(3);

        // Number of hexagons needed -- need to use fixed amount, too many and nexus will crash
        G4int nHole = 16;

        // Define the Stainless steel mesh cylinder to subtract hex pattern from
        G4Tubs *Mesh_Disk = new G4Tubs("Mesh_Disk", 0., EL_OD / 2.0, EL_mesh_thick / 2., 0.,
                                       twopi); // Use OD so mesh stays within the logical

        HexagonMeshTools::HexagonMeshTools *HexCreator; // Hexagonal Mesh Tool

        // Define a hexagonal prism
        G4ExtrudedSolid *HexPrism = HexCreator->CreateHexagon(EL_mesh_thick, hex_circumR);

        G4LogicalVolume *ELP_Disk_logic     = new G4LogicalVolume(Mesh_Disk, Steel, "ELP_Mesh_Logic");
        G4LogicalVolume *ELPP_Disk_logic    = new G4LogicalVolume(Mesh_Disk, Steel, "ELPP_Mesh_Logic");
        G4LogicalVolume *Cathode_Disk_logic = new G4LogicalVolume(Mesh_Disk, Steel, "Cathode_Mesh_Logic");
        G4LogicalVolume *EL_Hex_logic       = new G4LogicalVolume(HexPrism,  gxe,   "Mesh_Hex");

        //  ----------------------- Field Cage --------------------------------
        
        // FieldCage -- needs to be updated to rings and PEEK rods
        G4Tubs *FieldCage_Solid = new G4Tubs("FIELDCAGE", 0., Active_diam / 2., FielCageGap / 2, 0., twopi);
        G4LogicalVolume *FieldCage_Logic = new G4LogicalVolume(FieldCage_Solid, gxe, "FIELDCAGE");

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

        G4Tubs *PEEK_Rod_Solid_buffer_end = new G4Tubs("PEEK_Rod", 0., PEEK_Rod_OD / 2., 3.37 * cm / 2.0, 0., twopi);
        G4LogicalVolume *PEEK_logic_buffer_end = new G4LogicalVolume(PEEK_Rod_Solid_buffer_end, PEEK, "PEEK_Rod_BE");

        // FieldCage
        G4double FieldCagePos = -0 * cm;
        G4double EL_pos = -10.98 * cm;

        // --- Placement ---

        // Field Cage Rings
        G4VPhysicalVolume *FR_FC_10 = new G4PVPlacement(0, G4ThreeVector(0, 0, -FR_thick / 2.0 +
                                                                               5 * (FR_thick + PEEK_Rod_thick)),
                                                        FR_logic, FR_logic->GetName(), gas_logic, 0, 0, false);
        G4VPhysicalVolume *FR_FC_9 = new G4PVPlacement(0, G4ThreeVector(0, 0, -FR_thick / 2.0 +
                                                                              4 * (FR_thick + PEEK_Rod_thick)),
                                                       FR_logic, FR_logic->GetName(), gas_logic, 0, 0, false);
        G4VPhysicalVolume *FR_FC_8 = new G4PVPlacement(0, G4ThreeVector(0, 0, -FR_thick / 2.0 +
                                                                              3 * (FR_thick + PEEK_Rod_thick)),
                                                       FR_logic, FR_logic->GetName(), gas_logic, 0, 0, false);
        G4VPhysicalVolume *FR_FC_7 = new G4PVPlacement(0, G4ThreeVector(0, 0, -FR_thick / 2.0 +
                                                                              2 * (FR_thick + PEEK_Rod_thick)),
                                                       FR_logic, FR_logic->GetName(), gas_logic, 0, 0, false);
        G4VPhysicalVolume *FR_FC_6 = new G4PVPlacement(0, G4ThreeVector(0, 0, -FR_thick / 2.0 +
                                                                              1 * (FR_thick + PEEK_Rod_thick)),
                                                       FR_logic, FR_logic->GetName(), gas_logic, 0, 0, false);
        G4VPhysicalVolume *FR_FC_5 = new G4PVPlacement(0, G4ThreeVector(0, 0, -FR_thick / 2.0), FR_logic,
                                                       FR_logic->GetName(), gas_logic, 0, 0, false);
        G4VPhysicalVolume *FR_FC_4 = new G4PVPlacement(0, G4ThreeVector(0, 0, -FR_thick / 2.0 -
                                                                              1 * (FR_thick + PEEK_Rod_thick)),
                                                       FR_logic, FR_logic->GetName(), gas_logic, 0, 0, false);
        G4VPhysicalVolume *FR_FC_3 = new G4PVPlacement(0, G4ThreeVector(0, 0, -FR_thick / 2.0 -
                                                                              2 * (FR_thick + PEEK_Rod_thick)),
                                                       FR_logic, FR_logic->GetName(), gas_logic, 0, 0, false);
        G4VPhysicalVolume *FR_FC_2 = new G4PVPlacement(0, G4ThreeVector(0, 0, -FR_thick / 2.0 -
                                                                              3 * (FR_thick + PEEK_Rod_thick)),
                                                       FR_logic, FR_logic->GetName(), gas_logic, 0, 0, false);
        G4VPhysicalVolume *FR_FC_1 = new G4PVPlacement(0, G4ThreeVector(0, 0, -FR_thick / 2.0 -
                                                                              4 * (FR_thick + PEEK_Rod_thick)),
                                                       FR_logic, FR_logic->GetName(), gas_logic, 0, 0, false);


        // // EL Field Rings
        G4VPhysicalVolume *FR_EL_3 = new G4PVPlacement(0, G4ThreeVector(0, 0, -FR_thick / 2.0 -
                                                                              4 * (FR_thick + PEEK_Rod_thick) -
                                                                              2.5 * cm - 1.3 * cm - 0.7 * cm -
                                                                              1.3 * cm - 2 * cm - FR_thick -
                                                                              2 * (FR_thick + PEEK_Rod_thick)),
                                                       FR_logic, FR_logic->GetName(), gas_logic, 0, 0, false);
        G4VPhysicalVolume *FR_EL_2 = new G4PVPlacement(0, G4ThreeVector(0, 0, -FR_thick / 2.0 -
                                                                              4 * (FR_thick + PEEK_Rod_thick) -
                                                                              2.5 * cm - 1.3 * cm - 0.7 * cm -
                                                                              1.3 * cm - 2 * cm - FR_thick -
                                                                              1 * (FR_thick + PEEK_Rod_thick)),
                                                       FR_logic, FR_logic->GetName(), gas_logic, 0, 0, false);
        G4VPhysicalVolume *FR_EL_1 = new G4PVPlacement(0, G4ThreeVector(0, 0, -FR_thick / 2.0 -
                                                                              4 * (FR_thick + PEEK_Rod_thick) -
                                                                              2.5 * cm - 1.3 * cm - 0.7 * cm -
                                                                              1.3 * cm - 2 * cm - FR_thick), FR_logic,
                                                       FR_logic->GetName(), gas_logic, 0, 0, false);


        // PEEK Rods

        G4double x_rot_3 = 5 * std::sin(120 * deg) * cm;
        G4double y_rot_3 = -5 * std::cos(120 * deg) * cm;
        G4double x_rot_2 = 5 * std::sin(-120 * deg) * cm;
        G4double y_rot_2 = -5 * std::cos(-120 * deg) * cm;
        std::vector<G4double> x_rot_v = {0, x_rot_2, x_rot_3};
        std::vector<G4double> y_rot_v = {-5 * cm, y_rot_2, y_rot_3};
        G4double EL_PEEK_ROD_SHIFT = - 4 * (FR_thick + PEEK_Rod_thick) - 2.5 * cm - EL_thick - ElGap_ - EL_thick - 2 * cm - FR_thick;

        // Loop rotations
        for (int j = 0; j <=2; j++) {
            // The PEEK rods at the end
            new G4PVPlacement(0, G4ThreeVector(x_rot_v[j], y_rot_v[j], 1 * cm / 2.0 + 5 * (FR_thick + PEEK_Rod_thick)), PEEK_logic_cathode, PEEK_logic->GetName(), gas_logic, 0, 0, false);

            // The other FC PEEK rods
            for (int i = 4; i >= -4; i--) {
                new G4PVPlacement(0, G4ThreeVector(x_rot_v[j], y_rot_v[j], PEEK_Rod_thick / 2.0 + i * (FR_thick + PEEK_Rod_thick)), PEEK_logic, PEEK_logic->GetName(), gas_logic, 0, 0, false);
            }

            // PEEK EL
            new G4PVPlacement(0, G4ThreeVector(x_rot_v[j], y_rot_v[j], 1.1 * cm / 2.0       + EL_PEEK_ROD_SHIFT), PEEK_logic_buffer, PEEK_logic->GetName(), gas_logic, 0, 0, false);
            new G4PVPlacement(0, G4ThreeVector(x_rot_v[j], y_rot_v[j], PEEK_Rod_thick / 2.0 + EL_PEEK_ROD_SHIFT - 1 * (FR_thick + PEEK_Rod_thick)), PEEK_logic, PEEK_logic->GetName(), gas_logic, 0, 0, false);
            new G4PVPlacement(0, G4ThreeVector(x_rot_v[j], y_rot_v[j], PEEK_Rod_thick / 2.0 + EL_PEEK_ROD_SHIFT - 2 * (FR_thick + PEEK_Rod_thick)), PEEK_logic, PEEK_logic->GetName(), gas_logic, 0, 0, false);
            new G4PVPlacement(0, G4ThreeVector(x_rot_v[j], y_rot_v[j], 3.37 * cm / 2.0      + EL_PEEK_ROD_SHIFT - 2 * (FR_thick + PEEK_Rod_thick) - FR_thick - 3.37 * cm), PEEK_logic_buffer_end, PEEK_logic->GetName(), gas_logic, 0, 0, false);
        }


        // EL_Gap
        new G4PVPlacement(0, G4ThreeVector(0., 0., EL_pos), EL_logic, EL_solid->GetName(), gas_logic, 0, 0, false);

        G4VPhysicalVolume *EL_Ring_Plus = new G4PVPlacement(0, G4ThreeVector(0., 0., EL_thick / 2.0 - FR_thick -
                                                                                     4 * (FR_thick + PEEK_Rod_thick) -
                                                                                     2.5 * cm - EL_thick),
                                                            EL_ring_logic, EL_solid->GetName(), gas_logic, 0, 0, false);

        // Place the Mesh bits
        G4VPhysicalVolume *EL_Mesh_Plus_plus = new G4PVPlacement(rotateMesh, G4ThreeVector(0., 0.,
                                                                                           EL_thick / 2.0 - FR_thick -
                                                                                           4 *
                                                                                           (FR_thick + PEEK_Rod_thick) -
                                                                                           2.5 * cm - EL_thick -
                                                                                           EL_thick / 2.0),
                                                                 ELP_Disk_logic, ELP_Disk_logic->GetName(), gas_logic,
                                                                 0, 0, false);
        HexCreator->PlaceHexagons(nHole, EL_hex_size, EL_mesh_thick, ELP_Disk_logic, EL_Hex_logic);

        G4VPhysicalVolume *EL_Ring_Plus_plus = new G4PVPlacement(0, G4ThreeVector(0., 0., EL_thick / 2.0 - FR_thick -
                                                                                          4 *
                                                                                          (FR_thick + PEEK_Rod_thick) -
                                                                                          2.5 * cm - EL_thick - ElGap_ -
                                                                                          EL_thick), EL_ring_logic,
                                                                 EL_solid->GetName(), gas_logic, 0, 0, false);

        // Place the Mesh bits
        G4VPhysicalVolume *EL_Mesh_Plus = new G4PVPlacement(0, G4ThreeVector(0., 0., EL_thick / 2.0 - FR_thick -
                                                                                     4 * (FR_thick + PEEK_Rod_thick) -
                                                                                     2.5 * cm - EL_thick - ElGap_ -
                                                                                     EL_thick + EL_thick / 2.0),
                                                            ELPP_Disk_logic, ELPP_Disk_logic->GetName(), gas_logic, 0,
                                                            0, false);
        HexCreator->PlaceHexagons(nHole, EL_hex_size, EL_mesh_thick, ELPP_Disk_logic, EL_Hex_logic);


        // Cathode
        G4VPhysicalVolume *Cathode = new G4PVPlacement(0, G4ThreeVector(0., 0., EL_thick / 2.0 + 1 * cm +
                                                                                5 * (FR_thick + PEEK_Rod_thick)),
                                                       Cathode_ring_logic, EL_solid->GetName(), gas_logic, 0, 0, false);

        // Place the Mesh bits
        G4VPhysicalVolume *Cathode_EL_Mesh = new G4PVPlacement(rotateMesh, G4ThreeVector(0., 0.,
                                                                                         EL_thick / 2.0 + 1 * cm + 5 *
                                                                                                                   (FR_thick +
                                                                                                                    PEEK_Rod_thick) -
                                                                                         EL_thick / 2.0),
                                                               Cathode_Disk_logic, Cathode_Disk_logic->GetName(),
                                                               gas_logic, 0, 0, false);
        HexCreator->PlaceHexagons(nHole, EL_hex_size, EL_mesh_thick, Cathode_Disk_logic, EL_Hex_logic);

        // --- Optical ---
        new G4LogicalBorderSurface("SteelSurfaceFR1", gas_phys, FR_FC_1, OpSteelSurf);
        new G4LogicalBorderSurface("SteelSurfaceFR2", gas_phys, FR_FC_2, OpSteelSurf);
        new G4LogicalBorderSurface("SteelSurfaceFR3", gas_phys, FR_FC_3, OpSteelSurf);
        new G4LogicalBorderSurface("SteelSurfaceFR4", gas_phys, FR_FC_4, OpSteelSurf);
        new G4LogicalBorderSurface("SteelSurfaceFR5", gas_phys, FR_FC_5, OpSteelSurf);
        new G4LogicalBorderSurface("SteelSurfaceFR6", gas_phys, FR_FC_6, OpSteelSurf);
        new G4LogicalBorderSurface("SteelSurfaceFR7", gas_phys, FR_FC_7, OpSteelSurf);
        new G4LogicalBorderSurface("SteelSurfaceFR8", gas_phys, FR_FC_8, OpSteelSurf);
        new G4LogicalBorderSurface("SteelSurfaceFR9", gas_phys, FR_FC_9, OpSteelSurf);
        new G4LogicalBorderSurface("SteelSurfaceFR10", gas_phys, FR_FC_10, OpSteelSurf);

        new G4LogicalBorderSurface("SteelSurfaceFR_EL1", gas_phys, FR_EL_1, OpSteelSurf);
        new G4LogicalBorderSurface("SteelSurfaceFR_EL2", gas_phys, FR_EL_2, OpSteelSurf);
        new G4LogicalBorderSurface("SteelSurfaceFR_EL3", gas_phys, FR_EL_3, OpSteelSurf);

        new G4LogicalBorderSurface("SteelSurfaceELRing1", gas_phys, EL_Ring_Plus, OpSteelSurf);
        new G4LogicalBorderSurface("SteelSurfaceELRing2", gas_phys, EL_Ring_Plus_plus, OpSteelSurf);
        new G4LogicalBorderSurface("SteelSurfaceCathodeRing", gas_phys, Cathode, OpSteelSurf);


        new G4LogicalBorderSurface("SteelSurfaceELMesh", gas_phys, EL_Mesh_Plus_plus, OpSteelSurf);
        new G4LogicalBorderSurface("SteelSurfaceELMesh", gas_phys, EL_Mesh_Plus, OpSteelSurf);
        new G4LogicalBorderSurface("SteelSurfaceCathodeMesh", gas_phys, Cathode_EL_Mesh, OpSteelSurf);


        //  ----------------------- Needle Source -----------------------------
        // Source
        G4Tubs *SourceHolChamber_solid = new G4Tubs("SourceHolChamber", SourceEn_holedia / 2,
                                                    (SourceEn_diam / 2. + SourceEn_thickn),
                                                    (SourceEn_length / 2. + SourceEn_thickn), 0, twopi);
        G4LogicalVolume *SourceHolChamber_logic = new G4LogicalVolume(SourceHolChamber_solid, materials::Steel(),
                                                                      "SourceHolChamber_logic");

        G4Tubs *SourceHolChamberBlock_solid = new G4Tubs("SourceHolChBlock", 0, (SourceEn_holedia / 2),
                                                         (SourceEn_thickn / 2), 0., twopi);
        G4LogicalVolume *SourceHolChamberBlock_logic = new G4LogicalVolume(SourceHolChamberBlock_solid,
                                                                           materials::Steel(),
                                                                           "SourceHolChBlock_logic");

        /// Needle Source
        G4double NeedleyepRMin = 0;
        G4double NeedleyepRMax = (0.42) * mm;
        G4double NeedleyepDz = (2 / 2) * mm;
        G4double NeedleHalfLength = (2.56 / 2) * cm;
        G4double NeedleTailDiam = (0.6 / 2) * mm;
        G4double NeedleOffset = 1 * mm;

        G4Tubs *NeedleEye = new G4Tubs("NeedleEye", NeedleyepRMin, NeedleyepRMax, NeedleyepDz, 0., twopi);
        G4Tubs *NeedleTail = new G4Tubs("NeedleTail", NeedleyepRMin, NeedleTailDiam, NeedleHalfLength, 0., twopi);

        G4Tubs *Collimator = new G4Tubs("Collimator", (5.74 / 2) * mm, (11.5 / 2) * mm, 0.5 * cm, 0., twopi);
        G4Tubs *CollimatorBlock = new G4Tubs("CollimatorBlock", NeedleTailDiam, (5.74 / 2) * mm, 0.2 * cm, 0., twopi);

        // Combining them to create the Needle
        G4VSolid *Needle = new G4UnionSolid("Needle", NeedleEye, NeedleTail, 0, G4ThreeVector(0, 0, NeedleHalfLength));
        G4VSolid *CollimatorWithBlock = new G4UnionSolid("CollimatorWithBlock", Collimator, CollimatorBlock, 0,
                                                         G4ThreeVector(0, 0, 0.5 * cm - 0.2 * cm));

        G4LogicalVolume *Needle_Logic = new G4LogicalVolume(Needle, materials::Steel(), "Needle");
        G4LogicalVolume *Coll_Logic = new G4LogicalVolume(CollimatorWithBlock, materials::PEEK(),
                                                          "CollimatorWithBlock");


        // --- Placement ---

        // Source Holder
        G4VPhysicalVolume *Needle_Phys;
        if (!HideSourceHolder_) {

            // Needle Solid
            G4RotationMatrix *NeedleRotate = new G4RotationMatrix();
            NeedleRotate->rotateY(90. * deg);
            NeedleRotate->rotateX(90. * deg);
            NeedleRotate->rotateY(180. * deg);
            //NeedleRotate->rotateX(+10*deg);
            G4ThreeVector NeedlePos = {vtx_[0], vtx_[1] - NeedleOffset, vtx_[2] - FieldCagePos / 2};
            G4ThreeVector CollPosition = {NeedlePos[0], NeedlePos[1] - 5 * mm, NeedlePos[2]};

            Needle_Phys = new G4PVPlacement(NeedleRotate, NeedlePos, Needle_Logic, Needle->GetName(), gas_logic, true,
                                            0, false);

            // Collimator
            if (!HideCollimator_) {
                new G4PVPlacement(NeedleRotate, CollPosition, Coll_Logic, CollimatorWithBlock->GetName(), gas_logic,
                                  true, 0, false);
            }

        }

        // --- Optical ---
        if (!HideSourceHolder_ && !HideCollimator_) {
            new G4LogicalBorderSurface("SteelSurface_Needle", gas_phys, Needle_Phys, OpSteelSurf);
        }

        //  --------------------------- PMT -----------------------------------
        //Adding the PMTs in here
        // Krishan: Nexus complaining about PMTs being created twice, so commented out
        pmt1_ = new PmtR7378A();
        // pmt2_ = new PmtR7378A();
        // pmt1_->SetPMTName("S2");
        // pmt2_->SetPMTName("S1");
        pmt1_->Construct();
        // pmt2_->Construct();


        // Adding Logical Volumes for PMTs
        G4LogicalVolume *pmt1_logic = pmt1_->GetLogicalVolume();
        G4LogicalVolume *pmt2_logic = pmt1_->GetLogicalVolume();


        // PMT1 and PMT3

        // // PMTs
        G4double PMT_offset = 0.2 * cm;
        G4double PMT_pos =
                (chamber_length / 2) + chamber_thickn + (pmt1_->Length() / 2) + MgF2_window_thickness_ + PMT_offset;
        G4RotationMatrix *pmt1rotate = new G4RotationMatrix();
        pmt1rotate->rotateX(180. * deg);


        //// PMT Covering Tube ///
        G4double offset = 1.65 * cm;
        G4double PMT_Tube_Length1 =
                MgF2_window_thickness_ + (pmt1_->Length() + 0.5 * cm) / 2 + offset - PMT_offset - 0.05 * cm;
        G4double PMT_Tube_Length0 = (17 * cm + pmt1_->Length()) / 2 - 3.87 * cm - PMT_offset;
        G4double PMT_Tube_Block_Thickness = 0.2 * cm;
        G4double LongPMTTubeOffset = 7.5 * cm - 3.9 * cm;
        G4double PMTTubeDiam = 2.54 * cm;

        // Tube Away from EL
        G4Tubs *PMT_Tube_solid0 = new G4Tubs("PMT_TUBE0", (PMTTubeDiam / 2) + 0.5 * cm, (PMTTubeDiam / 2) + 0.7 * cm,
                                             PMT_Tube_Length0, 0, twopi);
        G4LogicalVolume *PMT_Tube_Logic0 = new G4LogicalVolume(PMT_Tube_solid0, materials::Steel(),
                                                               PMT_Tube_solid0->GetName());

        G4Tubs *PMT_Block_solid0 = new G4Tubs("PMT_TUBE_BLOCK0", 0, (PMTTubeDiam / 2 + 0.5 * cm),
                                              PMT_Tube_Block_Thickness, 0, twopi);
        G4LogicalVolume *PMT_Block_Logic0 = new G4LogicalVolume(PMT_Block_solid0, materials::Steel(),
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
        G4LogicalVolume *PMT_Tube_Logic1 = new G4LogicalVolume(PMT_Tube_solid1, materials::Steel(),
                                                               PMT_Tube_solid1->GetName());
        G4Tubs *PMT_Block_solid1 = new G4Tubs("PMT_TUBE_BLOCK1", 0, (PMTTubeDiam / 2 + 0.5 * cm),
                                              PMT_Tube_Block_Thickness, 0, twopi);
        G4LogicalVolume *PMT_Block_Logic = new G4LogicalVolume(PMT_Block_solid1, materials::Steel(),
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
                                                              lab_logic_volume, false, 0, false);
        G4VPhysicalVolume *PMT_Tube_Phys1 = new G4PVPlacement(0, G4ThreeVector(0, 0, -(PMT_pos - PMT_offset) - offset),
                                                              PMT_Tube_Logic1, PMT_Tube_Logic1->GetName(),
                                                              lab_logic_volume, false, 0, false);

        // PMT Tube Vacuum
        G4VPhysicalVolume *PMT_Tube_Vacuum_Phys0 = new G4PVPlacement(0,
                                                                     G4ThreeVector(0, 0, PMT_pos + LongPMTTubeOffset),
                                                                     InsideThePMT_Tube_Logic0, "PMT_TUBE_VACUUM0",
                                                                     lab_logic_volume, false, 0, false);
        G4VPhysicalVolume *PMT_Tube_Vacuum_Phys1 = new G4PVPlacement(0, G4ThreeVector(0, 0,
                                                                                      -(PMT_pos - PMT_offset) - offset),
                                                                     InsideThePMT_Tube_Logic1, "PMT_TUBE_VACUUM1",
                                                                     lab_logic_volume, false, 0, false);

        // PMT Tube Block
        new G4PVPlacement(0, G4ThreeVector(0, 0,
                                           PMT_pos - PMT_offset + PMT_Tube_Length0 - PMT_Tube_Block_Thickness / 2 +
                                           LongPMTTubeOffset), PMT_Block_Logic0, PMT_Block_Logic0->GetName(),
                          lab_logic_volume, false, 0, false);
        new G4PVPlacement(0, G4ThreeVector(0, 0,
                                           -(PMT_pos - PMT_offset + PMT_Tube_Length1 - PMT_Tube_Block_Thickness / 2) -
                                           offset), PMT_Block_Logic, PMT_Block_Logic->GetName(), lab_logic_volume,
                          false, 1, false);

        // PMTs
        // new G4PVPlacement(pmt1rotate,G4ThreeVector (0,0,((PMT3_Pos_)-pmt1_->Length()/2-PMT_Tube_Length1/2-MgF2_window_thickness_/2)),pmt1_logic,pmt1_->GetPMTName(),InsideThePMT_Tube_Logic0,true,0,false);
        new G4PVPlacement(0, G4ThreeVector(0, 0., (PMT1_Pos_ - pmt1_->Length() / 2 - MgF2_window_thickness_ / 2)),
                          pmt2_logic, "S1", InsideThePMT_Tube_Logic1, true, 0, false);


        // --- Optical ---        
        new G4LogicalBorderSurface("SteelSurface_PMT3_Enclosing", PMT_Tube_Vacuum_Phys0, PMT_Tube_Phys0, OpSteelSurf);
        new G4LogicalBorderSurface("SteelSurface_PMT1_Enclosing", PMT_Tube_Vacuum_Phys1, PMT_Tube_Phys1, OpSteelSurf);


        //  ------------------------ Camera -----------------------------------
        
        // Window
        G4double camHalfLength = 0.5 * mm;
        G4double camRadius = 12.7 * mm;
        G4VSolid *camSolid = new G4Tubs("camWindow", 0., camRadius, camHalfLength, 0., twopi);
        G4LogicalVolume *camLogical = new G4LogicalVolume(camSolid, MgF2, "camLogical");


        // --- Placement ---
        G4double ImageDist = 7.945 * cm; // Got from trial and error
        G4VPhysicalVolume *camPhysical = new G4PVPlacement(0, G4ThreeVector(0, 0, (chamber_length / 2 + chamber_thickn +
                                                                                   ImageDist) - PMT_pos -
                                                                                  LongPMTTubeOffset), camLogical,
                                                           "camPhysical", InsideThePMT_Tube_Logic0, false, 0, false);

        // Camera
        G4OpticalSurface *opXenon_Glass = new G4OpticalSurface("XenonCamSurface");
        opXenon_Glass->SetModel(glisur);                  // SetModel
        opXenon_Glass->SetType(dielectric_dielectric);   // SetType
        opXenon_Glass->SetFinish(ground);                 // SetFinish
        opXenon_Glass->SetPolish(0.0);
        new G4LogicalBorderSurface("XenonCamSurface", PMT_Tube_Vacuum_Phys0, camPhysical, opXenon_Glass);

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
                                                                bracket_logical, "bracketPhysical", gas_logic, false, 0,
                                                                false);
        G4VPhysicalVolume *bracketPhysical2 = new G4PVPlacement(rotateZ_120, G4ThreeVector(x_rot_2, y_rot_2, EL_pos),
                                                                bracket_logical, "bracketPhysical", gas_logic, false, 0,
                                                                false);
        G4VPhysicalVolume *bracketPhysical3 = new G4PVPlacement(rotateZ_m120, G4ThreeVector(x_rot_3, y_rot_3, EL_pos),
                                                                bracket_logical, "bracketPhysical", gas_logic, false, 0,
                                                                false);


        // ____________________________________________________________________
        // ================= Detector Properties  =============================

        //Construct a G4Region, connected to the logical volume in which you want to use the G4FastSimulationModel
        G4Region *regionGas = new G4Region("GasRegion");
        regionGas->AddRootLogicalVolume(gas_logic);


        G4SDManager *SDManager = G4SDManager::GetSDMpointer();
        IonizationSD* ionisd = new IonizationSD("/CRAB0/ACTIVE");
        SDManager->SetVerboseLevel(1);
        SDManager->AddNewDetector(ionisd);
        gas_logic->SetSensitiveDetector(ionisd);

        GarfieldHelper GH(chamber_diam/2.0/cm, chamber_length/cm, Active_diam/2.0/cm , FielCageGap/cm, gas_pressure_, ElGap_, fieldDrift_, fieldEL_);

        //These commands generate the four gas models and connect it to the GasRegion
        G4Region *region = G4RegionStore::GetInstance()->GetRegion("GasRegion");
        new DegradModel("DegradModel",region, GH);
        new GarfieldVUVPhotonModel("GarfieldVUVPhotonModel",region, GH);

        // Visuals
        AssignVisuals();

        return;

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
        G4VisAttributes *GasVa = new G4VisAttributes(nexus::WhiteAlpha());
        GasVa->SetForceCloud(true);
        Gas->SetVisAttributes(GasVa);

        //Source Enclosure Related
        G4LogicalVolume *SourceHolder = lvStore->GetVolume("SourceHolChamber_logic");
        G4LogicalVolume *Needle = lvStore->GetVolume("Needle");
        G4LogicalVolume *Collimator = lvStore->GetVolume("CollimatorWithBlock");
        G4VisAttributes *CollimatorVa = new G4VisAttributes(nexus::YellowAlpha());
        CollimatorVa->SetForceSolid(true);

        Collimator->SetVisAttributes(CollimatorVa);


        Needle->SetVisAttributes(ChamberVa);

        G4LogicalVolume *SourceHolderBlock = lvStore->GetVolume("SourceHolChBlock_logic");
        G4VisAttributes *SourceHolderVa = new G4VisAttributes(G4Colour(2, 2, 2));
        SourceHolderVa->SetForceSolid(true);

        // Flange
        G4LogicalVolume *flangeLog = lvStore->GetVolume("CHAMBER_FLANGE");
        G4VisAttributes flangeVis = nexus::DarkGreyAlpha();
        flangeVis.SetForceSolid(true);
        flangeLog->SetVisAttributes(ChamberVa);


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

        // Brackets
        G4LogicalVolume *BracketLog = lvStore->GetVolume("bracketLogical");
        G4VisAttributes BracketVis = nexus::DirtyWhiteAlpha();
        BracketVis.SetForceSolid(true);
        BracketLog->SetVisAttributes(BracketVis);


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
        G4VisAttributes PmttubeVacuumVis = nexus::DarkGreyAlpha();
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
        CAMLog->SetVisAttributes(CAMLog);

        // EL-Region
        G4LogicalVolume *ELLogic = lvStore->GetVolume("EL_GAP");
        G4VisAttributes ELVis = nexus::BlueAlpha();
        ELVis.SetForceCloud(true);
        ELLogic->SetVisAttributes(G4VisAttributes::GetInvisible());

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

        // FieldCage
        G4LogicalVolume *FieldCage = lvStore->GetVolume("FIELDCAGE");
        G4VisAttributes FielCageVis = nexus::Red();
        FielCageVis.SetForceCloud(true);
        FieldCage->SetVisAttributes(G4VisAttributes::GetInvisible());


        SourceHolder->SetVisAttributes(SourceHolderVa);
        SourceHolderBlock->SetVisAttributes(SourceHolderVa);
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

    else {
        G4Exception("[CRAB0]", "GenerateVertex()", FatalException,
        "Unknown vertex generation region!");
    }

    return vertex;
    }




}

