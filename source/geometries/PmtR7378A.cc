// ----------------------------------------------------------------------------
// nexus | PmtR7378A.cc
//
// Geometry of the Hamamatsu R7378A photomultiplier.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "PmtR7378A.h"
#include "FileHandling.h"
#include "SensorSD.h"
#include "OpticalMaterialProperties.h"
#include "MaterialsList.h"
#include "Visibilities.h"

#include <G4NistManager.hh>
#include <G4Tubs.hh>
#include <G4Sphere.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4Material.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4LogicalBorderSurface.hh>
#include <G4VisAttributes.hh>
#include <G4SDManager.hh>

#include <CLHEP/Units/SystemOfUnits.h>

namespace nexus {
  using namespace CLHEP;

  PmtR7378A::PmtR7378A(): GeometryBase()
  {

  }



  PmtR7378A::~PmtR7378A()
  {
  }



  void PmtR7378A::Construct()
  {
    // PMT BODY //////////////////////////////////////////////////////

    pmt_diam_   = 25.4 * mm;
    pmt_length_ = 43.0 * mm;
    PmtName= GetPMTName();
    G4Tubs* pmt_solid =
      new G4Tubs("PMT_R7378A", 0., pmt_diam_/2., pmt_length_/2., 0., twopi);

    G4Material* aluminum =
      G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");

    G4LogicalVolume* pmt_logic =
      new G4LogicalVolume(pmt_solid, aluminum, "PMT_R7378A");

    G4VisAttributes vis_body = nexus::TitaniumGreyAlpha();;
    vis_body.SetForceSolid(true);
    pmt_logic->SetVisAttributes(vis_body);

    this->SetLogicalVolume(pmt_logic);


    // PMT WINDOW ////////////////////////////////////////////////////

    G4double window_diam = pmt_diam_;
    G4double window_length = 6. * mm;

    G4Tubs* window_solid =
      new G4Tubs("PMT_WINDOW", 0., window_diam/2., window_length/2., 0., twopi);

    G4Material* quartz =
      G4NistManager::Instance()->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
    quartz->SetMaterialPropertiesTable(opticalprops::FusedSilica());
    G4String WindowName;
    G4String DetName;
      G4String SensDet;
    if(!PmtName){
        WindowName="PMT_WINDOW";
        DetName="PHOTOCATHODE";
        SensDet="/PMT_R7378A/Pmt";
    }else{
        WindowName=PmtName+"_WINDOW";
        DetName=PmtName+"_PHOTOCATHODE";
        SensDet="/PMT_R7378A/"+PmtName;

    }
    G4LogicalVolume* window_logic =
      new G4LogicalVolume(window_solid, quartz,WindowName );

    new G4PVPlacement(0, G4ThreeVector(0., 0., (pmt_length_-window_length)/2.),
		      window_logic,WindowName, pmt_logic, false, 0, false);

    G4VisAttributes wndw_col = nexus::Blue();
    window_logic->SetVisAttributes(wndw_col);


    // PHOTOCATHODE //////////////////////////////////////////////////

    G4double phcath_diam   = 22.0 * mm;
    G4double phcath_height =  4.0 * mm;
    G4double phcath_thickn =  0.1 * mm;
    G4double phcath_posz   =-16.0 * mm;

    G4double rmax =
      0.5 * (phcath_diam*phcath_diam / (4*phcath_height) + phcath_height);
    G4double rmin = rmax - phcath_thickn;
    G4double theta = asin(phcath_diam/(2.*rmax));

    G4Sphere* phcath_solid =
      new G4Sphere(DetName, rmin, rmax, 0, twopi, 0, theta);

    G4LogicalVolume* phcath_logic =
      new G4LogicalVolume(phcath_solid, aluminum, DetName);

    G4VisAttributes vis_solid;
    vis_solid.SetForceSolid(true);
    phcath_logic->SetVisAttributes(vis_solid);

    //G4PVPlacement* phcath_physi =
      new G4PVPlacement(0, G4ThreeVector(0.,0.,phcath_posz), phcath_logic,
			DetName, window_logic, false, 0, false);

    // Sensitive detector
    SensorSD* pmtsd = new SensorSD(SensDet);
    pmtsd->SetDetectorVolumeDepth(2);
    pmtsd->SetTimeBinning(100.*nanosecond);
    G4SDManager::GetSDMpointer()->AddNewDetector(pmtsd);
    phcath_logic->SetSensitiveDetector(pmtsd);

    // OPTICAL SURFACES //////////////////////////////////////////////

    // The values for the efficiency are chosen in order to match
    // the curve of the quantum efficiency provided by Hamamatsu:
    // http://sales.hamamatsu.com/en/products/electron-tube-division/detectors/photomultiplier-tubes/part-r7378a.php
    // The source of light is point-like, isotropic and it has been placed at a
    // distance of 25 cm from the surface of the PMT window.
    // The quantity to be compared with the Hamamatsu curve is:
    // number of detected photons/ number of photons that reach the PMT window.
    // The total number of generated photons is taken as
    // number of photons that reach the PMT window because
    // light is generated only in that fraction of solid angle that subtends the
    // window of the PMT.
    G4MaterialPropertiesTable* phcath_mpt = new G4MaterialPropertiesTable();
    FileHandling *f1=new FileHandling();
    std::vector<std::vector<G4double>> t;
    /// Note: In order this efficiency to work , you need to list it from low energy to high energy
    t=f1->GetData("data/PMTR7378A_Efficiency.txt",',',1);
    //t=f1->GetData("data/test.txt",',',1);
    const G4int entries = t[0].size();
    G4double ENERGIES[entries];
    G4double EFFICIENCY[entries];
    G4double REFLECTIVITY[entries];

  for(int kk=0;kk<t[0].size();kk++){

        ENERGIES[kk]=t[0][kk]*eV;
        EFFICIENCY[kk]=t[1][kk];
        REFLECTIVITY[kk]=0.;
        //G4cout<<ENERGIES[kk]<<" "<<EFFICIENCY[kk]<<G4endl;
    }

    /*
      const G4int entries = 31;
      G4double ENERGIES[entries] =
      {1.72194*eV, 1.77114*eV, 1.82324*eV, 1.87848*eV, 1.93719*eV,
       1.99968*eV,  2.06633*eV, 2.13759*eV, 2.21393*eV, 2.29593*eV,
       2.38423*eV, 2.47960*eV, 2.58292*eV, 2.69522*eV, 2.81773*eV,
       2.95190*eV, 3.0995*eV, 3.26263*eV, 3.44389*eV, 3.64647*eV,
       3.87438*eV, 4.13267*eV, 4.42786*eV, 4.76846*eV, 5.16583*eV,
       5.63545*eV, 6.19900*eV, 6.88778*eV,7.08*eV, 7.74875*eV, 8.85571*eV};


       G4double EFFICIENCY[entries] =
      { 0.00000, 0.00028, 0.00100, 0.00500, 0.00100,
    	0.02200, 0.04500, 0.07000, 0.11500, 0.16000,
    	0.20500, 0.23500, 0.27000, 0.29000, 0.31300,
    	0.35200, 0.38000, 0.38000, 0.37300, 0.37300,
    	0.37000, 0.36000, 0.35500, 0.33500, 0.31000,
    	0.29500, 0.27500, 0.23000,1, 0.52000, 0.00000};

       G4double REFLECTIVITY[entries] =
      { 0., 0., 0., 0., 0.,
	0., 0., 0., 0., 0.,
	0., 0., 0., 0., 0.,
	0., 0., 0., 0., 0.,
	0., 0., 0., 0., 0.,
	0., 0., 0., 0., 0.,0 };

    */


   phcath_mpt->AddProperty("EFFICIENCY", ENERGIES, EFFICIENCY, entries);
   phcath_mpt->AddProperty("REFLECTIVITY", ENERGIES, REFLECTIVITY, entries);


    G4OpticalSurface* phcath_opsur =new G4OpticalSurface(DetName, unified, polished, dielectric_metal);
    phcath_opsur->SetMaterialPropertiesTable(phcath_mpt);

    new G4LogicalSkinSurface(DetName, phcath_logic, phcath_opsur);

  }


} // end namespace nexus
