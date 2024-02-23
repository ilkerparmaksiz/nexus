// ----------------------------------------------------------------------------
// nexus | GarfieldHelper.cc
//
// This class provides tools for running Garfield
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "GarfieldHelper.h"

namespace nexus {


  GarfieldHelper::GarfieldHelper(G4double DetChamberR, G4double DetChamberL, G4double DetActiveR,
                                 G4double DetActiveL, G4double GasPressure, G4double gap_EL,
                                 G4double fieldDrift, G4double fieldEL, G4double v_drift,
                                 G4double v_drift_el, G4double e_lifetime, G4bool useCOMSOL, 
                                 G4bool useELFile,  G4bool useDEGRAD, G4String DetName){
    DetChamberR_ = DetChamberR;
    DetChamberL_ = DetChamberL;
    DetActiveR_  = DetActiveR;
    DetActiveL_  = DetActiveL;
    GasPressure_ = GasPressure;
    gap_EL_      = gap_EL;
    fieldDrift_  = fieldDrift;
    fieldEL_     = fieldEL;
    v_drift_     = v_drift;
    v_drift_el_  = v_drift_el;
    e_lifetime_  = e_lifetime;
    useCOMSOL_   = useCOMSOL;
    useELFile_   = useELFile;
    useDEGRAD_   = useDEGRAD;
    DetName_     = DetName;

  }

  GarfieldHelper::GarfieldHelper(){};


  void GarfieldHelper::DumpParams(){

    std::cout << 
    "\n\nPrinting Garfield geometry params\n"
    "Chamber Radius: " << DetChamberR_ << "\n" <<
    "Chamber Length: " << DetChamberL_ << "\n" <<
    "Active Radius: " << DetActiveR_ << "\n" <<
    "Chamber length: " << DetActiveL_ << "\n" <<
    "Gas Pressure: " << GasPressure_ << "\n" <<
    "EL Gap: " << gap_EL_ << "\n" <<
    "Drift Field: " << fieldDrift_ << "\n" <<
    "EL Field: " << fieldEL_ << "\n" <<
    "Drift Velocity: " << v_drift_ << "\n" <<
    "EL Drift Velocity: " << v_drift_el_ << "\n" <<
    "Electron Lifetime: " << e_lifetime_ << "\n" <<
    "COMSOL: " << useCOMSOL_ << "\n" <<
    "Generate from EL file: " << useELFile_ << "\n" <<
    "DEGRAD: " << useDEGRAD_ << "\n" <<
    std::endl;

  }


} // end namespace nexus
