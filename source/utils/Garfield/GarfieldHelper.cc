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
                                 G4double DetActiveL, G4double CathodePos, G4double ELPos, G4int nsides,
                                 G4double GasPressure, G4double gap_EL, G4ThreeVector origin,
                                 G4double fieldDrift, G4double fieldEL, G4double v_drift,
                                 G4double v_drift_el, G4double e_lifetime, G4bool useCOMSOL, 
                                 G4bool useELFile,  G4bool useDEGRAD, G4String DetName){
    DetChamberR_ = DetChamberR;
    DetChamberL_ = DetChamberL;
    DetActiveR_  = DetActiveR;
    DetActiveL_  = DetActiveL;
    CathodePos_  = CathodePos;
    ELPos_       = ELPos;
    nsides_      = nsides;
    GasPressure_ = GasPressure;
    gap_EL_      = gap_EL;
    origin_      = origin;
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
    "Chamber Radius: "        << DetChamberR_ << " cm"    << "\n" <<
    "Chamber Length: "        << DetChamberL_ << " cm"    << "\n" <<
    "Active Radius: "         << DetActiveR_  << " cm"    << "\n" <<
    "Chamber length: "        << DetActiveL_  << " cm"    << "\n" <<
    "Cathode Position: "      << CathodePos_  << " cm"    << "\n" <<
    "EL Position: "           << ELPos_       << " cm"    << "\n" <<
    "N polygon sides: "       << nsides_      << "   "    << "\n" <<
    "Gas Pressure: "          << GasPressure_ << " Pa"    << "\n" <<
    "EL Gap: "                << gap_EL_      << " cm"    << "\n" <<
    "OriginXYZ: "             << "( " << origin_.x()      << ", " << 
                              origin_.y()     << ", "     << origin_.z() 
                              << ") cm" <<"\n"<<
    "Drift Field: "           << fieldDrift_  << " V/cm"  << "\n" <<
    "EL Field: "              << fieldEL_     << " kV/cm" << "\n" <<
    "Drift Velocity: "        << v_drift_     << " mm/ns" << "\n" <<
    "EL Drift Velocity: "     << v_drift_el_  << " mm/ns" << "\n" <<
    "Electron Lifetime: "     << e_lifetime_  << " ns"    << "\n" <<
    "COMSOL: "                << useCOMSOL_   << "  "     << "\n" <<
    "Generate from EL file: " << useELFile_   << "  "     << "\n" <<
    "DEGRAD: "                << useDEGRAD_   << "  "     << "\n" <<
    std::endl;

  }


} // end namespace nexus
