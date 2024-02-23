// ----------------------------------------------------------------------------
// nexus | GarfieldHelper.h
//
// This class provides helpers for running garfield
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#ifndef GARFIELD_HELPER_H
#define GARFIELD_HELPER_H

#include <CLHEP/Units/SystemOfUnits.h>
#include "G4UnitsTable.hh"

using namespace CLHEP;

namespace nexus {

  class GarfieldHelper
  {
  public:

    /// Default constructor 
    GarfieldHelper();
    GarfieldHelper(G4double DetChamberR, G4double DetChamberL, G4double DetActiveR,
                   G4double DetActiveL, G4double GasPressure, G4double gap_EL,
                   G4double fieldDrift, G4double fieldEL, G4double v_drift,
                   G4double v_drift_el, G4double e_liftetime, G4bool useCOMSOL, 
                   G4bool useELFile, G4bool useDEGRAD, G4String DetName);

    void DumpParams();

    /// Destructor
    ~GarfieldHelper();

    // Detector geometry
    G4double DetChamberR_; // cm
    G4double DetChamberL_; // cm
    G4double DetActiveR_;  // cm
    G4double DetActiveL_;  // cm
    G4double GasPressure_; // bar

    G4double gap_EL_;      // cm
    G4double fieldDrift_;  // V/cm
    G4double fieldEL_;     // V/cm 
    G4double v_drift_;     // mm/us
    G4double v_drift_el_;  // mm/us
    G4double e_lifetime_;  // ms
    G4bool useCOMSOL_;
    G4bool useELFile_;
    G4bool useDEGRAD_;
    G4String DetName_; // Name of the detector

    G4double thermalE_{1.3*eV}; // eV


  };

  // inline methods ..................................................

  inline GarfieldHelper::~GarfieldHelper() { }

} // namespace nexus

#endif
