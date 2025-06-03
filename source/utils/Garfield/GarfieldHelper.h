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
    struct COMSOL_Variables{
        G4String Path,MeshFile,Data,Materialstxt;
        G4bool useCOMSOL;
        G4bool useOlderSimple;
        G4ThreeVector gGainReduction;
        G4ThreeVector gGainSdev;
        G4double ELYield;
        G4float stepsize;
    };
    /// Default constructor 
    GarfieldHelper();
    GarfieldHelper(G4double DetChamberR ,G4double DetChamberL_,G4double DetActiveR ,G4double DetActiveL ,G4double GasPressure ,G4double gap_EL ,G4double fieldDrift ,G4double fieldEL);


    void DumpParams();
    void SetGasFile(G4String g);
    void SetCOMSOLVariables(COMSOL_Variables *ComsolVariable );
    COMSOL_Variables * GetComsolVariables_();
    G4String GetGasFile();
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
    G4String GasFile_;

    G4double thermalE_{1.3*eV}; // eV

    COMSOL_Variables * ComsolVariable_;


  };

  // inline methods ..................................................

  inline GarfieldHelper::~GarfieldHelper() { }
  inline void GarfieldHelper::SetGasFile(G4String g) {GasFile_=g;}
  inline G4String GarfieldHelper::GetGasFile() { return GasFile_;}
  inline void GarfieldHelper::SetCOMSOLVariables(nexus::GarfieldHelper::COMSOL_Variables *ComsolVariable) {
      ComsolVariable_=ComsolVariable;
  }
  //GarfieldHelper::COMSOL_Variables * GarfieldHelper::GetComsolVariables_() { return ComsolVariable_;}


} // namespace nexus

#endif
