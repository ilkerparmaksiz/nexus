#include "GasModelParametersMessenger.h"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UIparameter.hh"
#include "GasModelParameters.h"
#include "DegradModel.h"

#include "G4Tokenizer.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
using namespace nexus;
GasModelParametersMessenger::GasModelParametersMessenger(GasModelParameters* gm)
    : fGasModelParameters(gm) {
  GasModelParametersDir = new G4UIdirectory("/gasModelParameters/");
  GasModelParametersDir->SetGuidance("GasModelParameters specific controls");
  DegradDir = new G4UIdirectory("/gasModelParameters/degrad/");
  DegradDir->SetGuidance("Degrad specific controls");

  thermalEnergyCmd = new G4UIcmdWithADoubleAndUnit("/gasModelParameters/degrad/thermalenergy",this);
  thermalEnergyCmd->SetGuidance("Set the thermal energy to be used by degrad");

  GeomDir = new G4UIdirectory("/gasModelParameters/degrad/");
  setComsolCmd = new G4UIcmdWithABool("/gasModelParameters/geometry/useComsol", this);
  setComsolCmd->SetDefaultValue(false);

  setEL_FileCmd = new G4UIcmdWithABool("/gasModelParameters/geometry/useEL_File", this);
  setEL_FileCmd->SetDefaultValue(false);

  setCOMSOL_Path = new G4UIcmdWithAString("/gasModelParameters/geometry/COMSOL_Path", this);
  
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GasModelParametersMessenger::~GasModelParametersMessenger() {
  delete GasModelParametersDir;
  delete DegradDir;
  delete GeomDir;
  delete thermalEnergyCmd;
  delete setComsolCmd;
  delete setEL_FileCmd;
  delete setCOMSOL_Path;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GasModelParametersMessenger::SetNewValue(G4UIcommand* command, G4String newValues) {
    if(command == thermalEnergyCmd){
      fGasModelParameters->SetThermalEnergy(thermalEnergyCmd->GetNewDoubleValue(newValues));
    }

    if (command == setComsolCmd)
      fGasModelParameters->SetComsol(setComsolCmd->GetNewBoolValue(newValues));

    if (command == setEL_FileCmd)
      fGasModelParameters->SetEL_file(setEL_FileCmd->GetNewBoolValue(newValues));

    if (command == setCOMSOL_Path)
      fGasModelParameters->SetCOMSOL_Path(newValues);

}
