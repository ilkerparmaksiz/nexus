#include "GasModelParameters.hh"
#include "DegradModel.hh"
#include "GasModelParametersMessenger.hh"
#include "DetectorConstruction.hh"
namespace nexus {
    GasModelParameters::GasModelParameters(){
        fMessenger = new GasModelParametersMessenger(this);
}
}