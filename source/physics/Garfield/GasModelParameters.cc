#include "GasModelParameters.h"
#include "DegradModel.h"
#include "GasModelParametersMessenger.h"
#include "DetectorConstruction.h"

using namespace nexus;

GasModelParameters::GasModelParameters(){
        fMessenger = new GasModelParametersMessenger(this);
}
GasModelParameters::~GasModelParameters{};
