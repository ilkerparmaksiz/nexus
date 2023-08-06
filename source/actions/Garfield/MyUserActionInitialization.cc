#include "MyUserActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "GarfieldRunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "GarfieldEventAction.hh"
#include "G4SDManager.hh"
#include "GasBoxSD.hh"
#include "GarfieldSteppingAction.hh"
#include "GarfieldTrackingAction.hh"

MyUserActionInitialization::MyUserActionInitialization(){}

MyUserActionInitialization::~MyUserActionInitialization(){}

void MyUserActionInitialization::Build() const {

        PrimaryGeneratorAction* primary = new PrimaryGeneratorAction();
	SetUserAction(primary);

	EventAction* evt = new EventAction();
	SetUserAction(evt);

	SteppingAction* stepAct = new SteppingAction(evt);
	SetUserAction(stepAct);

	SetUserAction(new RunAction());

	SetUserAction(new NESTStackingAction()); // comment to launch, e.g., opticalphotons as primaries. EC, 29-July-2022.
	TrackingAction* trackAct = new TrackingAction();
	SetUserAction(trackAct);
	
}

void MyUserActionInitialization::BuildForMaster() const {
	SetUserAction(new RunAction());
}
