#ifndef MyUserActionInitialization_hh
#define MyUserActionInitialization_hh

#include "G4VUserActionInitialization.hh"

#include "G4/NESTStackingAction.hh"

class DetectorConstruction;
class PhysicsList;

class MyUserActionInitialization : public G4VUserActionInitialization{
	public:
	MyUserActionInitialization();
	~MyUserActionInitialization();
	
	void Build() const;
	void BuildForMaster() const;
	
	private:
	
};

#endif
