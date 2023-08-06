//
//  PrimaryGeneratorAction.hh
//  Xenon
//
//  Created by Lennert De Keukeleere on 25/10/2018.
//

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4Navigator.hh"
#include "DetectorConstruction.hh"
#include "G4ParticleTable.hh"
#include "G4GeneralParticleSource.hh"
#include "PrimaryGenerator.hh"

class G4Event;


class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
    PrimaryGeneratorAction();
    ~PrimaryGeneratorAction();
    
    // you must define this method, it is called by the G4RunManager
    // run manager passes the pointer to an event object, it will be given
    // attributes from the Particle Gun

    G4GeneralParticleSource* GetParticleGun() {return fparticleGun;};
  //    PrimaryGenerator* GetPrimaryGenerator(){return fPrimaryGenerator;};

    void GeneratePrimaries(G4Event*);
    
private:
    
    G4GeneralParticleSource* fparticleGun;
    PrimaryGenerator* fPrimaryGenerator;
    
    
};

#endif /* PrimaryGeneratorAction_h */
