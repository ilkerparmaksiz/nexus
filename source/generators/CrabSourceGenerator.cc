//
// Created by ilker on 9/7/21.
//

#include "CrabSourceGenerator.h"
#include "GeometryBase.h"
#include "DetectorConstruction.h"
#include "FactoryBase.h"

#include <G4GenericMessenger.hh>
#include <G4RunManager.hh>
#include <G4ParticleDefinition.hh>
#include <G4PrimaryVertex.hh>
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"
#include "G4LogicalVolumeStore.hh"

#include "G4ParticleDefinition.hh"
#include "G4GenericMessenger.hh"


#include <math.h>

using namespace nexus;

REGISTER_CLASS(CrabSourceGenerator, G4VPrimaryGenerator)

//Generate
//////////////////////////////////////////////////////
CrabSourceGenerator::CrabSourceGenerator():
        G4VPrimaryGenerator(),
        atomic_number_(0),
        mass_number_(0),
        energy_level_(0.),
        decay_at_time_zero_(true),
        event_window_(0),
        N_Decays_per_s_(0.),
        region_(""),
        msg_(nullptr), geom_(nullptr),NDecays(1)
{
    msg_ = new G4GenericMessenger(this, "/Generator/CrabSourceGenerator/",

                                  "Control commands of the ion gun primary generator.");


    msg_->DeclareProperty("event_window", event_window_,"Event Window of the Decay");


    G4GenericMessenger::Command& atomic_number_cmd =
            msg_->DeclareProperty("atomic_number", atomic_number_, "Atomic number of the ion.");
    atomic_number_cmd.SetParameterName("atomic_number", false);
    atomic_number_cmd.SetRange("atomic_number > 0");

    G4GenericMessenger::Command& mass_number_cmd =
            msg_->DeclareProperty("mass_number", mass_number_, "Mass number of the ion.");
    mass_number_cmd.SetParameterName("mass_number", false);
    mass_number_cmd.SetRange("mass_number > 0");

    G4GenericMessenger::Command& energy_level_cmd =
            msg_->DeclarePropertyWithUnit("energy_level", "MeV", energy_level_,
                                          "Energy level of the ion.");
    energy_level_cmd.SetParameterName("energy_level", false);
    energy_level_cmd.SetRange("energy_level >= 0");

    msg_->DeclareProperty("decay_at_time_zero", decay_at_time_zero_,
                          "Set to true to make unstable ions decay at t=0.");

    msg_->DeclareProperty("region", region_,
                          "Region of the geometry where vertices will be generated.");

    msg_->DeclarePropertyWithUnit("decay_rate","nCi", N_Decays_per_s_,  "Rate of Decays in Bq");

    G4GenericMessenger::Command& NDecays_number_cmd =
            msg_->DeclareProperty("NDecays", NDecays, "Number of Decays");
    NDecays_number_cmd.SetParameterName("NDecays", false);
    NDecays_number_cmd.SetRange("NDecays > 0");
    // Load the detector geometry, which will be used for the generation of vertices
    const DetectorConstruction* detconst = dynamic_cast<const DetectorConstruction*>
    (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    if (detconst) geom_ = detconst->GetGeometry();
    else G4Exception("[CrabSourceGenerator]", "CrabSourceGenerator()", FatalException, "Unable to load geometry.");
}

//Destruct
////////////////////////////////////////////////////
CrabSourceGenerator::~CrabSourceGenerator(){
    delete msg_;
}
G4ParticleDefinition* CrabSourceGenerator::IonDefinition()
{
    G4ParticleDefinition* pdef =
            G4IonTable::GetIonTable()->GetIon(atomic_number_, mass_number_, energy_level_);

    if (!pdef) G4Exception("[CrabSourceGenerator]", "IonDefinition()",
                           FatalException, "Unable to find the requested ion.");

    // Unstable ions decay by default at a random time t sampled from an exponential
    // decay distribution proportional to their mean lifetime. This, even for
    // not so long-lived nuclides, pushes the global time of the event to scales
    // orders of magnitude longer than the nanosecond precision required in NEXUS
    // for the correct simulation of the ionization drift and photon tracing.
    // To prevent this behaviour, the lifetime of unstable isotopes is reset here
    // (unless the configuration variable 'decay_at_time_zero' has been set to false)
    // to an arbitrary, short value (1 ps).
    if (decay_at_time_zero_ && !(pdef->GetPDGStable())) pdef->SetPDGLifeTime(1.*ps);

    return pdef;
}
void CrabSourceGenerator::GeneratePrimaryVertex(G4Event* event)
{
    //Converting nCi to Bq 1nCi to 37 Bq
    G4int DecayRateinBq=(N_Decays_per_s_)*37*1e-6;
    G4int N_Decays=DecayRateinBq*event_window_;
    //Printing Some Values
    /*G4cout<<"-----Printing Setting Values----" <<G4endl;

    if(event_window_>0){
        G4cout<<"NumberOfDecays--> "<<N_Decays<<G4endl;
        G4cout<<"event_window(us)--> "<<event_window_<<G4endl;
        G4cout<<"N_Decays_per_s(nCi)--> "<<N_Decays_per_s_<<G4endl;
        G4cout<<"DecayRateInBq--> "<<DecayRateinBq<<G4endl;
    }
    else{
        G4cout<<"NumberOfDecays--> "<<NDecays<<G4endl;
    }*/
    //if time window is zero then run sim for specified number of Particles
    if(event_window_<=0){
        if (NDecays==0) NDecays=1;
        for(G4int i=0;i<NDecays;i++){
            decay_time_=event_window_;
            EventsWithWindow(event,decay_time_);
        }
    }else
    { //run sim for a given us time window

        for(G4int Decay=0;Decay<N_Decays;Decay++)
        {
            decay_time_=G4UniformRand()*event_window_;
            EventsWithWindow(event,decay_time_);
        }
    }

}
void CrabSourceGenerator::EventsWithWindow(G4Event*event,G4double decay_time){
    // Pointer declared as static so that it gets allocated only once
    // (i.e. the ion definition is only looked up in the first event).
    static G4ParticleDefinition* pdef = IonDefinition();
    // Create the new primary particle (i.e. the ion)
    G4PrimaryParticle* ion = new G4PrimaryParticle(pdef);

    // Generate an initial position for the ion using the geometry
    G4ThreeVector position = geom_->GenerateVertex(region_);
    // Ion generated at the start-of-event time
    // Create a new vertex
    //G4cout<<"ParticlePosition " <<position<< " DecayTime "<< decay_time<<G4endl;
    G4PrimaryVertex* vertex = new G4PrimaryVertex(position,decay_time);

    // Add ion to the vertex and this to the event
    vertex->SetPrimary(ion);
    event->AddPrimaryVertex(vertex);
}

