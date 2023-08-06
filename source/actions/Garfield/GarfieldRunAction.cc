//Runaction
#include "GarfieldRunAction.hh"
#include "G4Run.hh"
#include "Randomize.hh"
#include "PrimaryGeneratorAction.hh"
#include "Analysis.hh"
#include "GarfieldEventAction.hh"

#include "GasBoxSD.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "FactoryBase.h"
using namespace nexus;

REGISTER_CLASS(GarfieldRunAction,G4UserRunAction)

GarfieldRunAction::GarfieldRunAction(){
  G4cout << "Creating AnalysisManager" << G4endl;
  auto analysisManager = G4AnalysisManager::Instance();
//  analysisManager->SetNtupleMerging(true,0,0,10000000);
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(true);
  analysisManager->SetFileName("output.root");
  analysisManager->SetHistoDirectoryName("histo");
  analysisManager->SetNtupleDirectoryName("ntuple");

  analysisManager->SetNtupleActivation(true);

  G4cout << "Creating GarfieldRunAction" << G4endl;
}

GarfieldRunAction::~GarfieldRunAction() {
	G4cout << "Deleting GarfieldRunAction" << G4endl;
	delete G4AnalysisManager::Instance();
}


void GarfieldRunAction::BeginOfRunAction(const G4Run* run) {
  G4Random::showEngineStatus();

    G4cout << "### Run " << run->GetRunID() << " start." << G4endl;
  time_t currentTime;
  tm* ptm;
  time(&currentTime);
  ptm = localtime(&currentTime);
  G4cout << "Time: " << asctime(ptm) << G4endl;

  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->OpenFile();

  analysisManager->CreateNtuple("Camera", "Camera Hits");
  analysisManager->CreateNtupleDColumn("Event");     //column 0
  analysisManager->CreateNtupleDColumn("PID");       //column 1
  analysisManager->CreateNtupleDColumn("Time");      //column 2
  analysisManager->CreateNtupleDColumn("X");         //column 3
  analysisManager->CreateNtupleDColumn("Y");         //column 4
  analysisManager->CreateNtupleDColumn("Z");         //column 5
  analysisManager->CreateNtupleIColumn("Reflected"); //column 6
  analysisManager->CreateNtupleSColumn("Boundary");  //column 7
  analysisManager->CreateNtupleIColumn("SID");       //column 8

  analysisManager->FinishNtuple();

  analysisManager->CreateNtuple("S1", "S1 particles");
  analysisManager->CreateNtupleDColumn("Event");    //column 0
  analysisManager->CreateNtupleDColumn("PID");      //column 1
  analysisManager->CreateNtupleDColumn("Time");     //column 2
  analysisManager->CreateNtupleDColumn("X");        //column 3
  analysisManager->CreateNtupleDColumn("Y");        //column 4
  analysisManager->CreateNtupleDColumn("Z");        //column 5
  analysisManager->CreateNtupleSColumn("Process");  //column 6
  analysisManager->FinishNtuple();

  analysisManager->CreateNtuple("S2", "S2 particles");
  analysisManager->CreateNtupleDColumn("Event");    //column 0
  analysisManager->CreateNtupleDColumn("PID");      //column 1
  analysisManager->CreateNtupleDColumn("Time");     //column 2
  analysisManager->CreateNtupleDColumn("X");        //column 3
  analysisManager->CreateNtupleDColumn("Y");        //column 4
  analysisManager->CreateNtupleDColumn("Z");        //column 5
  analysisManager->CreateNtupleSColumn("Process");  //column 6
  analysisManager->FinishNtuple();

  analysisManager->CreateNtuple("Event", "Event stats");
  analysisManager->CreateNtupleDColumn("Event");    //column 0
  analysisManager->CreateNtupleDColumn("PPID");     //column 1 Primary PPartCode
  analysisManager->CreateNtupleDColumn("PKE");      //column 2 Primary KE
  analysisManager->CreateNtupleDColumn("PEDep");     //column 3 Primary final KE
  analysisManager->FinishNtuple();

  analysisManager->CreateNtuple("Lens", "Lens Hits");
  analysisManager->CreateNtupleDColumn("Event");     //column 0
  analysisManager->CreateNtupleDColumn("PID");       //column 1
  analysisManager->CreateNtupleDColumn("Time");      //column 2
  analysisManager->CreateNtupleDColumn("X");         //column 3
  analysisManager->CreateNtupleDColumn("Y");         //column 4
  analysisManager->CreateNtupleDColumn("Z");         //column 5
  analysisManager->FinishNtuple();

  analysisManager->CreateNtuple("PMT", "PMT Hits");
  analysisManager->CreateNtupleDColumn("Event");     //column 0
  analysisManager->CreateNtupleDColumn("PID");       //column 1
  analysisManager->CreateNtupleDColumn("Time");      //column 2
  analysisManager->CreateNtupleDColumn("X");         //column 3
  analysisManager->CreateNtupleDColumn("Y");         //column 4
  analysisManager->CreateNtupleDColumn("Z");         //column 5
  analysisManager->CreateNtupleIColumn("Reflected"); //column 6
  analysisManager->CreateNtupleSColumn("Boundary");  //column 7
  analysisManager->CreateNtupleIColumn("SID");       //column 8
  analysisManager->FinishNtuple();

  analysisManager->SetNtupleActivation(true);
}

void GarfieldRunAction::EndOfRunAction(const G4Run* run) {
    G4cout << "### Run " << run->GetRunID() << " end." << G4endl;
    auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetCompressionLevel(3);
  analysisManager->Write();
  analysisManager->CloseFile();

  G4cout << "End of run OK!" << G4endl;
  time_t currentTime;
  tm* ptm;
  time(&currentTime);
  ptm = localtime(&currentTime);
  G4cout << "Simulation finished." << G4endl << "Time: " << asctime(ptm)
         << G4endl;
  G4Random::showEngineStatus();


}
