// ----------------------------------------------------------------------------
// nexus | DefaultRunAction.cc
//
// This is the default run action of the NEXT simulations.
// A message at the beginning and at the end of the simulation is printed.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "DefaultRunAction.h"
#include "FactoryBase.h"
#include "IOUtils.h"
#include <G4Run.hh>
#include "PersistencyManager.h"
using namespace nexus;

REGISTER_CLASS(DefaultRunAction, G4UserRunAction)

DefaultRunAction::DefaultRunAction(): G4UserRunAction()
{
}



DefaultRunAction::~DefaultRunAction()
{
}



void DefaultRunAction::BeginOfRunAction(const G4Run* run)
{
    startTime=std::chrono::high_resolution_clock::now();
    Runtime=0;

    //PersistencyManager * pManger=dynamic_cast<PersistencyManager *> (PersistencyManager::GetPersistencyManager());
    //pManger->SetFstream(OpenFile(nullptr,"time.txt"));
    G4cout << "### Run " << run->GetRunID() << " start." << G4endl;

}


void DefaultRunAction::EndOfRunAction(const G4Run* run)
{

    // Calculate the duration
   //CloseFile();
    // Print the duration in seconds
  G4cout << "### Run " << run->GetRunID() << " end." << G4endl;
  endTime=std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration = endTime - startTime;
  Runtime=duration.count();

 // PersistencyManager * pManger=dynamic_cast<PersistencyManager *> (PersistencyManager::GetPersistencyManager());

  //CloseFile(pManger->GetFstream());

  std::cout << "Run time: " << Runtime<< " seconds." << std::endl;


}
