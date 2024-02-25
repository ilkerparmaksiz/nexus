#include "G4Electron.hh"
#include "G4SystemOfUnits.hh"
#include "DegradModel.h"
#include "G4Region.hh"
#include "G4ParticleDefinition.hh"
#include "G4UnitsTable.hh"
#include "G4Track.hh"
#include "Randomize.hh"
#include "G4UIcommand.hh"
#include <fstream>
#include "G4TransportationManager.hh"
#include "G4DynamicParticle.hh"
#include "G4RandomDirection.hh"
#include "G4VProcess.hh"
// #include "NESTProc.hh"
#include "IonizationElectron.h"
#include "G4GenericMessenger.hh"
#include <G4OpticalPhoton.hh>

using namespace nexus;

DegradModel::DegradModel(G4String modelName, G4Region* envelope, GarfieldHelper GH)
    : G4VFastSimulationModel(modelName, envelope){
    
    GH_ = GH;

    processOccured = false;


}

DegradModel::~DegradModel() {}

G4bool DegradModel::IsApplicable(const G4ParticleDefinition& particleType) {

    if (GH_.useDEGRAD_){
        if (particleType.GetParticleName()=="e-" || particleType.GetParticleName()=="gamma" || particleType.GetParticleName()=="e+")
            return true;
        else
            return false;
    }
    else {
        return false;
    }
}

G4bool DegradModel::ModelTrigger(const G4FastTrack& fastTrack) {
    G4int id = fastTrack.GetPrimaryTrack()->GetParentID();
    
    //  also require that only photoelectric effect electrons are tracked here.
    // Krishan turn this off since now we are running the gammas
    // if (id == 1){
    //     if ( (fastTrack.GetPrimaryTrack()->GetCreatorProcess()->GetProcessName().find("phot") != std::string::npos) ||
    //          (fastTrack.GetPrimaryTrack()->GetCreatorProcess()->GetProcessName().find("comp") != std::string::npos))
    //         return true;
    // }

    // Krishan: Degrad handles gammas/X-Rays but not for energies > 2 MeV
    // The photoelectron should be produced though from the >= 2 MeV gammas
    if (fastTrack.GetPrimaryTrack()->GetParticleDefinition()->GetParticleName() == "gamma" && fPrimKE>= 2/eV){
        std::cout << "Primary gamma energy larger than what Degrad can simulate, will use G4 generation" << std::endl;
        return false;
    }

    return true; // Now return true
    // return false;

}

void DegradModel::DoIt(const G4FastTrack& fastTrack, G4FastStep& fastStep) {

    // Here we start by killing G4's naive little one photo-electron.
    // Then we run degrad with a photon of desired energy. Then we read up all the electrons it produces.

    // To run degrad, we have to write the input parameters to a string
    // The string documentation is in the start of the degrad fortran file
    // The output file from degrad is then read back in

    fastStep.KillPrimaryTrack();
    if(!processOccured){
        
        // Initialization
        G4ThreeVector degradPos =fastTrack.GetPrimaryTrack()->GetVertexPosition();
        G4double degradTime = fastTrack.GetPrimaryTrack()->GetGlobalTime();
        fastStep.SetPrimaryTrackPathLength(0.0);
        G4cout<<"GLOBAL TIME "<<G4BestUnit(degradTime,"Time")<<" POSITION "<<G4BestUnit(degradPos,"Length")<<G4endl;
        
        
        // Input parameters
        G4int SEED=54217137*G4UniformRand();
        G4String seed = G4UIcommand::ConvertToString(SEED);

        // Gamma KE
        G4int KE = int(fPrimKE); // in eV
        G4String particleKE(","+std::to_string(KE));

        G4String particle_name = fastTrack.GetPrimaryTrack()->GetParticleDefinition()->GetParticleName();

        G4String degrad_mode = "0";
        if (particle_name == "gamma")
            degrad_mode = "3";
        // Assume its an electron then
        else
            degrad_mode = "2";
        
        std::cout << "The degrad mode is: " << degrad_mode << std::endl;

        // Xenon Pressure
        const static G4double torr = 1. / 750.062 * bar;
        G4int Press = GH_.GasPressure_/torr;  // Krishan: need to make sure the gas pressure is in the right units here -- Degrad takes torr units
        std::cout << "The gas pressure is:" << Press<< std::endl;
        G4String xenonP(","+std::to_string(Press));
        
        // Electric field
         G4int Efield = 500; // V/cm  
         G4String Efield_str(std::to_string(Efield));

        // Create the input card
        // Note the exact precision in below arguments. The integers gammaKE,xenonP in particular need a ".0" tacked on. 
        G4String degradString="printf \"1,1,"+degrad_mode+",2,"+seed+particleKE+".0,7.0,10000.0\n7,0,0,0,0,0\n100.0,0.0,0.0,0.0,0.0,0.0,20.0"+xenonP+".0\n"+Efield_str+".0,0.0,0.0,2,0\n100.0,0.5,0,0,1,1,1,1,1\n0,0,0,0,0,0\" > conditions_Degrad.txt";
        // G4String degradString="printf \"1,1,3,-1,"+seed+gammaKE+".0,7.0,0.0\n7,0,0,0,0,0\n100.0,0.0,0.0,0.0,0.0,0.0,20.0"+xenonP+".0\n"+Efield_str+".0,0.0,0.0,2,0\n100.0,0.5,1,1,1,1,1,1,1\n0,0,0,0,0,0\" > conditions_Degrad.txt";
        // Execute
        G4int stdout=system(degradString.data());
        G4cout << degradString << G4endl;
        const std::string degradpath = std::getenv("DEGRAD_HOME");
        G4cout << degradpath << G4endl;
        std::string exec = "/Degrad < conditions_Degrad.txt";
        std::string full_path = degradpath + exec;
        const char *mychar = full_path.c_str();
        G4cout << mychar << G4endl;
        stdout=system(mychar);
        stdout=system("./scripts/convertDegradFile.py");

        GetElectronsFromDegrad(fastStep,degradPos,degradTime);
        
        // We call Degrad only once, which now that we have the x,y,z location of our primary Xray interaction, re-simulates that interaction. 
        // Note the 5900 in the Degrad config file. The above system() line forces the single Degrad execution. EC, 2-Dec-2021.
        // Krishan - Turn this off so we can try to generate multiple particles above threshold
        // processOccured=true; 
    }

}

void DegradModel::GetElectronsFromDegrad(G4FastStep& fastStep,G4ThreeVector degradPos,G4double degradTime)
{
    G4int eventNumber,Nep, nline, i, electronNumber, S1Number; //Nep is the number of primary es that corresponds to what biagi calls "ELECTRON CLUSTER SIZE (NCLUS)
    G4double posX,posY,posZ,time,n;
    G4double  posXDegrad,posYDegrad,posZDegrad,timeDegrad;
    G4double  posXInitial=degradPos.getX();
    G4double  posYInitial=degradPos.getY();
    G4double  posZInitial=degradPos.getZ();
    G4double  timeInitial=degradTime;
    G4String line;
    std::vector<G4double> v;
    
    std::ifstream inFile;
    G4String fname= "DEGRAD.OUT";
    inFile.open(fname,std::ifstream::in);
    
    G4cout<< "Working in "<<fname<<G4endl;
    
    nline=1;
    electronNumber=0;
    S1Number=0;
    while (getline(inFile, line,'\n')) {
        
        std::istringstream iss(line);//stream de strings
        
        if (nline ==1) {
            while (iss >> n) {
                v.push_back(n);
            }
            
            eventNumber=v[0];
            Nep=v[1];
            //  Nexc=v[2];
            v.clear();
        }
        //Ionizations
        if (nline ==2)  {
            
            fastStep.SetNumberOfSecondaryTracks(1E6); // reasonable max # of electrons created by degrad
            
            while (iss >> n) {
                v.push_back(n); //o n é adicionado ao vector
            }
            
            for (i=0;i<v.size();i=i+7){
                posXDegrad=v[i];
                posYDegrad=v[i+1];
                posZDegrad=v[i+2];
                timeDegrad=v[i+3];
                //convert from um to mm in GEANT4
                //also Y and Z axes are swapped in GEANT4 and Garfield++ relatively to Degrad
                posX=posXDegrad*0.001+posXInitial;
                posY=posZDegrad*0.001+posYInitial;
                posZ=posYDegrad*0.001+posZInitial;
                // std::cout << "DegradModel::DoIt(): v[i-4]" << v[i] << "," << v[i+1] << "," << v[i+2] << "," << v[i+3] << "," << v[i+4]   << std::endl;
                //std::cout << "DegradModel::DoIt(): xinitial, poxXDegrad [mm]" << posXInitial << ", " << posXDegrad*0.001 << std::endl;
    
                //convert ps to ns
                time=timeDegrad*0.001+timeInitial;
                
                
                G4ThreeVector myPoint;
                myPoint.setX(posX);
                myPoint.setY(posY);
                myPoint.setZ(posZ);
                
                //Check in which Physical volume the point bellongs
                G4Navigator* theNavigator= G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
                
                G4VPhysicalVolume* myVolume = theNavigator->LocateGlobalPointAndSetup(myPoint);
              
                G4String solidName=myVolume->GetName();
                
                if (G4StrUtil::contains(solidName,"ACTIVE") ){

                    electronNumber++;
                    
                    // Create secondary electron
                    G4DynamicParticle electron(IonizationElectron::Definition(),G4RandomDirection(), 1.13*eV);
                    G4Track *newTrack=fastStep.CreateSecondaryTrack(electron, myPoint, time,false);

                    // }
                }
            }
            v.clear(); //Faz reset ao vector senão vai continuar a adicionar os dadosadicionar os dados
            // nline=0;
            
        }
        // Excitations
        if (nline == 3)  {
            
            while (iss >> n) {
                v.push_back(n); //o n é adicionado ao vector
            }

            std::cout << v.size() << std::endl;
            
            for (i=0;i<v.size();i=i+4){
                posXDegrad=v[i];
                posYDegrad=v[i+1];
                posZDegrad=v[i+2];
                timeDegrad=v[i+3];
                //convert from um to mm in GEANT4
                //also Y and Z axes are swapped in GEANT4 and Garfield++ relatively to Degrad
                posX=posXDegrad*0.001+posXInitial;
                posY=posZDegrad*0.001+posYInitial;
                posZ=posYDegrad*0.001+posZInitial;
                // std::cout << "DegradModel::DoIt(): v[i-4]" << v[i] << "," << v[i+1] << "," << v[i+2] << "," << v[i+3] << "," << v[i+4]   << std::endl;
                //std::cout << "DegradModel::DoIt(): xinitial, poxXDegrad [mm]" << posXInitial << ", " << posXDegrad*0.001 << std::endl;
    
                //convert ps to ns
                time=timeDegrad*0.001+timeInitial;
                
                
                G4ThreeVector myPoint;
                myPoint.setX(posX);
                myPoint.setY(posY);
                myPoint.setZ(posZ);
                
                //Check in which Physical volume the point bellongs
                G4Navigator* theNavigator= G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
                
                G4VPhysicalVolume* myVolume = theNavigator->LocateGlobalPointAndSetup(myPoint);
              
                G4String solidName=myVolume->GetName();
                
                if (G4StrUtil::contains(solidName,"ACTIVE")){

                    S1Number++;
                    
                    // Create secondary electron
                    auto* optphot = G4OpticalPhoton::OpticalPhotonDefinition();
                    G4DynamicParticle VUVphoton(optphot,G4RandomDirection(), 7.2*eV);
                    G4Track *newTrack=fastStep.CreateSecondaryTrack(VUVphoton, myPoint, time ,false);
                    newTrack->SetPolarization(G4ThreeVector(0.,0.,1.0)); // Needs some pol'n, else we will only ever reflect at an OpBoundary. EC, 8-Aug-2022.

                    // }
                }
            }
            v.clear(); //Faz reset ao vector senão vai continuar a adicionar os dadosadicionar os dados
            nline=0;
            
        }


        nline++;
        
        
    }
    inFile.close();
    G4cout << "Number of initial electrons: " << electronNumber << G4endl;
    G4cout << "Number of initial S1: " << S1Number << G4endl;
    
    
}


