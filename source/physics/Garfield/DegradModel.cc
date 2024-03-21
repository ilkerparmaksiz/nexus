#include "DegradModel.h"
#include "IonizationElectron.h"
#include "PhysicsUtils.h"

#include <G4TransportationManager.hh>
#include <G4GlobalFastSimulationManager.hh>
#include <G4OpticalPhoton.hh>
#include <G4RandomDirection.hh>
#include <G4RunManager.hh>

using namespace nexus;

DegradModel::DegradModel(G4String modelName, G4Region* envelope, GarfieldHelper GH)
    : G4VFastSimulationModel(modelName, envelope), theFastIntegralTable_(0){
    
    GH_ = GH;
    event_id_ = -1;

    BuildThePhysicsTable(theFastIntegralTable_);

}

DegradModel::~DegradModel() {}

G4bool DegradModel::IsApplicable(const G4ParticleDefinition& particleType) {

    if (GH_.useDEGRAD_){
        if (particleType.GetParticleName()== "e-" || particleType.GetParticleName()== "gamma" || particleType.GetParticleName()== "e+")
            return true;
        else
            return false;
    }
    else {
        return false;
    }
}

G4bool DegradModel::ModelTrigger(const G4FastTrack& fastTrack) {
    
    //  also require that only photoelectric effect electrons are tracked here.
    // Krishan turn this off since now we are running the gammas
    G4int id = fastTrack.GetPrimaryTrack()->GetParentID();
    // if (id == 1){
    //     if ( (fastTrack.GetPrimaryTrack()->GetCreatorProcess()->GetProcessName().find("phot") != std::string::npos) ||
    //          (fastTrack.GetPrimaryTrack()->GetCreatorProcess()->GetProcessName().find("comp") != std::string::npos))
    //         return true;
    // }

    // Set the kinetic energy
    fPrimKE = fastTrack.GetPrimaryTrack()->GetKineticEnergy()/eV;

    // Dont simulate anything larger than 4 MeV
    if (fPrimKE> 4/eV){
        std::cout << "Primary particle energy larger than chosen Degrad limit, will use G4 generation" << std::endl;
        return false;
    }

    // Krishan: Degrad handles gammas/X-Rays but not for energies > 2 MeV
    // The photoelectron should be produced though from the >= 2 MeV gammas
    if (fastTrack.GetPrimaryTrack()->GetParticleDefinition()->GetParticleName() == "gamma" && fPrimKE>= 2/eV){
        std::cout << "Primary gamma energy larger than what Degrad can simulate, will use G4 generation" << std::endl;
        return false;
    }

    // We should only simulate particles if the interaction was on the xenon
    G4String solidName = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()
                                  ->LocateGlobalPointAndSetup(fastTrack.GetPrimaryTrack()->GetVertexPosition())->GetName();

    if (solidName != "ACTIVE")
        return false;

    return true; // Now return true
    // return false;

}

void DegradModel::DoIt(const G4FastTrack& fastTrack, G4FastStep& fastStep) {

    // To run degrad, we have to write the input parameters to a string
    // The string documentation is in the start of the degrad fortran file
    // The output file from degrad is then read back in

    // Reset variables if we are on a new event
    G4int evt =  G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    if (event_id_ !=  evt){
        event_id_  = evt;
        Reset();
    }

    // Get the current track ID
    G4int trk_id = fastTrack.GetPrimaryTrack()->GetTrackID();
    AddTrack(trk_id);

    // Start by killing G4's primary before it deposits any energy.
    fastStep.KillPrimaryTrack();

    // Initialization
    G4ThreeVector degradPos  = fastTrack.GetPrimaryTrack()->GetVertexPosition();
    G4double      degradTime = fastTrack.GetPrimaryTrack()->GetGlobalTime();
    fastStep.SetPrimaryTrackPathLength(0.0);
    G4cout << "GLOBAL TIME " << G4BestUnit(degradTime,"Time") << " POSITION " << G4BestUnit(degradPos,"Length") << G4endl;

    // Set the scintialltion components based on the xenon
    G4Material* mat = fastTrack.GetPrimaryTrack()->GetMaterial();
    G4MaterialPropertiesTable* mpt = mat->GetMaterialPropertiesTable();

    slow_comp_ = mpt->GetConstProperty("SCINTILLATIONTIMECONSTANT1");
    slow_prob_ = mpt->GetConstProperty("SCINTILLATIONYIELD1");
    fast_comp_ = mpt->GetConstProperty("SCINTILLATIONTIMECONSTANT2");
    fast_prob_ = mpt->GetConstProperty("SCINTILLATIONYIELD2");
    spectrum_integral = (G4PhysicsOrderedFreeVector*)(*theFastIntegralTable_)(mat->GetIndex());

    // Input parameters
    G4int    SEED = CLHEP::HepRandom::getTheSeed() + evt;
    G4String seed = G4UIcommand::ConvertToString(SEED);

    // Gamma KE
    std::cout <<"The particle energy is: " << fPrimKE*eV << "MeV" << std::endl;
    G4int KE = int(fPrimKE); // in eV
    G4String particleKE(","+std::to_string(KE));

    G4String particle_name = fastTrack.GetPrimaryTrack()->GetParticleDefinition()->GetParticleName();

    G4String degrad_mode = "0";
    if (particle_name == "gamma"){
        degrad_mode = "3";
        std::cout << "The particle is a gamma"<< std::endl;
    }
    // Assume its an electron then
    else{
        degrad_mode = "2";
        std::cout << "The particle is an electron/positron"<< std::endl;
    }
    
    std::cout << "The degrad mode is: " << degrad_mode << std::endl;

    // Xenon Pressure
    const static G4double torr = 1. / 750.062 * bar;
    G4int Press = GH_.GasPressure_/torr; // Use torr in Degrad
    std::cout << "The gas pressure is:" << Press<< std::endl;
    G4String xenonP(","+std::to_string(Press));
    
    // Electric field
    G4int Efield = GH_.fieldDrift_; // V/cm  
    G4String Efield_str(std::to_string(Efield));

    // Create the input card
    // Note the exact precision in below arguments. The integers gammaKE,xenonP in particular need a ".0" tacked on. 
    G4String degradString="printf \"1,1,"+degrad_mode+",2,"+seed+particleKE+".0,7.0,10000.0\n7,0,0,0,0,0\n100.0,0.0,0.0,0.0,0.0,0.0,20.0"+xenonP+".0\n"+Efield_str+".0,0.0,0.0,2,0\n100.0,0.5,0,0,1,1,1,1,1\n0,0,0,0,0,0\" > conditions_Degrad.txt";
    G4int stdout=system(degradString.data());

    // Execute degrad
    std::string degrad_exec = std::string(std::getenv("DEGRAD_HOME")) + "/Degrad < conditions_Degrad.txt"+ " > degrad_print.txt";
    const char *degrad_exec_str = degrad_exec.c_str();

    // Sometimes degrad just fails, so try running a couple times if it does
    // Usually 2 times is enough
    G4int return_status, tries{0};
    do {
        // Execute the command and get the return status
        return_status = stdout=system(degrad_exec_str);

        // If return status is not zero, print a message
        if (return_status != 0) {
            std::cout << "Command failed with return status: " << return_status << std::endl;
        }
        tries++;
    } while (return_status != 0 && tries < 4);  // Try three times before giving up

    std::cout << "Degrad return status: " << return_status << " number of tries: "<< tries << std::endl;
    stdout=system("cat degrad_print.txt");
    AddTrackLength(trk_id);

    
    // Convert file format
    std::string conv_path = "python3 " + std::string(std::getenv("NEXUSDIR")) + "/scripts/convertDegradFile.py";
    const char *conv_path_str = conv_path.c_str();
    stdout=system(conv_path_str);

    // Load the file and pop ie- and S1 to the stack
    GetElectronsFromDegrad(fastStep,degradPos,degradTime,trk_id);

}

void DegradModel::GetElectronsFromDegrad(G4FastStep& fastStep,G4ThreeVector degradPos,G4double degradTime, G4int trk_id)
{
    G4int     eventNumber, Nep, Nexc, nline, i, electronNumber, S1Number; // Nep is the number of primary es that corresponds to what biagi calls "ELECTRON CLUSTER SIZE (NCLUS)
    G4double  posX, posY, posZ, time, n;
    G4double  posXDegrad, posYDegrad, posZDegrad, timeDegrad;
    G4double  posXInitial{degradPos.getX()};
    G4double  posYInitial{degradPos.getY()};
    G4double  posZInitial{degradPos.getZ()};
    G4double  timeInitial{degradTime};
    G4double  Fluorescence{0}, PairProd{0}, Brems{0};
    G4String  line;
    std::vector<G4double> v;

    std::ifstream inFile;
    G4String fname= "DEGRAD.OUT";
    inFile.open(fname,std::ifstream::in);

    // Get the track index
    G4int trk_index = GetCurrentTrackIndex(trk_id);
    
    nline=1;
    electronNumber=0;
    S1Number=0;
    while (getline(inFile, line,'\n')) {
        
        std::istringstream iss(line);
        
        if (nline ==1) {
            while (iss >> n) {
                v.push_back(n);
            }

            eventNumber = v[0];
            Nep =v[1];
            Nexc = v[2];

            SetNioni(Nep, trk_index);
            G4cout << "Total ie-: " << Nep << G4endl;
            G4cout << "Total S1: " << Nexc << G4endl;
            v.clear();
        }
        // Ionizations
        if (nline ==2)  {
            
            fastStep.SetNumberOfSecondaryTracks(2E6); // reasonable max # of ie-/excitations created by degrad
            
            while (iss >> n) {
                v.push_back(n);
            }
            
            for (i=0;i<v.size();i=i+7){
                posXDegrad   = v[i];
                posYDegrad   = v[i+1];
                posZDegrad   = v[i+2];
                timeDegrad   = v[i+3];
                Fluorescence = v[i+4]; // 0 to N where N is the ie- for each N absorbed fluorescence photon in the event
                PairProd     = v[i+5]; // 0: not from pair prod, 1: produced from electron track,     2: produced from positron track
                Brems        = v[i+6]; // 0: not from brem,      1: produced from remaining electrons, 2: produced from brem gamma
                
                // Convert from um to mm in GEANT4
                // also Y and Z axes are swapped in GEANT4 and Garfield++ relatively to Degrad
                posX=posXDegrad*0.001+posXInitial;
                posY=posZDegrad*0.001+posYInitial;
                posZ=posYDegrad*0.001+posZInitial;
                time=timeDegrad*0.001+timeInitial; // Convert ps to ns

                // std::cout << "DegradModel::DoIt(): v[i-4]" << v[i] << "," << v[i+1] << "," << v[i+2] << "," << v[i+3] << "," << v[i+4]   << std::endl;
                // std::cout << "DegradModel::DoIt(): xinitial, poxXDegrad [mm]" << posXInitial << ", " << posXDegrad*0.001 << std::endl;
                
                G4ThreeVector myPoint(posX, posY, posZ);
                        
                // Set the track end position, the 3ns timing is in case there was a brem. May need to adjust
                if (Fluorescence == 0 && PairProd == 0 && (Brems == 0 || Brems == 1) && time < 3){
                    SetTrackEndPoint(myPoint, time, trk_index);
                }

                // std::cout << posX << ", " << posY << ", " << posZ << ", " << time << ", " << v[i+4] << ", " << v[i+5] << ", " << v[i+6]   << std::endl;

                electronNumber++;
                
                // Create secondary electron
                G4DynamicParticle electron(IonizationElectron::Definition(),G4RandomDirection(), 1.13*eV);
                G4Track *newTrack=fastStep.CreateSecondaryTrack(electron, myPoint, time,false);

            }
            v.clear();
        }
        // Excitations
        if (nline == 3)  {
            
            while (iss >> n) {
                v.push_back(n);
            }
            
            for (i=0;i<v.size();i=i+4){
                posXDegrad=v[i];
                posYDegrad=v[i+1];
                posZDegrad=v[i+2];
                timeDegrad=v[i+3];
                
                // Convert from um to mm in GEANT4
                // also Y and Z axes are swapped in GEANT4 and Garfield++ relatively to Degrad
                posX=posXDegrad*0.001+posXInitial;
                posY=posZDegrad*0.001+posYInitial;
                posZ=posYDegrad*0.001+posZInitial;
                time=timeDegrad*0.001+timeInitial +  GetScintTime(); // Convert ps to ns, add delay due to scintillation timing
                
                G4ThreeVector myPoint(posX, posY, posZ);
                
                S1Number++;
                
                // Create secondary electron
                auto* optphot = G4OpticalPhoton::OpticalPhotonDefinition();

                G4ThreeVector momentum, polarization;
                GetPhotonPol(momentum, polarization);

                // Get Photon energy
                G4double sc_max = spectrum_integral->GetMaxValue();
                G4double sc_value = G4UniformRand()*sc_max;
                G4double sampled_energy = spectrum_integral->GetEnergy(sc_value);
                
                G4DynamicParticle VUVphoton(optphot, momentum, sampled_energy);

                G4Track *newTrack=fastStep.CreateSecondaryTrack(VUVphoton, myPoint, time ,false);
                newTrack->SetPolarization(polarization);

            }
            v.clear();
            
        }

        nline++;

    }

    inFile.close();
    G4cout << "Number of active ie-: " << electronNumber << G4endl;
    G4cout << "Number of active S1: " << S1Number << G4endl;
    
}

void DegradModel::SetTrackEndPoint(G4ThreeVector pos, G4double time, G4int trk_index){

    if (time > end_times_[trk_index]){
        end_times_[trk_index] = time;
        track_end_pos_[trk_index] = pos;
    }

}

void DegradModel::SetNioni(G4int Ne, G4int trk_index){

    N_ioni_[trk_index] =  Ne;

}

G4ThreeVector DegradModel::GetTrackEndPoint(G4int trk_id){
    G4int trk_index = GetCurrentTrackIndex(trk_id);
    return G4ThreeVector(track_end_pos_[trk_index]);
}

G4double DegradModel::GetTrackEndTime(G4int trk_id){
    G4int trk_index = GetCurrentTrackIndex(trk_id);
    return end_times_[trk_index];
}

G4double DegradModel::GetTrackLength(G4int trk_id){
    G4int trk_index = GetCurrentTrackIndex(trk_id);
    return trk_len_vec_[trk_index];
}

void DegradModel::Reset(){
    end_times_.clear();
    track_end_pos_.clear();
    track_ids_.clear();
    ke_vec_.clear();
    N_ioni_.clear();
    trk_len_vec_.clear();

}

void DegradModel::AddTrack(G4int trk_id){

    // Check if track exists in vec, if not then add
    if (std::find(track_ids_.begin(), track_ids_.end(), trk_id) == track_ids_.end()) {
        track_ids_.push_back(trk_id);
        end_times_.push_back(-1);
        track_end_pos_.push_back(G4ThreeVector(0,0,0));
        ke_vec_.push_back(fPrimKE);
        N_ioni_.push_back(0);
        trk_len_vec_.push_back(0.);
    }
}

G4int DegradModel::GetCurrentTrackIndex(G4int trk_id){

    G4int index =-1;

    // Search for searchValue in vec
    auto it = std::find(track_ids_.begin(), track_ids_.end(), trk_id);

    // Check if value was found
    if (it != track_ids_.end()) {
        // Calculate the index
        index = std::distance(track_ids_.begin(), it);
    }

    // if (index == -1)
        // std::cout <<"DegradModel: Warning negative index " << trk_id << std::endl;

    return index;

}

G4double DegradModel::GetAvgIoniEnergy(G4int trk_id){

    G4int trk_index = GetCurrentTrackIndex(trk_id);

    return ke_vec_[trk_index]/N_ioni_[trk_index];

}

G4int DegradModel::GetTotIonizations(G4int trk_id){

    G4int trk_index = GetCurrentTrackIndex(trk_id);

    return N_ioni_[trk_index];

}

void DegradModel::AddTrackLength(G4int trk_id){

    // Please note, it appears the track length in degrad is actually the dist
    // to the track end point from start rather than the prim length. 

    std::cout << "Adding Track Length from Degrad" << std::endl;

    G4int trk_index = GetCurrentTrackIndex(trk_id);

    G4double length = 0;

    // Open the file for reading
    std::ifstream inputFile("degrad_print.txt");
    if (!inputFile.is_open()) {
        std::cerr << "Failed to open the file!" << std::endl;
        return;
    }

    // Search for the line containing "AV.MAX. RANGE IN XYZ="
    G4String line;
    G4String targetString = "AV.MAX. RANGE IN XYZ=";
    while (std::getline(inputFile, line)) {
        
        size_t pos = line.find(targetString);
        
        if (pos != G4String::npos) {
            // Extract the value after the "=" sign
            G4String value = line.substr(pos + targetString.length());
            // Replace "D" with "e"
            size_t pos = value.find('D');
            if (pos != G4String::npos) {
                value.replace(pos, 1, "e");
            }

            // Convert the modified string to a double
            length = std::stod(value);
            break;  // Stop searching after finding the line
        }
    }

    trk_len_vec_[trk_index] = length*0.001; // degrad length is in um

    // Close the file
    inputFile.close();

}

G4double DegradModel::GetScintTime(){

  G4double scint_time = 0;

  // Generate a random number to determine which distribution to sample from
  double randomNumber = G4UniformRand();

  // Fast Component - 4.5 ns
  if (randomNumber < slow_prob_) {
      scint_time = G4RandExponential::shoot(slow_comp_);
  // Slow Component - 100 ns
  } else {
      scint_time = G4RandExponential::shoot(fast_comp_);
  }

  return scint_time;
}
