#include "G4Electron.hh"
#include "G4SystemOfUnits.hh"
#include "GarfieldVUVPhotonModel.h"
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
#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4OpticalPhoton.hh"
#include "G4ProcessManager.hh"
#include "G4EventManager.hh"
#include <G4GlobalFastSimulationManager.hh>

#include "globals.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/ComponentUser.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/Medium.hh"
#include "Garfield/SolidTube.hh"
#include "Garfield/ComponentComsol.hh"
#include "config.h"
#include "Trajectory.h"
#include "TrajectoryMap.h"
#include "XenonProperties.h"
#include "PhysicsUtils.h"


namespace nexus{

const G4double res(0.01); // Estimated fluctuations in EL yield - high? EC, 21-June-2022.

GarfieldVUVPhotonModel::GarfieldVUVPhotonModel(G4String modelName,G4Region* envelope, GarfieldHelper GH,  IonizationSD* ionisd) : G4VFastSimulationModel(modelName, envelope), fGarfieldSD(ionisd), theFastIntegralTable_(0) {
    
    GH_ = GH;
    GH_.DumpParams();
    ionMobFile = "IonMobility_Xe+_P12_Xe.txt";
    gasFile = "Xenon_10Bar.gas";
    InitialisePhysics();
    BuildThePhysicsTable(theFastIntegralTable_);
    dm_ = (DegradModel*)(G4GlobalFastSimulationManager::GetInstance()->GetFastSimulationModel("DegradModel"));
    event_id_ = -1;
    theNavigator_ = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  
}

G4bool GarfieldVUVPhotonModel::IsApplicable(const G4ParticleDefinition& particleType) {
  //  std::cout << "GarfieldVUVPhotonModel::IsApplicable() particleType is " << particleType.GetParticleName() << std::endl;
  
  if (particleType.GetParticleName()=="ie-") // || particleType.GetParticleName()=="opticalphoton")
    return true;
  return false;
}

G4bool GarfieldVUVPhotonModel::ModelTrigger(const G4FastTrack& fastTrack){
 
  G4double ekin = fastTrack.GetPrimaryTrack()->GetKineticEnergy();
  //  std::cout << "GarfieldVUVPhotonModel::ModelTrigger() thermalE, ekin is " << GH_.thermalE_/eV << ",  "<< ekin / eV << std::endl;
  
  G4String particleName = fastTrack.GetPrimaryTrack()->GetParticleDefinition()->GetParticleName();

  // std::cout << "The particle name is: " << particleName <<  std::endl;

   if (ekin < GH_.thermalE_ && particleName=="ie-")
    {
      return true;
    }
  return false;

} 
    
void GarfieldVUVPhotonModel::DoIt(const G4FastTrack& fastTrack, G4FastStep& fastStep) 
{

  /* 
     This tracks all the ionization/conversion electrons created by G4/NEST/Degrad in its simulation.
     Each such electron is drifted in the E-field. We create S2 photons if the 
     electrons reach the EL region. These will be tracked in the normal way with G4 or via optics. 

   */
    fastStep.SetNumberOfSecondaryTracks(1E3);

    fastStep.KillPrimaryTrack(); //KILL NEST/DEGRAD/G4 TRACKS

    // Reset variables if we are on a new event
    G4int evt =  G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    if (event_id_ !=  evt){
        event_id_  = evt;
        Reset();
    }

    G4int parent_id = fastTrack.GetPrimaryTrack()->GetParentID();
    AddTrack(parent_id);    
    garfPos         = fastTrack.GetPrimaryTrack()->GetVertexPosition();
    garfTime        = fastTrack.GetPrimaryTrack()->GetGlobalTime();
    G4int parent_track_idx = GetCurrentTrackIndex(parent_id);

    //G4cout<<"GLOBAL TIME "<<G4BestUnit(garfTime,"Time")<<" POSITION "<<G4BestUnit(garfPos,"Length")<<G4endl;
    
    N_ioni_[parent_track_idx]++;
        
    // Print how many of each type we have
    if (!(N_ioni_[parent_track_idx]%1000)){
      G4cout << "GarfieldVUV: ie-: " << N_ioni_[parent_track_idx] << " for track: " << parent_id << G4endl;
      
    }

    if (!(N_S2_[parent_track_idx]%1000) and (N_S2_[parent_track_idx]>0)){
      G4cout << "GarfieldVUV: S2: "  << N_S2_[parent_track_idx]   << " for track: " << parent_id << G4endl;
    }

    G4double mean_ioni_E = 22.4*eV;

    // Add the energy deposited to the trajectory so the event gets stored
    Trajectory* trj = (Trajectory*) TrajectoryMap::Get(parent_id);
    if (trj) {
    
      // Check if the parent was made by degrad and set mean energy deposit
      if (dm_ && dm_->GetCurrentTrackIndex(parent_id) != -1){
        mean_ioni_E = dm_->GetAvgIoniEnergy(parent_id)*eV;
        trj->SetEnergyDeposit(mean_ioni_E * dm_->GetTotIonizations(parent_id)); // This overwrites the same track
      }
     
    }

    // Set the scintillation timing
    if (N_ioni_[parent_track_idx] == 1){
      // Krishan: need to throw an exception if the medium is not GXe

      G4Material* mat = fastTrack.GetPrimaryTrack()->GetMaterial();
      G4MaterialPropertiesTable* mpt = mat->GetMaterialPropertiesTable();
      // std::cout << "The material is: " << mat->GetName()<< std::endl;

      slow_comp_ = mpt->GetConstProperty("SCINTILLATIONTIMECONSTANT1");
      slow_prob_ = mpt->GetConstProperty("SCINTILLATIONYIELD1");
      fast_comp_ = mpt->GetConstProperty("SCINTILLATIONTIMECONSTANT2");
      fast_prob_ = mpt->GetConstProperty("SCINTILLATIONYIELD2");
      spectrum_integral = (G4PhysicsOrderedFreeVector*)(*theFastIntegralTable_)(mat->GetIndex());
      // std::cout << "Printing timing properties " << slow_comp_ << ", " << fast_comp_ << ", " << slow_prob_ << std::endl;

    }

    // Drift and make S2
    GenerateVUVPhotons(fastTrack,fastStep,garfPos,garfTime, parent_id, mean_ioni_E);

}

void GarfieldVUVPhotonModel::GenerateVUVPhotons(const G4FastTrack& fastTrack, G4FastStep& fastStep, G4ThreeVector garfPos, G4double garfTime, G4int trk_id, G4double mean_ioni_E) {

  G4double x0=garfPos.getX()/cm; //Garfield length units are in cm
  G4double y0=garfPos.getY()/cm;
  G4double z0=garfPos.getZ()/cm;
  G4double t0=garfTime;

  // std::cout << "GVUVPM: positions are " << x0<<"," <<y0<<","<<z0 <<"," <<t0<< std::endl;

  // Check in which Physical volume the point bellongs
  G4ThreeVector myPoint(x0, y0, z0);

  //Check in which Physical volume the point bellongs
  G4VPhysicalVolume* myVolume = theNavigator_->LocateGlobalPointAndSetup(myPoint);
  G4String solidName = myVolume->GetName();

  // Ensure that we only drift electrons in the active volume
  if (solidName != "ACTIVE")
    return;

  // Insert hits
  if (dm_ && dm_->GetCurrentTrackIndex(trk_id) != -1){
    InsertHits(x0, y0, z0, t0, trk_id, mean_ioni_E);
  }
  // return;

  // Print Electric Field
  // PrintElectricField(x0, y0, z0);

  G4double xi,yi,zi,ti;
  
  // AvalancheMC drift for drift, and then create excitations when reaching the EL.
  fAvalancheMC->DriftElectron(x0,y0,z0,t0);
  size_t n =fAvalancheMC->GetElectrons().at(0).path.size();
  auto GetDriftLines=fAvalancheMC->GetElectrons().at(0).path;

  // Get zi when in the beginning of the EL region
  for(G4int i=0;i<n;i++){

    xi=GetDriftLines.at(i).x;
    yi=GetDriftLines.at(i).y;
    zi=GetDriftLines.at(i).z;
    ti=GetDriftLines.at(i).t;
    // std::cout << "GVUVPM: positions are " << xi<<"," <<yi<<","<<zi <<"," <<ti<< std::endl;

    G4bool b_xy_bounds;

    // Do a check against a polygon
    if (GH_.nsides_> 1){
      std::vector<G4double> point = {xi,yi};
      b_xy_bounds = CheckXYBoundsPolygon(point);
    }
    // Check against cylinder
    else {
       b_xy_bounds = std::sqrt(xi*xi + yi*yi) < GH_.DetActiveR_;
    }

    // Drift line point entered EL region
    if (zi < GH_.ELPos_ && b_xy_bounds ){
      
      // Check if attachment killed the ie-
      if (!GetAttachment(ti))
        return;
      
      break; 
    }
    
    // No drift line point meets criteria, so return
    else if (i==n-1)
      return;

  }

  // Generate the El photons from a microphys model ran externally in Garfield
  // We sample the output file which contains the timing profile of emission and diffusion
  if (GH_.useELFile_)
      MakeELPhotonsFromFile(fastStep, xi, yi, zi, ti, trk_id);
  // Use a simpler model
  else
      MakeELPhotonsSimple(fastStep, xi, yi, zi, ti, trk_id);

}


// Selection of Xenon exitations and ionizations
void GarfieldVUVPhotonModel::InitialisePhysics(){

    // Set the path
    char* path = std::getenv("NEXUSDIR");
    if (path == nullptr) {
      G4Exception("[GarfieldVUVPhotonModel]", "InitializePhysics()", FatalException,
              "Environment variable NEXUSDIR not defined!");
    }

    G4String nexus_path(path);

    // Set the gas Properties
    fMediumMagboltz = new Garfield::MediumMagboltz();
    fMediumMagboltz->SetComposition("Xe", 100.);
    fMediumMagboltz->DisableDebugging();
    
    //  --- Load in Ion Mobility file --- 
    const std::string garfpath = getenv("GARFIELD_HOME");
    
    if(ionMobFile!="")
      fMediumMagboltz->LoadIonMobility(garfpath + "/Data/" + ionMobFile);
    
    if(gasFile!=""){
      fMediumMagboltz->LoadGasFile(nexus_path + "/data/" + gasFile.c_str());
      std::cout << "Loaded gasfile "<< gasFile << std::endl;
    }

    // Initialize the gas
    fMediumMagboltz->Initialise(true);

    // Print the gas properties
    // fMediumMagboltz->PrintGas();

    std::cout << "Detector Dimentions: "<< GH_.DetChamberR_ << " " << GH_.DetChamberL_ << "  " << GH_.DetActiveR_ << "  " << GH_.DetActiveL_ << std::endl; 

    fSensor = new Garfield::Sensor();

    if (!GH_.useCOMSOL_){

        // The drift tube is a polygon not cylinder
        if (GH_.nsides_ > 1)
          InitalizePolygon();

        Garfield::ComponentUser* componentDriftLEM = SetComponentField();
        fSensor->AddComponent(componentDriftLEM);
        
        // Set the region where the sensor is active -- based on the gas volume
        // fSensor->SetArea(-GH_.DetChamberR_, -GH_.DetChamberR_, -GH_.DetChamberL_, GH_.DetChamberR_, GH_.DetChamberR_, GH_.CathodePos_); // cm

    }
    else {
        std::cout << "Initialising Garfiled with a COMSOL geometry" << std::endl;
        G4String home = nexus_path + "/data/Garfield/";
        std::string meshfile   = GH_.DetName_ + "_Mesh.mphtxt";
        std::string fieldfile  = GH_.DetName_ + "_Data.txt";
        std::string fileconfig = GH_.DetName_ + "MaterialProperties.txt";

        // Setup the electric potential map
        Garfield::ComponentComsol* fm = new Garfield::ComponentComsol(); // Field Map
        fm->Initialise(home + meshfile ,home + fileconfig, home + fieldfile, "cm");
        
        // Print some information about the cell dimensions.
        fm->PrintRange();

        // Associate the gas with the corresponding field map material.
        fm->SetGas(fMediumMagboltz); 
        fm->PrintMaterials();
        fm->Check();

        fSensor->AddComponent(fm);
        // fSensor->SetArea(-DetChamberR, -DetChamberR, -DetChamberL/2.0, DetChamberR, DetChamberR, DetChamberL/2.0); // cm

    }

        
    fAvalancheMC = new Garfield::AvalancheMC(); // drift, not avalanche, to be fair.
    fAvalancheMC->SetSensor(fSensor);
    fAvalancheMC->SetTimeSteps(0.05);      // nsec, per example
    fAvalancheMC->SetDistanceSteps(2.e-2); // cm, 10x example
    fAvalancheMC->EnableDebugging(false);  // way too much information. 
    fAvalancheMC->DisableAttachment();     // Currently getting warning messages about the attachment. You can supress those by switching this on.
    fAvalancheMC->EnableDriftLines();

    // Load in the events
    if (GH_.useELFile_)
        GetTimeProfileData(nexus_path + "/data/Garfield/" + GH_.DetName_ + "_Profiles_Rotated.csv", EL_profiles, EL_events);
    
}

void GarfieldVUVPhotonModel::Reset(){

  std::cout <<"Resetting GarfieldVUV" << std::endl;
  fSensor->ClearSignal();

  track_ids_.clear();
  N_ioni_.clear();
  N_S2_.clear();
}


Garfield::ComponentUser* GarfieldVUVPhotonModel::SetComponentField(){

    //  ---- Create the Garfield Field region --- 
    Garfield::GeometrySimple* geo = new Garfield::GeometrySimple();

    // Tube oriented in z-axis (0.,0.,1.,) The addition of the 1 cm is for making sure it doesnt fail on the boundary
    Garfield::SolidTube* tube = new Garfield::SolidTube(GH_.origin_.x(), GH_.origin_.y() ,GH_.origin_.z(), GH_.DetChamberR_+1, GH_.DetChamberL_, 0.,0.,1.);

    // Add the solid to the geometry, together with the medium inside
    geo->AddSolid(tube, fMediumMagboltz);

    // Make a component with analytic electric field
    Garfield::ComponentUser* componentDriftLEM = new Garfield::ComponentUser();
    componentDriftLEM->SetGeometry(geo);

    // Set the electric field calculation function within the custom bounds for the component
    componentDriftLEM->SetElectricField([&](double x, double y, double z, double& ex, double& ey, double& ez) {
    
    // Only want Ez component to the field
    ex = ey = 0.;

    // Define a field region for the whole gas region

    // Set ez for regions outside of the radius of the FC
    if ( std::sqrt(x*x + y*y) > GH_.DetActiveR_){
        ez = -GH_.fieldDrift_; // Negative field will send them away from the LEM region
    }

    // Field past the cathode drift them away from the LEM with negative field
    // These make no photons and we kill these anyway
    // Garfield does not like any points with zero field
    if (z > GH_.CathodePos_)
        ez = -GH_.fieldDrift_;

    // Drift region
    if (z <= GH_.CathodePos_)
        ez = GH_.fieldDrift_;

    // EL region
    if (z <= GH_.ELPos_ && z > GH_.ELPos_-GH_.gap_EL_)
        ez = GH_.fieldEL_;

    // Drift towards the end cap
    if (z <= GH_.ELPos_ - GH_.gap_EL_)
        ez = GH_.fieldDrift_; 
  });

  // Printing pressure and temperature
  std::cout << "GarfieldVUVPhotonModel::buildBox(): Garfield mass density [g/cm3], pressure [Torr], temp [K]: " <<
        geo->GetMedium(0.,0.,0.)->GetMassDensity() << ", " << geo->GetMedium(0.,0.,0.)->GetPressure()<< ", " 
        << geo->GetMedium(0.,0.,0.)->GetTemperature() << std::endl;

  return componentDriftLEM;

}


void GarfieldVUVPhotonModel::MakeELPhotonsFromFile( G4FastStep& fastStep, G4double xi, G4double yi, G4double zi, G4double ti, G4int trk_id){

    // Here we get the photon timing profile from a file
    G4int EL_event  = round(G4UniformRand()* (EL_events.size() - 1) );
    
    // Get the right EL array
    std::vector<std::vector<G4double>> EL_profile = EL_profiles[EL_event];

    // Now loop over and make the photons
    G4int el_gain = EL_profile.size();
    // el_gain=1 ;

    G4double tig4(0.);
    
    for (G4int i=0;i<el_gain;i++){
      
      // std::cout << xi << ", " <<  EL_profile[i][0]  << std::endl;
      G4ThreeVector fakepos ( (xi+ EL_profile[i][0])*cm, (yi+ EL_profile[i][1])*cm, (zi+ EL_profile[i][2])*cm ); // 0 = x, 1 = y, 2 = z. Garfield units are cm, x10 for mm G4 units
    
      auto* optphot = G4OpticalPhoton::OpticalPhotonDefinition();

      G4ThreeVector momentum, polarization;
      GetPhotonPol(momentum, polarization);

      // Get Photon energy
      G4double sc_max = spectrum_integral->GetMaxValue();
      G4double sc_value = G4UniformRand()*sc_max;
      G4double sampled_energy = spectrum_integral->GetEnergy(sc_value);

      
      G4DynamicParticle VUVphoton(optphot, momentum, sampled_energy);

      G4double drift_time = EL_profile[i][3];
      G4bool el_survived = GetAttachment(GH_.gap_EL_*cm/GH_.v_drift_el_);

      // The electron did not survive drifting in the EL gap due to attachment
      if (!el_survived)
        break;
      
      tig4 = ti + drift_time + GetScintTime(); // in nsec, t = index 3 in vector. Units are ns, so just add it on
      // std::cout <<  "fakepos,time is " << fakepos[0] << ", " << fakepos[1] << ", " << fakepos[2] << ", " << tig4 << std::endl;
      
      G4Track *newTrack=fastStep.CreateSecondaryTrack(VUVphoton, fakepos, tig4 ,false);
      newTrack->SetPolarization(polarization);
    
      N_S2_[GetCurrentTrackIndex(trk_id)]++;
    }
}


void GarfieldVUVPhotonModel::MakeELPhotonsSimple(G4FastStep& fastStep, G4double xi, G4double yi, G4double zi, G4double ti, G4int trk_id){

    // std::cout << "Generating Photons"<< std::endl;
    // EL field in G4 units (V/mm), Pressure in G4 units (MeV/mm3), gap in G4 units (mm)
    const G4int el_gain = XenonELLightYield(GH_.fieldEL_*kilovolt/cm, GH_.GasPressure_)*GH_.gap_EL_*cm; // E [V/mm], P [MeV/mm3], EL gap [mm]

    // std::cout << " Yield is "<< el_gain <<" Field " <<GH_.fieldEL_<< " Pressure  " << GH_.GasPressure_/bar<< " EL  " << GH_.gap_EL_*cm << std::endl;
    
    G4double tig4(0.);

    for (G4int i=0;i<el_gain;i++){

      G4ThreeVector fakepos (xi*cm,yi*cm,zi*cm); /// ignoring diffusion in small LEM gap, EC 17-June-2022.     
      
      auto* optphot = G4OpticalPhoton::OpticalPhotonDefinition();

      G4ThreeVector momentum, polarization;
      GetPhotonPol(momentum, polarization);

      // Get Photon energy
      G4double sc_max = spectrum_integral->GetMaxValue();
      G4double sc_value = G4UniformRand()*sc_max;
      G4double sampled_energy = spectrum_integral->GetEnergy(sc_value);
      
      G4DynamicParticle VUVphoton(optphot, momentum, sampled_energy);
    
      // in nsec (gap_EL_ is in cm). Still ignoring diffusion in small LEM.
      // Also add in the xenon scintillation timing delays and attachment
      G4double drift_time = G4float(i)/G4float(el_gain)*GH_.gap_EL_*cm/GH_.v_drift_el_;
      G4bool el_survived = GetAttachment(GH_.gap_EL_*cm/GH_.v_drift_el_);

      // The electron did not survive drifting in the EL gap
      if (!el_survived)
        break;

      tig4 = ti + drift_time + GetScintTime(); 

      // std::cout <<  "pos,time is " << fakepos[0] << ", " << fakepos[1] << ", " << fakepos[2] << ", " << tig4 << std::endl;
      
      G4Track *newTrack=fastStep.CreateSecondaryTrack(VUVphoton, fakepos, tig4 ,false);
      newTrack->SetPolarization(polarization);
      
      N_S2_[GetCurrentTrackIndex(trk_id)]++;
    }


}

void GarfieldVUVPhotonModel::PrintElectricField(G4double x,G4double y, G4double z){
    std::array<G4double, 3> ef{0,0,0};
    std::array<G4double, 3> bf{0,0,0};
    std::vector<G4double> vf{0,0,0};
    Garfield::Medium* medium = nullptr;
    G4int status(0);
    fSensor->ElectricField(x, y, z, ef[0], ef[1], ef[2], medium, status);
    std::cout << "GVUVPM: E field in medium " << medium << " at " << x<<","<<y<<","<<z << " is: " << ef[0]<<","<<ef[1]<<","<<ef[2] << std::endl;
}

void GarfieldVUVPhotonModel::InsertHits(G4double x,G4double y, G4double z, G4double t, G4int trk_id, G4double mean_ioni_E){
  G4ThreeVector ie_pos(x*cm,y*cm,z*cm);
  IonizationHit* ie_hit = new IonizationHit();
  ie_hit->SetTrackID(trk_id);
  ie_hit->SetTime(t);
  ie_hit->SetEnergyDeposit(mean_ioni_E);
  ie_hit->SetPosition(ie_pos);
  fGarfieldSD->InsertIonizationHit(ie_hit);

}

G4double GarfieldVUVPhotonModel::GetScintTime(){

  G4double scint_time = 0;

  // Define decay constants and probabilities
  const G4double decayConstant1  = 4.5;   // ns
  const G4double probability1    = 0.1;   // 10%
  const G4double decayConstant2  = 100.0; // ns
  const G4double probability2    = 0.9;   // 90%

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

G4bool GarfieldVUVPhotonModel::GetAttachment(G4double t){

  G4double lifetime = GH_.e_lifetime_;

  G4double survivalProb = exp(-t / lifetime);

  // Generate a random number
  G4double randNum = G4UniformRand();
  
  // Determine if survival occurred based on the random number
  return randNum < survivalProb;
}

G4bool GarfieldVUVPhotonModel::CheckXYBoundsPolygon(std::vector<G4double> point){

  // This function checks if a point is inside a polygon (needed for checking if
    // inside the drift tube which has an octagonal-like shape).
    // It works by counting intersections of the point in the horizontal plane
    // as described in 
    // https://www.geeksforgeeks.org/how-to-check-if-a-given-point-lies-inside-a-polygon/
    

    G4double x = point[0], y = point[1];
    G4bool inside = false;

    // Store the first point in the polygon and initialize the second point
    std::vector<G4double> p1 = {polygon_[0][0],polygon_[0][1]};

    // Loop through each edge in the polygon
    for (G4int i = 1; i <= GH_.nsides_; i++) {
      
      // Get the next point in the polygon
      std::vector<G4double> p2 ={polygon_[i % GH_.nsides_][0], polygon_[i % GH_.nsides_][1]};

      // 1. Check if the point is above the minimum y coordinate of the edge
      // 2. Check if the point is below the maximum y coordinate of the edge
      // 3. Check if the point is to the left of the maximum x coordinate of the edge
      if ( (y > std::min(p1[1], p2[1])) && (y <= std::max(p1[1], p2[1])) && (x <= std::max(p1[0], p2[0])) ) {

        // Calculate the x-intersection of the line connecting the point to the edge
        G4double x_intersection = (y - p1[1]) * (p2[0] - p1[0]) / (p2[1] - p1[1]) + p1[0];

        // Check if the point is on the same line as the edge or to the left of
        // the x-intersection
        if (p1[0] == p2[0] || x <= x_intersection) {
          // Flip the inside flag
          inside = !inside;
        }
      }

        // Store the current point as the first point for the next iteration
        p1 = p2;
    }
   
    // Return the value of the inside flag
    return inside;
}

void GarfieldVUVPhotonModel::InitalizePolygon(){

  std::cout << "Initializing polygon with " << GH_.nsides_ << " sides" << std::endl;

  // Define the number of sides, the radius, and the center 
  G4int n = GH_.nsides_;
  G4double radius = GH_.DetActiveR_;

  // Plot the regular polygon
  std::vector<G4double> theta(n);
  for (G4int i = 0; i < n; ++i) {
      theta[i] = 2 * M_PI * i / n;
  }

  for (int i = 0; i < n; ++i) {
      polygon_.push_back({radius * cos(theta[i]) + GH_.origin_.x(), radius * sin(theta[i]) + GH_.origin_.y()});
  }

}

void GarfieldVUVPhotonModel::AddTrack(G4int trk_id){

    // Check if track exists in vec, if not then add
    if (std::find(track_ids_.begin(), track_ids_.end(), trk_id) == track_ids_.end()) {
        std::cout <<"Adding Track!" << std::endl;
        track_ids_.push_back(trk_id);
        N_ioni_.push_back(0);
        N_S2_.push_back(0);
    }
}

G4int GarfieldVUVPhotonModel::GetCurrentTrackIndex(G4int trk_id){

    G4int index =-1;

    // Search for searchValue in vec
    auto it = std::find(track_ids_.begin(), track_ids_.end(), trk_id);

    // Check if value was found
    if (it != track_ids_.end()) {
        // Calculate the index
        index = std::distance(track_ids_.begin(), it);
    }

    return index;

}


}