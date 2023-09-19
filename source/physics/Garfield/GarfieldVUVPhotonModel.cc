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
#include "G4OpticalPhoton.hh"
#include "G4ProcessManager.hh"
#include "G4EventManager.hh"
#include "Garfield/ComponentComsol.hh"
#include "S2Photon.h"
#include "config.h"
#include "G4AnalysisManager.hh"
#include "G4Event.hh"

namespace nexus{

const G4double res(0.01); // Estimated fluctuations in EL yield - high? EC, 21-June-2022.

GarfieldVUVPhotonModel::GarfieldVUVPhotonModel(G4String modelName,G4Region* envelope, GarfieldHelper GH) : G4VFastSimulationModel(modelName, envelope) {
    
    GH_ = GH;
    GH_.DumpParams();
    InitialisePhysics();

    G4OpBoundaryProcess* fBoundaryProcess   = new G4OpBoundaryProcess();
    G4OpAbsorption*      fAbsorptionProcess = new G4OpAbsorption();
    G4OpWLS*             fTheWLSProcess     = new G4OpWLS();

   
}

G4bool GarfieldVUVPhotonModel::IsApplicable(const G4ParticleDefinition& particleType) {
  //  std::cout << "GarfieldVUVPhotonModel::IsApplicable() particleType is " << particleType.GetParticleName() << std::endl;
  if (particleType.GetParticleName()=="thermalelectron") // || particleType.GetParticleName()=="opticalphoton")
    return true;
  return false;
}

G4bool GarfieldVUVPhotonModel::ModelTrigger(const G4FastTrack& fastTrack){
 
  G4double ekin = fastTrack.GetPrimaryTrack()->GetKineticEnergy();
  //  std::cout << "GarfieldVUVPhotonModel::ModelTrigger() thermalE, ekin is " << GH_.thermalE_ << ",  "<< ekin / MeV << std::endl;
  
  G4String particleName = fastTrack.GetPrimaryTrack()->GetParticleDefinition()->GetParticleName();

   if (ekin < GH_.thermalE_ && particleName=="thermalelectron")
    {
      return true;
    }
  return false;

} 
    
void GarfieldVUVPhotonModel::DoIt(const G4FastTrack& fastTrack, G4FastStep& fastStep) 
{

  /* 
     This tracks all the ionization/conversion electrons created by Degrad in its simulation of the primary gamma's photoelectric effect on Xe.
     Each such electron is drifted in the E-field and avalanched, as appropriate. That creates excited Xe atoms. We put one ELM photon
     per excitation of 172 (7.2) nm (eV) onto the optical photon stack. Geant4 will track those in the normal way. 


     Note, the weirdness of userHandle, which seems to be called at each excitation by the AvalancheMicroscopic model and fills our
     GarfieldExcitationHitCollection. And which somehow we're allowed to grab here.

     EC, 2-Dec-2021.

   */
    fastStep.SetNumberOfSecondaryTracks(1E3);
  
    // G4cout<<"HELLO Garfield"<<G4endl;
    ////The details of the Garfield model are implemented here
     fastStep.KillPrimaryTrack();//KILL NEST/DEGRAD/G4 TRACKS
    //  fastStep.ProposeTrackStatus(fSuspend);

     garfPos = fastTrack.GetPrimaryTrack()->GetVertexPosition();
     garfTime = fastTrack.GetPrimaryTrack()->GetGlobalTime();
     //G4cout<<"GLOBAL TIME "<<G4BestUnit(garfTime,"Time")<<" POSITION "<<G4BestUnit(garfPos,"Length")<<G4endl;
     counter[1]++; // maybe not threadsafe
     if (!(counter[1]%1000))
       G4cout << "GarfieldVUV: actual NEST thermales: " << counter[1] << G4endl;

      if (!(counter[3]%1000) and (counter[3]>0))
         G4cout << "GarfieldVUV: S2 OpticalPhotons: " << counter[3] << G4endl;


    //     if (!(counter[1]%1000)) // uncomment!
       GenerateVUVPhotons(fastTrack,fastStep,garfPos,garfTime);

}

// GarfieldExcitationHitsCollection *garfExcHitsCol;

void GarfieldVUVPhotonModel::GenerateVUVPhotons(const G4FastTrack& fastTrack, G4FastStep& fastStep,G4ThreeVector garfPos,G4double garfTime)
{

     // Drift in main region, then as LEM region is approached and traversed, avalanche multiplication and excitations  will occur.
  
    G4double x0=garfPos.getX()*0.1;//Garfield length units are in cm
    G4double y0=garfPos.getY()*0.1;
    G4double z0=garfPos.getZ()*0.1;
    G4double t0=garfTime;
    

    G4double ekin = fastTrack.GetPrimaryTrack()->GetKineticEnergy();
    G4ThreeVector dir = fastTrack.GetPrimaryTrack()->GetMomentumDirection();
    G4String particleName = fastTrack.GetPrimaryTrack()->GetParticleDefinition()->GetParticleName();

    // std::cout << "GVUVPM: start positions are " << particleName << " " << x0<<"," <<y0<<","<<z0 <<"," <<t0<< std::endl;

    int status(0);

    // Debug the electric field
    // std::array<double, 3> ef{0,0,0};
    // std::array<double, 3> bf{0,0,0};
    // std::vector<double> vf{0,0,0};
    // Garfield::Medium* medium = nullptr;
    // fSensor->ElectricField(x0,y0,-12, ef[0], ef[1], ef[2], medium, status);                                        
    // std::cout << "GVUVPM: E field in medium " << medium << " at " << x0<<","<<y0<<","<<z0 << " is: " << ef[0]<<","<<ef[1]<<","<<ef[2] << std::endl;
    double xi,yi,zi,ti;
    // Need to get the AvalancheMC drift at the High-Field point in z, and then call fAvalanche-AvalancheElectron() to create excitations/VUVphotons.
    fAvalancheMC->DriftElectron(x0,y0,z0,t0);

    size_t n =fAvalancheMC->GetElectrons().at(0).path.size();
    auto GetDriftLines=fAvalancheMC->GetElectrons().at(0).path;

    // std::cout << "Avalanche end points: " << n<< std::endl;
    // Get zi when in the beginning of the EL region
    for(unsigned int i=0;i<n;i++){

      xi=GetDriftLines.at(i).x;
      yi=GetDriftLines.at(i).y;
      zi=GetDriftLines.at(i).z;
      ti=GetDriftLines.at(i).t;
      // std::cout << "GVUVPM: positions are " << xi<<"," <<yi<<","<<zi <<"," <<ti<< std::endl;
      
      // Drift line point entered LEM
      if (zi < ELPos_ && ( std::sqrt(xi*xi + yi*yi) < GH_.DetActiveR_) )
        break; 
      
      // No drift line point meets criteria, so return
      else if (i==n-1)
        return;

    }  // pts in driftline

    // Generate the El photons from a microphys model ran externally in Garfield
    // We sample the output file which contains the timing profile of emission and diffusion
    G4bool use_ELFile = false;
    if (use_ELFile)
        MakeELPhotonsFromFile(fastStep, xi, yi, zi, ti);
    // Use a simpler model
    else
        MakeELPhotonsSimple(fastStep, xi, yi, zi, ti);


    // std::cout << "GVUVPM: Avalanching in high field starting at: "  << xi<<"," <<yi<<","<<zi <<"," <<ti << std::endl;

        
}


// Selection of Xenon exitations and ionizations
void GarfieldVUVPhotonModel::InitialisePhysics(){
    // Set the gas Properties
    fMediumMagboltz = new Garfield::MediumMagboltz();
    fMediumMagboltz->SetComposition("Xe", 100.);
    fMediumMagboltz->DisableDebugging();
    
    //  --- Load in Ion Mobility file --- 
    ionMobFile = "IonMobility_Ar+_Ar.txt";
    const std::string path = getenv("GARFIELD_HOME");
    
    if(ionMobFile!="")
      fMediumMagboltz->LoadIonMobility(path + "/Data/" + ionMobFile);
    
    if(gasFile!=""){
      
      std::cout << "Loaded gasfile." << std::endl;
    }

    gasFile = "data/Xenon_10Bar.gas";
    G4cout << gasFile << G4endl;
    fMediumMagboltz->LoadGasFile(gasFile.c_str());
    std::cout << "Finished Loading in the gas file" << std::endl;

    // Initialize the gas
    fMediumMagboltz->Initialise(true);

    // Print the gas properties
    // fMediumMagboltz->PrintGas();

    ELPos_ =  - GH_.DetActiveL_/2.0;
    FCTop_ =  + GH_.DetActiveL_/2.0;

    std::cout << "Detector Dimentions: "<< GH_.DetChamberR_ << " " << GH_.DetChamberL_ << "  " << GH_.DetActiveR_ << "  " << GH_.DetActiveL_ << std::endl; 

    fSensor = new Garfield::Sensor();

    G4bool use_COMSOL = false;

    if (!use_COMSOL){
        Garfield::ComponentUser* componentDriftLEM = CreateSimpleGeometry();
        fSensor->AddComponent(componentDriftLEM);
        
        // Set the region where the sensor is active -- based on the gas volume
        fSensor->SetArea(-GH_.DetChamberR_, -GH_.DetChamberR_, -GH_.DetChamberL_/2.0, GH_.DetChamberR_, GH_.DetChamberR_, GH_.DetChamberL_/2.0); // cm

    }
    else {
        std::cout << "Initialising Garfiled with a COMSOL geometry" << std::endl;

        G4String home = "./data/Garfield/";
        std::string meshfile   = "CRAB_Mesh.mphtxt";
        std::string fieldfile  = "CRAB_Data.txt";
        std::string fileconfig = "CRABMaterialProperties.txt";

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

    G4bool use_ELFile = false;

    // Load in the events
    if (use_ELFile)
        GetTimeProfileData("data/Garfield/CRAB_Profiles_Rotated.csv", EL_profiles, EL_events);
    
}

void GarfieldVUVPhotonModel::Reset()
{
  fSensor->ClearSignal();
  counter[1] = 0;
  counter[3] = 0;
}


Garfield::ComponentUser* GarfieldVUVPhotonModel::CreateSimpleGeometry(){

    //  ---- Create the Garfield Field region --- 
    Garfield::GeometrySimple* geo = new Garfield::GeometrySimple();

    // Tube oriented in Y'axis (0.,1.,0.,) The addition of the 1 cm is for making sure it doesnt fail on the boundary
    Garfield::SolidTube* tube = new Garfield::SolidTube(0.0, 0.0 ,0.0, GH_.DetChamberR_+1, GH_.DetChamberL_*0.5, 0.,0.,1.);

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
      if ( std::sqrt(x*x + y*y) > GH_.DetActiveR_/2.0){
          ez = -GH_.fieldDrift_; // Negative field will send them away from the LEM region
      }

      // Field past the cathode drift them away from the LEM with negative field
      if (z > FCTop_)
          ez = -GH_.fieldDrift_;

      // Drift region
      if (z <= FCTop_)
          ez = GH_.fieldDrift_;

      // EL region
      if (z <= ELPos_ && z > ELPos_-GH_.gap_EL_)
          ez = GH_.fieldEL_;

      // Drift towards the end cap
      if (z <= ELPos_ - GH_.gap_EL_)
          ez = GH_.fieldDrift_; 
    });

    // Printing pressure and temperature
    std::cout << "GarfieldVUVPhotonModel::buildBox(): Garfield mass density [g/cm3], pressure [Torr], temp [K]: " <<
         geo->GetMedium(0.,0.,0.)->GetMassDensity() << ", " << geo->GetMedium(0.,0.,0.)->GetPressure()<< ", " 
         << geo->GetMedium(0.,0.,0.)->GetTemperature() << std::endl;

    return componentDriftLEM;

}


void GarfieldVUVPhotonModel::MakeELPhotonsFromFile( G4FastStep& fastStep, G4double xi, G4double yi, G4double zi, G4double ti){

    // Here we get the photon timing profile from a file
    G4int EL_event  = round(G4UniformRand()* (EL_events.size() - 1) );
    
    // Get the right EL array
    std::vector<std::vector<G4double>> EL_profile = EL_profiles[EL_event];

    // Now loop over and make the photons
    G4int colHitsEntries = EL_profile.size();
    // std::cout <<  colHitsEntries<< std::endl;
    // colHitsEntries=1 ;

    G4double tig4(0.);
    
    for (G4int i=0;i<colHitsEntries;i++){
      
        // std::cout << xi << ", " <<  EL_profile[i][0]  << std::endl;
      G4ThreeVector fakepos ( (xi+ EL_profile[i][0])*10., (yi+ EL_profile[i][1])*10., (zi+ EL_profile[i][2])*10. ); // 0 = x, 1 = y, 2 = z
      
      if (i % (colHitsEntries/colHitsEntries ) == 0){ // 50. Need to uncomment this condition, along with one in degradmodel.cc. EC, 2-Dec-2021.
      
        auto* optphot = S2Photon::OpticalPhotonDefinition();
        
        G4DynamicParticle VUVphoton(optphot,G4RandomDirection(), 7.2*eV);
       
        tig4 = ti + EL_profile[i][3]; // in nsec, t = index 3 in vector. Units are ns, so just add it on
        // std::cout <<  "fakepos,time is " << fakepos[0] << ", " << fakepos[1] << ", " << fakepos[2] << ", " << tig4 << std::endl;
        
        G4Track *newTrack=fastStep.CreateSecondaryTrack(VUVphoton, fakepos, tig4 ,false);
        newTrack->SetPolarization(G4ThreeVector(0.,0.,1.0)); // Needs some pol'n, else we will only ever reflect at an OpBoundary. EC, 8-Aug-2022.
        //	G4ProcessManager* pm= newTrack->GetDefinition()->GetProcessManager();
        //	G4ProcessVectorfAtRestDoItVector = pm->GetAtRestProcessVector(typeDoIt);
      }
      counter[3]++;
    }
}


void GarfieldVUVPhotonModel::MakeELPhotonsSimple(G4FastStep& fastStep, G4double xi, G4double yi, G4double zi, G4double ti){

  //std::cout << "Generating Photons"<< std::endl;
    
    G4int colHitsEntries= 0.0; //garfExcHitsCol->entries();
    //	G4cout<<"GarfExcHits entries "<<colHitsEntries<<G4endl; // This one is not cumulative.

    const G4double YoverP = 140.*GH_.fieldEL_/((GH_.GasPressure_/bar)*1000) - 116.; // yield/cm/bar, with P in Torr ... JINST 2 p05001 (2007).
    colHitsEntries = YoverP * (GH_.GasPressure_/bar) * (GH_.gap_EL_/10); // with P in bar this time.
    // colHitsEntries*=2; // Max val before G4 cant handle the memory anymore
    //std::cout<<" Yield is "<<colHitsEntries <<" Field " <<GH_.fieldEL_<< " Pressure  " << GH_.GasPressure_/bar<< " EL  " << GH_.gap_EL_/10<<std::endl;
     //colHitsEntries=1; // This is to turn down S2 so the vis doesnt get overwelmed

    colHitsEntries *= (G4RandGauss::shoot(1.0,res));
    
    G4double tig4(0.);
    const G4double vd(2.4); // mm/musec, https://arxiv.org/pdf/1902.05544.pdf. Pretty much flat at our E/p..
    for (G4int i=0;i<colHitsEntries;i++){

      G4ThreeVector fakepos (xi*10,yi*10.,zi*10.-10*GH_.gap_EL_*float(i)/float(colHitsEntries)); /// ignoring diffusion in small LEM gap, EC 17-June-2022.     
      
      if (i % (colHitsEntries/colHitsEntries ) == 0){ // 50. Need to uncomment this condition, along with one in degradmodel.cc. EC, 2-Dec-2021.
      
        auto* optphot = S2Photon::OpticalPhotonDefinition();
        
        G4DynamicParticle VUVphoton(optphot,G4RandomDirection(), 7.2*eV);
       

        /// std::cout <<  "fakepos,time is " << fakepos[0] << ", " << fakepos[1] << ", " << fakepos[2] << ", " << ti << std::endl;
       

        tig4 = ti + float(i)/float(colHitsEntries)*GH_.gap_EL_*10./vd*1E3; // in nsec (gap_EL_ is in cm). Still ignoring diffusion in small LEM.
        //	    std::cout << "VUV photon time [nsec]: " << tig4 << std::endl;
        
        G4Track *newTrack=fastStep.CreateSecondaryTrack(VUVphoton, fakepos, tig4 ,false);
        newTrack->SetPolarization(G4ThreeVector(0.,0.,1.0)); // Needs some pol'n, else we will only ever reflect at an OpBoundary. EC, 8-Aug-2022.
        //	G4ProcessManager* pm= newTrack->GetDefinition()->GetProcessManager();
        //	G4ProcessVectorfAtRestDoItVector = pm->GetAtRestProcessVector(typeDoIt);
      }
      counter[3]++;
    }


}
}