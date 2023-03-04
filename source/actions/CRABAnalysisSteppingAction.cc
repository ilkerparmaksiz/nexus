//
// Created by ilker on 6/28/22.
//

#include "CRABAnalysisSteppingAction.h"
#include "FactoryBase.h"

#include <G4Step.hh>
#include <G4SteppingManager.hh>
#include <G4ProcessManager.hh>
#include <G4OpticalPhoton.hh>
#include <G4OpBoundaryProcess.hh>
#include <G4VPhysicalVolume.hh>
#include "FileHandling.h"
#include "IonizationElectron.h"

using namespace nexus;

REGISTER_CLASS(CRABAnalysisSteppingAction, G4UserSteppingAction)

CRABAnalysisSteppingAction::CRABAnalysisSteppingAction(): G4UserSteppingAction(),msg_(0),filePath_(""),SavetoFile_(false),NumEvents_(0),FileName_("photoncount.txt"),isLead210Decay_(false)
{
    msg_=new G4GenericMessenger(this,"/Actions/CRABAnalysisSteppingAction/");
    msg_->DeclareProperty("FilePath",filePath_,"This is the path for saving some counts to text file..");
    msg_->DeclareProperty("FileSave",SavetoFile_,"Save Counts to File");
    msg_->DeclareProperty("FileName",FileName_,"FileName To save count information");
    msg_->DeclareProperty("isLead210",isLead210Decay_,"Extra Analysis for Lead210 Decay");

    Names={"Pb210","Bi210[46.539]","Bi210","Po210"};
    TotalPhotons=0;
    TotalIonizationElectron=0;

}



CRABAnalysisSteppingAction::~CRABAnalysisSteppingAction()
{
    FileHandling *f1=new FileHandling();
    G4int total_counts = 0;

    G4String str;
    G4String Path;
    if(filePath_.empty())
        Path=FileName_;
    else Path=filePath_+"/"+FileName_;
    if(!isLead210Decay_){
        detectorCounts::iterator it = my_counts_.begin();
        DepositedEnergy ::iterator depit=PhotonEnergy_.begin();
        G4cout<<"#### Detected Photons ####"<<G4endl;

        while (it != my_counts_.end()) {
            G4cout << "Detector " << it->first << ": " << it->second << " counts" << G4endl;
            if(SavetoFile_){

                str=it->first + ",None," + to_string(it->second) ;
                f1->SaveToTextFile(Path,"Detector,Source,Counts",str);
            }
            if(DisplayPhotonEnergy_){
                for(int i=0;i<PhotonEnergy_[it->first].size();i++){
                    G4cout << "Detector " << it->first << ": " << PhotonEnergy_[it->first][i] << " Energy" << G4endl;
                }
            }

            total_counts += it->second;
            it ++;
        }


    } else{
        G4cout<<"#### Detailed info about  Photons ####"<<G4endl;

        DetailedPhotonCounts ::iterator DDetectors=DPhotonCounts.begin();
        while (DDetectors!=DPhotonCounts.end()){
            detectorCounts ::iterator DPhotons=DDetectors->second.begin();
            while(DPhotons!=DDetectors->second.end()){
                G4cout<< DDetectors->first<<","<<DPhotons->first<<" :"<<DPhotons->second<<G4endl;
                if(SavetoFile_){

                    str=DDetectors->first + "," +DPhotons->first +","+to_string(DPhotons->second) ;
                    f1->SaveToTextFile(Path,"Detector,Source,Counts",str);
                }
                total_counts=total_counts+DPhotons->second;

                DPhotons++;
            }

            DDetectors++;
        }
    }




    G4int TotalAbsPhotons=0;
    if(filePath_.empty())
        Path="Extra_"+FileName_;
    else Path=filePath_+"/"+"Extra_"+FileName_;
    // For Absorbed Photons
    if(!isLead210Decay_){

        detectorCounts::iterator Absit = ObservedPhotons.begin();

        G4cout <<"------------ Absorbed Photons ----"<<G4endl;

        while (Absit != ObservedPhotons.end()) {
            G4cout << "Detector " << Absit->first << ": " << Absit->second << " counts" << G4endl;
            if(SavetoFile_){
                str=Absit->first +  ",None," + to_string(Absit->second) ;
                f1->SaveToTextFile(Path,"Detector,Source,Counts",str);
            }
            TotalAbsPhotons += Absit->second;
            Absit ++;
        }
    }else{
        G4cout<<"#### Absorbed Photons ####"<<G4endl;

        DetailedPhotonCounts ::iterator DAbsDetectors=DAbsPhotonCounts.begin();
        while (DAbsDetectors!=DAbsPhotonCounts.end()){
            detectorCounts ::iterator DAbsPhotons=DAbsDetectors->second.begin();
            while(DAbsPhotons!=DAbsDetectors->second.end()){
                G4cout<<DAbsDetectors->first<<","<<DAbsPhotons->first<<" :"<<DAbsPhotons->second<<G4endl;
                if(SavetoFile_){
                    str=DAbsDetectors->first + "," +DAbsPhotons->first +","+to_string(DAbsPhotons->second) ;
                    f1->SaveToTextFile(Path,"Detector,Source,Counts",str);
                }
                TotalAbsPhotons=TotalAbsPhotons+DAbsPhotons->second;
                DAbsPhotons++;
            }

            DAbsDetectors++;
        }
        // Total Produced Photons
        ProducedPhotonCount ::iterator ProducedPhotons=DTotalPhotonCounts.begin();
        while (ProducedPhotons!=DTotalPhotonCounts.end()){
            G4cout<<ProducedPhotons->first<<","<<ProducedPhotons->first<<" :"<<ProducedPhotons->second<<G4endl;
            if(SavetoFile_){
                str="All_p," +ProducedPhotons->first +","+to_string(ProducedPhotons->second) ;
                f1->SaveToTextFile(Path,"Detector,Source,Counts",str);
            }
            ProducedPhotons++;
        }
        // Total Produced Electrons

        G4cout<<"#### Electrons ####"<<G4endl;
        DetailedElectronCounts::iterator Die=DElectronCounts.begin();
        while (Die!=DElectronCounts.end()){
            G4cout<<Die->first<<" :"<<Die->second<<G4endl;
            if(SavetoFile_){
                str="All_e," +Die->first +","+to_string(Die->second) ;
                f1->SaveToTextFile(Path,"Detector,Source,Counts",str);
            }
            Die++;
        }
    }



    if(SavetoFile_){
        str="Total_Detected,None," + to_string(total_counts);
        f1->SaveToTextFile(Path,"Detector,Source,Counts",str);
        str="Absorption_Photons,None,"+ to_string(TotalAbsPhotons);
        f1->SaveToTextFile(Path,"Detector,Source,Counts",str);
        str="Total_Produced_Photons,None," + to_string(TotalPhotons);
        f1->SaveToTextFile(Path,"Detector,Source,Counts",str);
        str="Total_Produced_ie,None," + to_string(TotalIonizationElectron);
        f1->SaveToTextFile(Path,"Detector,Source,Counts",str);
    }

    G4cout << "Abosrved Photons " <<TotalAbsPhotons<<G4endl;
    G4cout << "Total_Detected: " << total_counts << G4endl;
    G4cout << "Total Produced Photons " <<TotalPhotons<<G4endl;
    G4cout << "Total Produced iElectrons " <<TotalIonizationElectron<<G4endl;





}



void CRABAnalysisSteppingAction::UserSteppingAction(const G4Step* step)
{
    G4ParticleDefinition* pdef = step->GetTrack()->GetDefinition();


    if(std::find(Names.begin(),Names.end(),pdef->GetParticleName())!=Names.end() ){
        TempTrackID=step->GetTrack()->GetTrackID();
        TempName=pdef->GetParticleName();
        return;
    }

    if(pdef == IonizationElectron::Definition()){
        DetailedElectronCounts::iterator Die=DElectronCounts.find(TempName);
        G4Track* track = step->GetTrack();

        G4String volumePost = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName();

        if(volumePost=="EL_GAP"){
            TotalIonizationElectron++;
            Die!=DElectronCounts.end() ? DElectronCounts[Die->first]+=1:DElectronCounts[TempName]+=1;
            // G4cout<<"Time -> "<<track->GetLocalTime()<<G4endl;
        }
        if (volumePost!="FIELDCAGE" and volumePost!="EL_GAP"){
            track->SetTrackStatus(fStopAndKill);
            TotalIonizationElectron++;
            return;
        }
        return;
    }


    //Check whether the track is an optical photon
    if (pdef != G4OpticalPhoton::Definition()) return;

    /*
    // example of information one can access about optical photons

    G4Track* track = step->GetTrack();
    G4int pid = track->GetParentID();
    G4int tid = track->GetTrackID();
    G4StepPoint* point1 = step->GetPreStepPoint();
    G4StepPoint* point2 = step->GetPostStepPoint();
    G4TouchableHandle touch1 = point1->GetTouchableHandle();
    G4TouchableHandle touch2 = point2->GetTouchableHandle();
    G4String vol1name = touch1->GetVolume()->GetName();
    G4String vol2name = touch2->GetVolume()->GetName();

    G4String proc_name = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    G4int copy_no = step->GetPostStepPoint()->GetTouchable()->GetReplicaNumber(1);
    */
    G4Track* track = step->GetTrack();
    G4int tid = track->GetTrackID();
    G4int pid = track->GetParentID();


    /*PhotonCount ::iterator tit = Photons.find(tid);
    PhotonCount ::iterator pit = Parents.find(pid);


    if(tit==Photons.end()) {
        Photons[tid]=1;
        TotalPhotons++;

    } else {
        Photons[tid]+=1;
    }*/

    if(track->GetTrackStatus()==fStopAndKill){

        ProducedPhotonCount::iterator PhotonCount=DTotalPhotonCounts.find(TempName);

        PhotonCount!=DTotalPhotonCounts.end() ? DTotalPhotonCounts[PhotonCount->first]+=1 : DTotalPhotonCounts[TempName]=1;

        TotalPhotons++;
    }





    // Retrieve the pointer to the optical boundary process.
    // We do this only once per run defining our local pointer as static.
    static G4OpBoundaryProcess* boundary = 0;
    if (!boundary) { // the pointer is not defined yet
        // Get the list of processes defined for the optical photon
        // and loop through it to find the optical boundary process.
        G4ProcessVector* pv = pdef->GetProcessManager()->GetProcessList();
        for (size_t i=0; i<pv->size(); i++) {
            if ((*pv)[i]->GetProcessName() == "OpBoundary") {
                boundary = (G4OpBoundaryProcess*) (*pv)[i];
                break;
            }
        }
    }

    if (step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {
        if (boundary->GetStatus() == Detection ){
            G4String detector_name = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName();
            //G4cout << "##### Sensitive Volume: " << detector_name << G4endl;
            detectorCounts::iterator it = my_counts_.find(detector_name);
            DepositedEnergy ::iterator deit = PhotonEnergy_.find(detector_name);
            DetailedPhotonCounts ::iterator DPhoton=DPhotonCounts.find(detector_name);
            if(DisplayPhotonEnergy_){
                if(deit!=PhotonEnergy_.end())   PhotonEnergy_[it->first].push_back(step->GetTotalEnergyDeposit()/CLHEP::eV);
                else  PhotonEnergy_[detector_name].push_back(step->GetTotalEnergyDeposit()/CLHEP::eV);
            }


            //PhotonCounting
            if (it != my_counts_.end())  my_counts_[it->first] += 1;

            else my_counts_[detector_name] = 1;

            if(DPhoton!=DPhotonCounts.end()) {
                detectorCounts ::iterator ddphoton=DPhoton->second.find(TempName);
                ddphoton!=DPhoton->second.end() ?  DPhoton->second[ddphoton->first]+=1 :  DPhoton->second[TempName]=1;

            }else{
                detectorCounts tempP;
                tempP[TempName]=1;
                DPhotonCounts[detector_name]=tempP;
            }


        } else if(boundary->GetStatus()==Absorption ) {
            G4String detector_name = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName();

            if(!(detector_name=="Anode_PHOTOCATHODE" || detector_name=="Cathode_PHOTOCATHODE")){
                return;
            }
            detectorCounts::iterator Absit=ObservedPhotons.find(detector_name);
            if (Absit != ObservedPhotons.end())  ObservedPhotons[Absit->first] += 1;

            else ObservedPhotons[detector_name] = 1;

            DetailedPhotonCounts ::iterator DAbsPhoton=DAbsPhotonCounts.find(detector_name);

            if(DAbsPhoton!=DAbsPhotonCounts.end()) {
                detectorCounts ::iterator ddAbsphoton=DAbsPhoton->second.find(TempName);
                ddAbsphoton!=DAbsPhoton->second.end() ?  DAbsPhoton->second[ddAbsphoton->first]+=1 :  DAbsPhoton->second[TempName]=1;

            }else{
                detectorCounts tempAbsP;
                tempAbsP[TempName]=1;
                DAbsPhotonCounts[detector_name]=tempAbsP;
            }



        }

    }

}