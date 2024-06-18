// ----------------------------------------------------------------------------
// nexus | PersistencyManager.h
//
// This class writes all the relevant information of the simulation
// to an ouput file.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#ifndef PERSISTENCY_MANAGER_H
#define PERSISTENCY_MANAGER_H

#include "PersistencyManagerBase.h"

#include <G4VPersistencyManager.hh>
#include <map>
#include <vector>
#include "fstream"
#include "G4ThreeVector.hh"

class G4GenericMessenger;
class G4TrajectoryContainer;
class G4HCofThisEvent;
class G4VHitsCollection;

namespace nexus {
  class HDF5Writer;
  class IonizationHit;
}

namespace nexus {

    struct OpticalHit{
        G4String name;
        G4double time;
        G4ThreeVector position;
    };
  class PersistencyManager: public PersistencyManagerBase
  {
  public:
    PersistencyManager();
    ~PersistencyManager();
    PersistencyManager(const PersistencyManager&);

    /// Create the singleton instance of the persistency manager
   // static void Initialize(G4String init_macro, std::vector<G4String>& macros,
    //                       std::vector<G4String>& delayed_macros);

    /// Set whether to store or not the current event
    void StoreCurrentEvent(G4bool);
    void InteractingEvent(G4bool);
    void StoreSteps(G4bool);
    void SaveNumbOfInteractingEvents(G4bool);
    void SaveTimeInfo();
      ///
    virtual G4bool Store(const G4Event*);
    virtual G4bool Store(const G4Run*);
    virtual G4bool Store(const G4VPhysicalVolume*);

    virtual G4bool Retrieve(G4Event*&);
    virtual G4bool Retrieve(G4Run*&);
    virtual G4bool Retrieve(G4VPhysicalVolume*&);

  public:
    void StoreOpticksHits();
    void OpenFile();
    void CloseFile();


    void SetCompletionTime(G4double);
    void AddPhotons(int64_t);
    void SetFstream(std::fstream *fstream);
    std::fstream *GetFstream();

      void AddOpticalHit(G4String name,G4ThreeVector position,G4double time);


  private:
    void StoreTrajectories(G4TrajectoryContainer*);
    void StoreHits(G4HCofThisEvent*);
    void StoreIonizationHits(G4VHitsCollection*);
    void StoreSensorHits(G4VHitsCollection*);
    void StoreSteps();
    void StoreOpticalHits();




    void SaveConfigurationInfo(G4String history);


  private:
    G4GenericMessenger* msg_; ///< User configuration messenger

    std::vector<G4String> secondary_macros_;

    G4String output_file_; ///< Path of output file
    G4bool ready_;     ///< Is the PersistencyManager ready to go?
    G4bool store_evt_; ///< Should we store the current event?
    G4bool store_steps_; ///< Should we store the steps for the current event?
    G4bool interacting_evt_; ///< Has the current event interacted in ACTIVE?
    G4bool save_ie_numb_; ///< Should we save the number of interacting events in the configuration table?

    G4String event_type_; ///< event type: bb0nu, bb2nu, background or not set

    int64_t saved_evts_; ///< number of events to be saved
    int64_t interacting_evts_; ///< number of events interacting in ACTIVE
    G4double pmt_bin_size_, sipm_bin_size_; ///< bin width of sensors

    int64_t nevt_; ///< Event ID
    int64_t start_id_; ///< ID for the first event in file
    G4bool first_evt_; ///< true only for the first event of the run

    HDF5Writer* h5writer_;  ///< Event writer to hdf5 file

    std::vector<G4int>* ihits_;
    std::map<G4int, std::vector<G4int>* > hit_map_;
    std::vector<G4int> sns_posvec_;

    std::map<G4String, G4double> sensdet_bin_;

    // Studying timing
    G4double EventCompletionTime;
    int64_t photonCount;
    std::fstream *fstream_;


      std::vector<OpticalHit*> AllOpticalHits;
  };


  // INLINE DEFINITIONS //////////////////////////////////////////////

  inline void PersistencyManager::StoreCurrentEvent(G4bool sce)
  { store_evt_ = sce; }
  inline void PersistencyManager::StoreSteps(G4bool ss)
  { store_steps_ = ss; }
  inline void PersistencyManager::InteractingEvent(G4bool ie)
  { interacting_evt_ = ie; }
  inline void PersistencyManager::SaveNumbOfInteractingEvents(G4bool sie)
  {save_ie_numb_ = sie;}
  inline G4bool PersistencyManager::Store(const G4VPhysicalVolume*)
  { return false; }
  inline G4bool PersistencyManager::Retrieve(G4Event*&)
  { return false; }
  inline G4bool PersistencyManager::Retrieve(G4Run*&)
  { return false; }
  inline G4bool PersistencyManager::Retrieve(G4VPhysicalVolume*&)
  { return false; }
  inline void PersistencyManager::SetFstream(std::fstream *fstream) {
      fstream_=fstream;
  }
  inline std::fstream * PersistencyManager::GetFstream(){
    return fstream_;
  }


    inline void PersistencyManager::AddPhotons(int64_t  ph) {
        photonCount+=ph;
    }
    inline void PersistencyManager::SetCompletionTime(G4double t) {
        EventCompletionTime=t;
    }
    inline void PersistencyManager::AddOpticalHit(G4String name,G4ThreeVector position,G4double time){
        OpticalHit * ahit= new OpticalHit();
        ahit->name=name;
        ahit->time=time;
        ahit->position=position;
        AllOpticalHits.push_back(ahit);
  }

} // namespace nexus

#endif
