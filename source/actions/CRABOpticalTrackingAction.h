//
// Created by ilker on 7/19/22.
//

#ifndef NEXUS_CRABOPTICALTRACKINGACTION_H
#define NEXUS_CRABOPTICALTRACKINGACTION_H


#include <G4UserTrackingAction.hh>
#include "G4String.hh"
class G4Track;
class G4GenericMessenger;


namespace nexus {

    // Optical-checking user tracking action

    class CRABOpticalTrackingAction: public G4UserTrackingAction
    {
    public:
        /// Constructor
        CRABOpticalTrackingAction();
        /// Destructor
        virtual ~CRABOpticalTrackingAction();

        virtual void PreUserTrackingAction(const G4Track*);
        virtual void PostUserTrackingAction(const G4Track*);
    private:
        G4GenericMessenger *msg_;
        G4String Particle;
    };

}
#endif //NEXUS_CRABOPTICALTRACKINGACTION_H
