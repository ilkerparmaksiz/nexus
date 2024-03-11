//
// Created by ilker parmaksiz  on 2/29/24.
//

#ifndef NEXUS_OPTICKSSTEPPINGACTION_HH
#define NEXUS_OPTICKSSTEPPINGACTION_HH
#include <G4UserSteppingAction.hh>
#include <globals.hh>
#include <map>

class G4Step;


namespace nexus {
    class OpticksSteppingAction : public G4UserSteppingAction{
    public:
        /// Constructor
        OpticksSteppingAction();

        ///Destructor
        ~OpticksSteppingAction();
        virtual void UserSteppingAction(const G4Step*);
    private:
    };
}

#endif //NEXUS_OPTICKSSTEPPINGACTION_HH
