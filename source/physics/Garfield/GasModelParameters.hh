#ifndef GasModelParameters_hh
#define GasModelParameters_hh

#include "G4SystemOfUnits.hh"
#include "G4String.hh"
#include <map>

class DegradModel;
class GasModelParametersMessenger;
class DetectorConstruction;
class G4String;
namespace nexus{
    class GasModelParameters{
        public:

        GasModelParameters();
        ~GasModelParameters();

        /*Getters and Setters*/
        inline void SetThermalEnergy(G4double d){thermalE=d;}
        inline G4double GetThermalEnergy(){return thermalE;};
        inline void SetComsol(G4bool b){useComsol_=b;};
        inline void SetEL_file(G4bool b){useEL_File_=b;};
        inline void SetCOMSOL_Path(G4String s){COMSOL_Path_=s;};
        inline bool GetbComsol(){return useComsol_;};
        inline bool GetbEL_File(){return useEL_File_;};
        inline G4String GetCOMSOL_Path(){return COMSOL_Path_;};


        private:
        GasModelParametersMessenger* fMessenger;
        G4double thermalE;

        G4String COMSOL_Path_;
        G4bool 	useEL_File_;
        G4bool 	useComsol_;

    };
}
#endif
