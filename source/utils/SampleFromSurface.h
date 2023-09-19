//
// Created by Ilker Parmaksiz on 8/15/23.
//

#ifndef SAMPLEFROMSURFACE_H
#define SAMPLEFROMSURFACE_H
#include "G4UnionSolid.hh"
#include "G4TessellatedSolid.hh"
#include "G4UnitsTable.hh"
#include "IOUtils.h"

namespace nexus{
    using namespace std;
    class SampleFromSurface {
    public:
        SampleFromSurface(G4String name);
        ~SampleFromSurface();

        // Random Sampling from needle
        void SampleFromFacet(G4String n1,G4TessellatedSolid * solid1);
        void SampleFromFacet(G4TessellatedSolid * solid1);
        void SampleFromFacet(G4String n1,G4TessellatedSolid * solid1,G4String n2,G4TessellatedSolid * solid2);
        void SampleFromFacet(G4String n1,G4TessellatedSolid * solid1,G4String n2,G4TessellatedSolid * solid2,G4String n3,G4TessellatedSolid * solid3);
        void FaceTransform(const G4VPhysicalVolume* tr,const G4VPhysicalVolume * Mother);
        void SaveAllPointsToOneFile();

        std::map<G4String,std::vector<G4ThreeVector>> * getRawPoints();
        std::map<G4String,std::vector<G4ThreeVector>> * getTranslatedPoints();
        // This class for determining what percent of the surface we would like to sample
        // For the needle 15% of the geometry gives me the best needle
        void SetPercent(G4double);
        G4double GetPercent();
        G4String CRABPATH;
    private:
        std::map<G4String,std::vector<G4ThreeVector>> *SamplePoints_;
        std::map<G4String,std::vector<G4ThreeVector>> *TranslatedSamplePoints_;
        G4double frac_;
        G4String name_;
        G4String AllFilePath;
        //G4String SingleFilePath;
        bool OverRide;

    };
    // Get Possible Points
    inline std::map<G4String,std::vector<G4ThreeVector>>* SampleFromSurface::getRawPoints(){
        return SamplePoints_;

    };
    // Transform the points
    inline std::map<G4String,std::vector<G4ThreeVector>>* SampleFromSurface::getTranslatedPoints(){
        return TranslatedSamplePoints_;

    };
    // Get Percent Value
    inline G4double SampleFromSurface::GetPercent() { return frac_ ;}
}
#endif
