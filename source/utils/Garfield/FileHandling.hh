//
// Created by ilker on 6/27/22.
// This class is created for Importing data tables from a txt,and csv files for material properties
//

#ifndef NEXUS_IMPORTFROMAFILE_H
#define NEXUS_IMPORTFROMAFILE_H
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "G4UnitsTable.hh"
using namespace std;

namespace filehandler{
    class FileHandling {
        public:
            //Construct
            FileHandling();
            //Destruct
            ~FileHandling();
            // Reading 2 Column data like ev and Absorbtion Length
            vector<vector<G4double>> GetData(string file, char del,G4int SkipRow);

            // Load in data for photon yields
            void GetTimeProfileData(string filename, vector<vector<vector<G4double>>> &data, vector<G4double> &events);

            // Get event data from csv
            void GetEvent(string filename, vector<vector<G4double>> &data);


            // This is for wrting detector counts to a text file
            void SaveToTextFile(string file,string labels, char del, std::vector<vector<G4double>>data);
            void SaveToTextFile(string file,string labels, G4String data);

        private:
        ifstream ifsfile;
        fstream fstfile;
    };
}

#endif //NEXUS_IMPORTFROMAFILE_H
