//
// Created by ilker on 6/27/22.
//

#include <G4SIunits.hh>
#include "FileHandling.hh"
#include "G4Exception.hh"
namespace filehandler{
    using namespace CLHEP;
    // construct
    FileHandling::FileHandling() {

    }
    //Destruct
    FileHandling::~FileHandling() {
        if (fstfile.is_open()) fstfile.close();
        if(ifsfile.is_open()) ifsfile.close();
    }

    // This function imports data with two column with any delimeter
    vector<vector<G4double>>  FileHandling::GetData(string file, char del,G4int SkipRow=1)
    {
        string str;
        ifsfile =ifstream (file);

        if (!ifsfile.is_open()) {
            G4Exception("FileHandling","[Get2ColumnData]",FatalException,"Could not open the file!");
        }
        vector<vector<G4double>> Data2d;
        vector<G4double>c1,c2;

        G4int SkipCount=0;

        while(getline(ifsfile,str)){
            string val;
            stringstream sline(str);
            if(SkipRow!=0 and SkipCount<SkipRow) {
                G4cout<<"Skipping following lines "<<G4endl;
                G4cout<<str<< " is Skipped "<<G4endl;
                SkipCount++;
            } else{
                int Counter=0;
                while (getline(sline,val,del)){
                    if(Counter==0)
                        c1.push_back(stof(val));
                    else
                        c2.push_back(stof(val));
                    Counter++;
                }

            }

        }
        Data2d.push_back(c1);
        Data2d.push_back(c2);
        return Data2d;

    }
    void FileHandling::SaveToTextFile(string file, string labels, char del , std::vector<vector<G4double>>data) {
        string val;
        fstfile = fstream (file);
        if (!fstfile.is_open()) G4Exception("FileHandling","[SaveToTextFile]",FatalException,"Couldnt open the file!");


        stringstream sline(labels);
        G4int counter=0;

        // Checking to see if something written in the file, if it is skip labels
        if(fstfile.peek()==std::fstream::traits_type::eof()) {

            // These are needed to point to the begining of the empty page
            fstfile.seekp(std::ios::beg);
            fstfile.seekg(std::ios::beg);

            //////////////////////////////////////////////////////////////////////////


            // Just a Warning
            G4Exception("FileHandling","[SaveToTextFile]",JustWarning,"Empty File!");

            // If the labels are not present,  just add them
            if(!labels.empty()){
                fstfile << labels +"\n";
            }
        }

        if (data.size()==0) G4Exception("FileHandling","[SaveToTextFile]",FatalException,"Data array is empty!");

        for (int i=0; i<data.size();i++){
            string str;

            for (int k=0;k<data.at(i).size();k++){
                if(k==0)
                    str=to_string(data.at(i).at(k));
                else if(k==data.at(i).size()-1)
                    str=str+del+to_string(data.at(i).at(k)) +"\n";
                else
                    str=str+del+ to_string(data.at(i).at(k));
            }
            fstfile << str;
        }

        fstfile.close();
    }

    void FileHandling::SaveToTextFile(string file, string labels, G4String data) {
        string val;
        if(!fstfile.is_open()) fstfile=fstream (file.c_str(),std::fstream::in | std::fstream::out | std::fstream::app);

        if (!fstfile.is_open()) G4Exception("FileHandling","[SaveToTextFile]",FatalException,"Couldnt open the file!");


        stringstream sline(labels);
        G4int counter=0;

        // Checking to see if File is Empty
        if(fstfile.peek()==std::fstream::traits_type::eof()) {

            // These are needed to point to the begining of the empty page
            fstfile.seekp(std::ios::beg);
            fstfile.seekg(std::ios::beg);

            //////////////////////////////////////////////////////////////////////////


            // Just a Warning
            G4Exception("FileHandling","[SaveToTextFile]",JustWarning,"Empty File!");

            // If the labels are not present,  just add them
            if(!labels.empty()){
                fstfile << labels +"\n";
            }
        }


        if (data.size()==0) G4Exception("FileHandling","[SaveToTextFile]",FatalException,"Data array is empty!");
        fstfile << data <<"\n";

        G4cout<<data << " is written to "<< file <<G4endl;
        fstfile.close();
    }

    void FileHandling::GetTimeProfileData(string filename, vector<vector<vector<G4double>>> &data, vector<G4double> &events) {

        std::cout << "[File Handling] Loading EL timing profiles..." << std::endl;
        

        // Open the file
        std::ifstream FileIn_(filename);

        // Check if file has opened properly
        if (!FileIn_.is_open()){
        G4Exception("[FileHandling]", "GetTimeProfileData()",
                    FatalException, " could not read in the CSV file ");
        }

        // Read the Data from the file as strings
        std::string s_event, s_x, s_y, s_z, s_t;

        G4double event, x, y, z, t;

        std::vector<G4double> temp_data;
        std::vector<std::vector<G4double>> temp_profile;

        G4double temp_event = -1;

        // Loop over the lines in the file and add the values to a vector
        while (FileIn_.peek()!=EOF) {

            std::getline(FileIn_, s_event, ',');
            std::getline(FileIn_, s_x, ',');
            std::getline(FileIn_, s_y, ',');
            std::getline(FileIn_, s_z, ',');
            std::getline(FileIn_, s_t, '\n');

            // Set it to the first event
            if (temp_event == -1)
                temp_event = stod(s_event);

            if (temp_event != stod(s_event)){
                events.push_back(stod(s_event));
                temp_event = stod(s_event);
                data.push_back(temp_profile);
                temp_profile.clear();
            }

            temp_data.push_back(stod(s_x));
            temp_data.push_back(stod(s_y));
            temp_data.push_back(stod(s_z));
            temp_data.push_back(stod(s_t));

            temp_profile.push_back(temp_data);
            temp_data.clear();

        } // END While

        FileIn_.close();

        std::cout << "[File Handling] Finished loading EL timing profiles..." << std::endl;


    }

    void FileHandling::GetEvent(string filename, vector<vector<G4double>> &data) {

        std::cout << "[File Handling] Loading event data..." << std::endl;
        

        // Open the file
        std::ifstream FileIn_(filename);

        // Check if file has opened properly
        if (!FileIn_.is_open()){
        G4Exception("[FileHandling]", "GetEvent()",
                    FatalException, " could not read in the CSV file ");
        }

        // Read the Data from the file as strings
        std::string s_event, s_x, s_y, s_z, s_e;

        std::vector<G4double> temp_data;

        // Loop over the lines in the file and add the values to a vector
        while (FileIn_.peek()!=EOF) {

            std::getline(FileIn_, s_event, ',');
            std::getline(FileIn_, s_x, ',');
            std::getline(FileIn_, s_y, ',');
            std::getline(FileIn_, s_z, ',');
            std::getline(FileIn_, s_e, '\n');

            temp_data.push_back(stod(s_x));
            temp_data.push_back(stod(s_y));
            temp_data.push_back(stod(s_z));
            temp_data.push_back(stod(s_e));

            data.push_back(temp_data);

            temp_data.clear();

        } // END While

        FileIn_.close();

        std::cout << "[File Handling] Finished loading event data..." << std::endl;


    }


} // Namespace nexus is closed