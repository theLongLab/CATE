#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <filesystem>
#include <list>
#include <set>
#include <sstream>
#include <vector>

using namespace std;

class vcf_splitter
{

private:
    char folder[256];
    string vcf_Folder;
    string output_Path;
    string population_File;

    int sample_Column;int pop_Column;

    set<string> super_pops;
    vector<pair<string, string>> pop_Index;
    
    int REF;
    int ALT;
    int snp_Count;

public:
    vcf_splitter(char folder[], string vcf_Folder, string population_File, string output_Path, int REF, int ALT, int snp_Count,int sample_Column,int pop_Column);
    void read_File();
    void index_population();
    void find_VCF(string &vcf_Folder_Path, list<string> &vcf_Files);
    string get_Country(string &ID);
    void write_File_SNP_only(vector<string> &patient_Coutry, string &file, string &file_path);

    void split(vector<string> &line_Data, string &line);
    void write_Header(vector<string> &patient_Coutry, string &file, vector<string> &line_Data_header, string &first_SNP_size);
    string split_check(vector<string> &line_Data, string &line);
};