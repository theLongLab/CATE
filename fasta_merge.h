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
#include <algorithm>

using namespace std;

class fasta_merge
{

private:
    string fasta_Folder;
    string output_FASTA;

public:
    fasta_merge(string fasta_Folder, string output_FASTA);
    
    void ingress();
    vector<string> index_FASTA_folder();
};