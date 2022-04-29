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

class fasta_splitter
{

private:
    string fasta_File;
    string output_Folder;
    string fasta_Name;

public:
    fasta_splitter(string fasta_File, string output_Path, string fasta_Name);
    void split_all();
    void split_select();
    void ingress();
};