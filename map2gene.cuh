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

class map2gene
{
private:
    string input_map_File_path;
    string ouput_Path;

    string pre_fix;

    vector<string> chromosomes;
    vector<int> count;

public:
    map2gene(string input_map_File_path, string ouput_Path, string pref_fix);

    void ingress();

    int chr_Index(string chr);
};