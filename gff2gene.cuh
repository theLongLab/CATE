#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <filesystem>
#include <list>
#include <set>
#include <sstream>

using namespace std;

class gff2gene
{
private:
    string input_File;
    string output_Path;

public:
    gff2gene(string input_File, string output_Path);
    void ingress();
};