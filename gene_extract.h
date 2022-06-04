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

class gene_extract
{

private:
    string gene_List;
    string reference_Path;
    string output_Path;
    string intermediate_Path;

public:
    gene_extract(string gene_List, string reference_Path, string output_Path, string intermediate_Path);

    void ingress();
};