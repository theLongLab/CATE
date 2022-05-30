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

class print_param
{
    private:
    
    string file_Name;

    public:

    print_param(string file_Name);
    void ingress();

};