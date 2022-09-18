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
    /**
     * Prints a default parameter.json file for the user.
     * User can then customize and use this parameter file.
     **/

private:
    /**
     * @param file_Name is sued to print the parameter file to a specific location.
     * This is the only instance where the paramater file location is used to point to a location where the ile will be printed to
     * and does not already exist.
     **/
    string file_Name;

public:
    /**
     * Constructor Function assigns passed variables to the private variables and initiates the function.
     **/
    print_param(string file_Name);

    /**
     * Execution function prints the paramater file to the user defined location in *.json format.
     **/
    void ingress();
};