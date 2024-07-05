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

#include "functions_library.cuh"

using namespace std;

class parameter_load
{
private:
    string file_location;

public:
    parameter_load(string file_Location);
    parameter_load();

    void get_parameters(int &CUDA_device_ID, string &parent_SEQ_folder,
                        float &mean_rep_time, float &standard_deviation_rep_time,
                        float &mean_days_host, float &standard_deviation_host_time);

    vector<string> get_parameters(string file_Location, vector<string> &parameters_List);

    int get_INT(string value);
    string get_STRING(string value);
    float get_FLOAT(string value);

    vector<pair<string, string>> get_block_from_File(string &file_parameter_Location, string block_Header);
    vector<pair<string, string>> get_block_from_block(vector<pair<string, string>> &block, string block_Header);

    vector<pair<string, string>> check_block_from_block(vector<pair<string, string>> &block, string block_Header);

    int get_INT(vector<pair<string, string>> block, string value);
    string get_STRING(vector<pair<string, string>> block, string value);
    float get_FLOAT(vector<pair<string, string>> block, string value);

    vector<string> clean_Line(string line, functions_library &function);

    float **get_Profile_Array(string profile_File_location, int &num_Sites, functions_library &function);
    vector<pair<int, vector<float>>> get_recombination_Hotspot_Parameters(string parameter, string hotspot, string path_to_file, functions_library &function, int &total);

};