#pragma once
#include <iostream>
#include <future>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <cmath>
#include <thread>
// #include <mutex>
// #include <shared_mutex>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

using namespace std;

class hap_extract
{
private:
    string input_Folder;
    string output_Path;
    string intermediate_Path;

    string gene_List;
    string reference_File;

    int ploidy;
    int N;

    int tot_Blocks;
    int tot_ThreadsperBlock;

    char *cuda_reference;
    int reference_size;

public:
    hap_extract(string gene_List, string input_Folder, string output_Path, int cuda_ID, string intermediate_Path, int ploidy, string reference_File);
    void set_Values(string gene_List, string input_Folder, string output_Path, int cuda_ID, string intermediate_Path, int ploidy);

    void ingress();

    void hap_extraction(vector<string> &write_Lines, vector<string> &write_Sequences, vector<string> &total_Segregrating_sites, vector<pair<int, int>> &pos_INDEX, string gene_Name, string chr, int start_Pos, int end_Pos);
};