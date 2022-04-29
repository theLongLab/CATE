#pragma once
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <iostream>
#include <cmath>
#include <cstring>
#include <bits/stdc++.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

using namespace std;

class ehh
{

private:
    string range_Mode;
    string file_Mode_path;

    int fixed_Mode_value;
    char fixed_mode_add_OR_minus;

    string input_Folder;
    string ouput_Path;
    string intermediate_Path;

    int ploidy;

    int tot_Blocks;
    int tot_ThreadsperBlock;

public:
    ehh(string range_Mode, string file_Mode_path, string fixed_Mode_value, string input_Folder, string output_Path, int cuda_ID, string intermediate_Path, int ploidy);

    void ingress();
    void process_EHH(vector<string> &total_Segregrating_sites, vector<int> &core_OR_ext, int core_Count, float N, vector<int> &core_Haplotype_Collection, vector<int> &extended_Haplotype_Sums);
};