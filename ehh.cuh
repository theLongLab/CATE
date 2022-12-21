#pragma once
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include <thread>
#include <mutex>
#include <shared_mutex>
#include <set>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <bits/stdc++.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

using namespace std;

class ehh
{

private:
    shared_mutex g_mutex;
    mutex gi_mutex;

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

    int default_SNP_count = 100;
    double EHH_cutoff = 0.05;

    int CPU_cores = 1;

    // Multithread global vectors

    vector<string> EHH_0_up;
    vector<string> EHH_1_up;

    vector<string> EHH_0_down;
    vector<string> EHH_1_down;

    // vector<string> EHH_values;

    vector<string> positions_Collect;

public:
    ehh(string range_Mode, string file_Mode_path, string fixed_Mode_value, string input_Folder, string output_Path, int cuda_ID, string intermediate_Path, int ploidy, int default_SNP_count, int EHH_CPU_cores);

    void ingress();

    void process_EHH(vector<string> &total_Segregrating_sites, vector<int> &core_OR_ext, int core_Count, float N, vector<int> &core_Haplotype_Collection, vector<int> &extended_Haplotype_Sums);

    void process_SNP_EHH(int SNP_position, vector<pair<string, string>> &folder_Index, int segment_Position, int SNP_Index_in_file, vector<string> &segment_Segregrating_sites, float N, int num_top, int num_bottom);

    void calc_EHH_0_1_up(vector<string> Haplotypes, int position_of_core_SNP, int augment_start, int augment_stop, int combo_zero_one, int zero_one);
    void calc_EHH_0_1_down(vector<string> Haplotypes, int position_of_core_SNP, int augment_start, int augment_stop, int combo_zero_one, int zero_one);

    // void calc_EHH_0_1(vector<string> Haplotypes, int position_of_core_SNP, char up_down, int augment, int one_zero, int combo_one_zero);

    void pos_Construct(char *full_Char, int *pos_start_Index, int *pos_end_Index, int gen_position_of_Core,int start_pos_Array,int stop_pos_Array);
};