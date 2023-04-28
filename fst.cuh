#pragma once
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <cstring>
#include <set>
#include <algorithm>
#include <iterator>
#include <thread>
#include <future>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "fst_test_pop.cuh"

using namespace std;
class fst
{

private:
    string calc_Mode = "FILE";
    int window_Size = 0;
    int step_Size = 0;

    vector<string> super_Pop_Unique_vec;

    string gene_List;
    string input_Folder;
    string output_Path;
    string intermediate_Path;
    string pop_Index_path;
    string pop_List;

    int ploidy;

    int tot_Blocks;
    int tot_ThreadsperBlock;

    int max_Location_Size;
    int *locations_Size, *cuda_locations_Size;
    int *pop_seqeunce_Size_Array, *cuda_pop_seqeunce_Size_Array;
    int **sample_Location_array, **cuda_sample_Location_array;

public:
    // GENE FILE NORMAL MODE
    fst(string gene_List, string input_Folder, string output_Path, int cuda_ID, string intermediate_Path, int ploidy, string pop_Index_path, string pop_List);
    // WINDOW MODE
    fst(string calc_Mode, int window_Size, int step_Size, string input_Folder, string output_Path, int cuda_ID, int ploidy, string pop_Index_path, string pop_List);

    void set_Values(string gene_List, string input_Folder, string output_Path, int cuda_ID, string intermediate_Path, int ploidy);

    void ingress();
    void population_Processing(vector<string> &test_Pops, vector<vector<string>> &super_Pop_per_ID_full, vector<vector<string>> &sample_IDs_per_ID_full);
    void sample_location_Index(int num_Pop_Ids, fst_test_pop pop_IDs[]);

    void process_FST(fst_test_pop pop_IDs[], int num_Pop_Ids, int &Segs_count_All, float &Fst_All, float &Avg_Fst, float &numerator_Avg, float &denominator_Avg, float &ratio_of_Avg);
};