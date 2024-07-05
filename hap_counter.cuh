#pragma once
#include <iostream>

#include "functions_library.cuh"
#include "parameter_load.h"

#include <cstdlib>

#include <curand.h>
#include <curand_kernel.h>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <thrust/system_error.h>
#include <thrust/system/cuda/error.h>

#include <sstream>

#include <algorithm>
#include <random>
#include <chrono>
#include <iomanip>
#include <string>
#include <map>

#include <string>
#include <vector>
#include <queue>

#include <thread>
#include <mutex>
#include <shared_mutex>

using namespace std;

class hap_counter
{
private:
    shared_mutex g_mutex;
    mutex gi_mutex;

    string output_Folder_location;
    string intermediate_Folder_location;

    string intermediary_Sequence_location;
    string intermediary_Index_location;

    string node_Master_location;

    int CUDA_device_number;

    int *CUDA_device_IDs;
    int num_Cuda_devices;
    int gpu_Limit;
    int *tot_Blocks;
    int *tot_ThreadsperBlock;

    int CPU_cores;
    int max_Cells_at_a_time = 0;

    string multi_Read;

    string nodes_to_Analyse = "";

    vector<pair<int, string>> all_Hap_Count;
    vector<pair<int, string>> all_Hap_Alive_Count;
    vector<pair<int, string>> all_Hap_Parent_Count;

    int all_Hap_Total = 0;
    int all_Hap_Alive_Total = 0;
    int all_Hap_Parent_Total = 0;

public:
    hap_counter(string parameter_Master_Location);

    void ingress();

    void all_Haplotype_Counter(vector<pair<string, string>> &line_Data);
    void all_Haplotype_Alive_Counter(vector<pair<string, string>> &line_Data, functions_library &functions);
    void all_Haplotype_Parent_Counter(vector<pair<string, string>> &line_Data, functions_library &functions,
                                      vector<vector<vector<int>>> &tissue_generation_Sequence,
                                      int tissue, int generation, int &track_Seq);

    void write_Files(string tissue_Name, int generation, vector<pair<int, string>> &Hap_count, string location_Frequencies, string location_Summaries,int total);
};