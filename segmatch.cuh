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

class segmatch
{
private:
    shared_mutex g_mutex;

    vector<pair<int, int>> positions_start_end;
    vector<string> sequences_to_Check_ID;
    vector<string> sequences_to_Check;

    int total_Bases = 0;

    string node_ID = "";
    int node_Index = -1;
    float cutoff = 0.000;
    int cutoff_Count = 0;
    string tissue_Name = "";
    int tissue_Index = -1;

    string output_Folder_location;
    string intermediate_Folder_location;

    string node_intermediary_location;
    string node_results_Location;

    int CPU_cores = -1;
    string Multi_Read = "YES";

public:
    segmatch(string parameter_Master_Location);

    void ingress();

    void get_Match(int start, int stop, string &sequence_Query, string header,string generation);
};