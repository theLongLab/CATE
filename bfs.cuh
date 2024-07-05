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

class bfs
{
private:
    string output_Folder_location;
    string intermediate_Folder_location;
    string pedigree_Folder_location;

    string print_Sequences = "NO";

    // string current_node_ID = "";
    vector<string> tissue_Names;

    int generation = -1;
    vector<string> search_sequence_IDs;

    vector<pair<int, string>> node_Indexes;

    string node_main_Index = "";
     int tissue_main_Index = -1;
    int re_tar_sequence_Folder;
    int re_tar_Tissue_folder;
    int re_tar_Generation;

public:
    bfs(string parameter_Master_Location);

    void ingress();

    int check_Tar_Folder(string location);
};