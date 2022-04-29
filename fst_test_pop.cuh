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

using namespace std;

class fst_test_pop
{
private:
    string Pop_ID;
    vector<string> super_Pops;
    vector<string> sample_IDs;
    vector<string> countries_Folder;

    vector<vector<pair<string, string>>> folder_Index_Super_Pops;

    vector<int> sample_Location;
    int samples;
    int ploidy;
    int N;
    float N_float;
    long int combinations;

    // clear these before next gene
    vector<vector<string>> file_List_All;
    vector<vector<pair<int, string>>> collect_Segregrating_sites_All;
    vector<pair<int, string>> final_Seg_Collection;

public:
    fst_test_pop();
    fst_test_pop(string Pop_ID, vector<string> super_Pops, vector<string> sample_IDs, vector<string> countries_Folder, int ploidy);

    void index_Population();
    void index_Samples();

    void folder_Search(int start_Co, int end_Co);
    void seg_Retrival(int start_Co, int end_Co);
    void seg_Retrival_with_Found(int start_Co, int end_Co, vector<string> super_Pops_FOUND, vector<vector<pair<int, string>>> collect_Segregrating_sites_FOUND);

    void combine_Segs();

    vector<pair<int, string>> return_Seg_site(string super_Pop);

    vector<pair<int, string>> get_final_Seg_Collection();
    void clear_final_Seg_Collection();
    vector<int> get_sample_Location();
    void clear_sample_Location();
    int get_Sample_Size();
    int get_Sequence_Size();
};