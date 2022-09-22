#pragma once
#include <iostream>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <thread>
#include <future>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

using namespace std;

class functions
{
private:
public:
    functions();

    int calc_Pairwise(string &line, int N, int tot_Blocks, int tot_ThreadsperBlock);
    void split_Convert(int *line_temp, string line, string delim);

    /**
     * Administrative function responsible for calculating the values needed for Tajima's D.
     * Function directly calls upon the GPU function.
     **/
    void process_Seg_sites_tajima(vector<string> &total_Segregrating_sites, int N, int &real_segregrating_Sites, float &tot_pairwise_Differences, int tot_Blocks, int tot_ThreadsperBlock);
    /**
     * Administrative function responsible for calculating the values needed for Fu and Li values.
     * Function directly calls upon the GPU function.
     **/
    void process_Seg_sites_fu_li(vector<string> &total_Segregrating_sites, float N, int &real_segregrating_Sites, float &tot_pairwise_Differences, int &singletons_ne, int &singletons_ns, int tot_Blocks, int tot_ThreadsperBlock);

    vector<string> get_Countries(string &input_Folder);
    vector<pair<string, string>> index_Folder(string &country);

    vector<string> compound_interpolationSearch(vector<pair<string, string>> &folder_Index, int &start_Co, int &end_Co);
    void forward_Search(int pos, vector<pair<string, string>> folder_Index, int start_Co, int end_Co, promise<vector<int>> &forward_Found);
    void backward_Search(int pos, vector<pair<string, string>> folder_Index, int start_Co, int end_Co, promise<vector<int>> &backward_Found);
    // vector<string> forward_Search_Only(int pos, vector<pair<string, string>> &folder_Index, int start_Co, int end_Co);

    void split(vector<string> &line_Data, string line, char delim);
    void split_space(vector<string> &line_Data, string line, string delim);
    void split_getPos_ONLY(vector<string> &line_Data, string line, char delim);
    // void split_getPos(vector<string> &line_Data, string line, string delim);
    int getN_Split(string file);
    int get_Valid(string line);
    void split_to_MA(vector<string> &line_Data, string line, char delim);

    long int fact_half(int count);
    long int combos_N(int count);
    float add(int N, float *array);
    string roundoff(float value, unsigned char prec);

    void createFile(string path, string text);
    void createFile(string path);
};