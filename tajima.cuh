#pragma once
#include <iostream>
#include <future>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <thread>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

using namespace std;

class tajima
{
private:
    string calc_Mode = "FILE";
    int window_Size = 0;
    int step_Size = 0;

    string gene_List;
    string input_Folder;
    string ouput_Path;
    string intermediate_Path;

    int ploidy;

    int tot_Blocks;
    int tot_ThreadsperBlock;

    string prometheus_Activate = "NO";
    string Multi_read = "NO";
    int number_of_genes = 0;
    int CPU_cores = 0;
    int SNPs_per_Run = 0;

public:
    tajima();
    // Normal GENE mode
    tajima(string gene_List, string input_Folder, string ouput_Path, int cuda_ID, string intermediate_Path, int ploidy);
    // PROMETHEUS GENE mode
    tajima(string gene_List, string input_Folder, string ouput_Path, int cuda_ID, string intermediate_Path, int ploidy, string prometheus_Activate, string Multi_read, int number_of_genes, int CPU_cores, int SNPs_per_Run);
    // WINDOW mode PROMETHEUS
    tajima(string calc_Mode, int window_Size, int step_Size, string input_Folder, string ouput_Path, int cuda_ID, int ploidy, string prometheus_Activate, string Multi_read, int number_of_genes, int CPU_cores, int SNPs_per_Run);

    void set_Values(string gene_List, string input_Folder, string ouput_Path, int cuda_ID, string intermediate_Path, int ploidy);
    void ingress();

    vector<string> get_Countries();
    vector<pair<string, string>> index_Folder(string &country);

    vector<string> compound_interpolationSearch(vector<pair<string, string>> &folder_Index, int &start_Co, int &end_Co);
    vector<int> forward_Search(int pos, vector<pair<string, string>> folder_Index, int start_Co, int end_Co);
    vector<int> backward_Search(int pos, vector<pair<string, string>> folder_Index, int start_Co, int end_Co);

    int getN_Split(string file);
    void calc_Pre(int &N_tot, float &a1, float &e1, float &e2); // calc a1 to e2
    int calc_Pairwise(string &line, int pair_Count);

    void split(vector<string> &line_Data, string line, string delim);
    void split_Convert(int *line_temp, string line, string delim);
    void split_getPos(vector<string> &line_Data, string line, string delim);
    long int fact_half(int count);
    long int combos_N(int count);
    void createFile(string path, string text);
    void createFile(string path);
};