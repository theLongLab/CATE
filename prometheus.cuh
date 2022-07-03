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
#include <mutex>
#include <shared_mutex>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
// #include "functions.cuh"

using namespace std;

class prometheus
{
private:
    shared_mutex g_mutex;
    mutex gi_mutex;

    int tot_Blocks;
    int tot_ThreadsperBlock;
    int CPU_cores;
    int SNPs_per_Run;
    // int number_of_genes;

    long int combinations;
    int N;
    float N_float;
    // TAJIMA
    float an, e1, e2;
    // FU LI
    float vd, ud, vd_star, ud_star, uf, vf, uf_star, vf_star;
    // FAY WU
    float bn, bn_plus1;

    vector<pair<string, string>> folder_Index;
    string Multi_read;

    // clear these after every set
    int gene_Size = 0;

    vector<int> all_start_Co;
    vector<int> all_end_Co;

    vector<int> catch_Point;
    // vector<string> catch_files_index;
    vector<vector<int>> catch_forward_index;
    vector<vector<int>> catch_back_index;
    set<int> all_Files_index;

    // set<string> all_Files;

    vector<string> all_Lines;

    int tot_Segs = 0;
    vector<pair<int, int>> start_stop;
    vector<pair<int, int>> position_index_Segs;
    vector<int> seg_catch_points_ALL;

    vector<int> seg_catch_index_ALL;
    vector<vector<int>> seg_backward_index_ALL;
    vector<vector<int>> seg_forward_index_ALL;

    vector<string> write_Lines;

    // vector<string> concat_Segs;
    // vector<int *> all_site_Index;

    vector<int *> all_end_Index;
    vector<char **> all_segs_Matrix;
    vector<int> max_Track;

public:
    // TAJIMA
    prometheus(vector<pair<string, string>> folder_Index, string Multi_read, int tot_Blocks, int tot_ThreadsperBlock, long int combinations, float an, float e1, float e2, int N, int CPU_cores, int SNPs_per_Run, int number_of_genes);
    // FU_LI
    prometheus(vector<pair<string, string>> folder_Index, string Multi_read, int tot_Blocks, int tot_ThreadsperBlock, int CPU_cores, int SNPs_per_Run, int number_of_genes, int N, long int combinations, float an, float vd, float ud, float vd_star, float ud_star, float uf, float vf, float uf_star, float vf_star);
    // FAY_WU
    prometheus(vector<pair<string, string>> folder_Index, string Multi_read, int tot_Blocks, int tot_ThreadsperBlock, int CPU_cores, int SNPs_per_Run, int number_of_genes, int N, long int combinations, float an, float bn, float bn_plus1);
    // NEUTRALITY ALL
    prometheus(vector<pair<string, string>> folder_Index, string Multi_read, int tot_Blocks, int tot_ThreadsperBlock, int CPU_cores, int SNPs_per_Run, int number_of_genes, int N, long int combinations, float an, float e1, float e2, float vd, float ud, float vd_star, float ud_star, float uf, float vf, float uf_star, float vf_star, float bn, float bn_plus1);

    // vector<string> collection_Engine(vector<string> &gene_Collect);
    vector<string> collection_Engine(vector<string> &gene_Collect, string test_Type);

    void get_Gene_info_and_catch(string gene_Combo, int gene_ID);
    void forward_Search(int pos, int start_Co, int end_Co, int gene_ID);
    void backward_Search(int pos, int start_Co, int end_Co, int gene_ID);
    void compile(int gene_ID);

    // vector<string> compound_interpolationSearch(int &start_Co, int &end_Co);

    // Multi_read
    void file_Reader_multi(string files);
    // Single_read
    void file_Reader_single();
    // void set_values(int gene_ID, string gene_Name, string chromosome, string start_Co, string end_Co);

    void process_Tajima();
    void process_Fu_Li();
    void process_Fay_Wu();
    void process_Neutrality();

    // void seg_Concat(int round_ID, int start_Seg, int stop_Seg);
    void seg_Concat_New(int round_ID, int start_Seg, int stop_Seg);
    void seg_Search_catch_point(int gene_ID);
    void seg_Search_forward(int gene_ID);
    void seg_Search_backward(int gene_ID);
    void seg_Indexer(int start_Seg, int stop_Seg, char **full_Char, int *VALID_or_NOT, int *pos_start_Index, int *pos_end_Index, int start);

    void calc_Tajima_Segs(int gene_ID, int *MA_Count);
    void calc_Fu_Li_Segs(int gene_ID, int *MA_Count, int *ne, int *ns);
    void calc_Fay_Wu_Segs(int gene_ID, int *MA_Count, int *Theta_partials);
    void calc_Neutrality_Segs(int gene_ID, int *MA_Count, int *ne, int *ns, int *Theta_partials);

    void erase();
    void set_Values(vector<pair<string, string>> folder_Index, int tot_Blocks, int tot_ThreadsperBlock, string Multi_read, int CPU_cores, int SNPs_per_Run, int number_of_genes);
};