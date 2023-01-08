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

class vcf_splitter_2
{
private:
    shared_mutex g_mutex;
    mutex gi_mutex;

    string input_vcf_Folder;
    string output_Folder;

    string population_File_path;
    int column_Sample_ID;
    int column_Population_ID;

    // int max_SNPs_per_file;

    // 0 = equal;
    // 1 = greater than
    // 2 = less than
    // 10 = greater than or equal
    // 20 = less than or equal
    int logic_MAF = 0;
    double MAF = 0.0;

    int summary_Individuals = 0;

    int tot_Blocks;
    int tot_ThreadsperBlock;

    int ploidy;
    int hap_Size;

    int allele_Count_REF = 1;
    int allele_Count_ALT = 1;

    int SNP_count_per_File = 1000;

    int cores = 1;
    int SNPs_per_time_CPU = 1000;
    int SNPs_per_time_GPU = 100;

    // CLEAR ALL below

    /**
     * @param concat_Segs stores the concatenated Seg sites for GPU processing
     * @param all_site_Index stores the index data for the concatenated Seg array.
     **/
    vector<string> concat_Segs;
    vector<int *> all_site_Index;

    vector<pair<int, int>> start_stop;

    vector<string> all_Lines;

    vector<int> VALID_List;

    vector<string> unique_CHRs;

    vector<vector<int>> CHR_collected_SEGs;
    vector<pair<string, int>> CHR_unbound;

    vector<string> write_CHR;

    // string trim_Seg_full;
    vector<string> seg_Collection;

    // clear after whole file is read;
    vector<string> CHR_individuals;
    vector<vector<int>> CHR_individuals_COUNT;
    vector<int> CHR_full_Segs;

    vector<string> header_Data;

    vector<string> population_Unique_IDs;

    int *sample_ID_population_ID, *cuda_sample_ID_population_ID;
    vector<vector<int>> population_sample_IDs;

    vector<string> pop_Header;

    int *MAF_count_per_Population, *cuda_MAF_count_per_Population;

    string file_CHR_value;

    vector<vector<pair<int, string>>> position_Complete_segs;
    vector<int> POS_valid;

public:
    // Split by chromosome
    vcf_splitter_2(int cuda_ID, string input_vcf_Folder, string output_Folder, int cores, int SNPs_per_time_CPU, int SNPs_per_time_GPU, int allele_Count_REF, int allele_Count_ALT, int ploidy, int summary_Individuals);
    vcf_splitter_2(int cuda_ID, string input_vcf_Folder, string output_Folder, string population_File, int sampled_ID_col, int pop_ID_column, int cores, int SNPs_per_time_CPU, int SNPs_per_time_GPU, int ploidy, int max_SNPs_per_file, int logic_MAF, double MAF);

    void cuda_Set_device(int cuda_ID);

    void ingress_chr_Split();
    void ingress_file_hierarchy();

    void process_SNPs_CHR(int N, int augment, string &output_vcf_Folder, string &head_Line);

    void process_SNPs_Hierarchy(int N, int num_Populations, int augment, string output_vcf_Folder);

    /**
     * Concatenates the Segregating sites for processing by each GPU round.
     **/
    void seg_Concat(int round_ID, int start_Seg, int stop_Seg);

    void seg_VALID_OR_NOT_list(int start_Seg, int stop_Seg, int *VALID_or_NOT, int N_individuals, string seg_Full, int augment);
    void seg_CHR_get(int start_Seg, int stop_Seg, char *full_Char, int *CHR_start_Index, int *CHR_end_Index);

    void concat_ALL(string CHR_name, int start_Seg, int stop_Seg, vector<int> collected_Segs, char *full_Char, int *pos_start_Index, int *pos_end_Index, int *ID_start_Index, int *ID_end_Index, int *REF_start, int *REF_stop, int *ALT_start, int *ALT_stop, int *six_8_start_Index, int *six_8_stop_Index);
    void concat_ALL_hierarchy(int start_Seg, int stop_Seg, vector<int> VALID_Segs, char *full_Char, int *chr_start_Index, int *chr_end_Index, int *pos_start_Index, int *pos_end_Index, int *ID_start_Index, int *ID_end_Index, int *REF_start, int *REF_stop, int *ALT_start, int *ALT_stop, int *six_9_start_Index, int *six_9_stop_Index);

    void summary_Individuals_process(int CHR_ID, int start_Individual, int stop_Individual, int **sample_sequence_Tracker, vector<int> collected_Segs);

    void individual_Map(int start_N, int stop_N, vector<string> header_Data, vector<pair<string, string>> sample_population, vector<string> population_Unique_IDs);

    void complete_Pop_segs(int pop_ID, int start_Seg, int stop_Seg, vector<int> VALID_Segs, int **VALID_or_NOT_populations, vector<int> population_sample_IDs);

    void write_segs_to_Hierarchy(int pop_ID, string output_vcf_Folder);
};