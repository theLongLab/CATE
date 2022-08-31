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

class mk_test
{
private:
    string gene_List;
    string input_Folder;
    string ouput_Path;
    string intermediate_Path;
    string reference_Path;
    string alignment_Path;

    string primary_Intermediate_Path;

    int ploidy;

    string mode;
    string ORF_mode = "NO";

    string genetic_Code;
    string start_Codons;
    string stop_Codons;
    vector<string> start_Codons_list;
    vector<string> stop_Codons_list;

    int size_of_genetic_Code;
    char *index_Gen_code;

    int tot_Blocks;
    int tot_ThreadsperBlock;

    char *cuda_stop_Codons;
    int stop_Codon_size;
    char *cuda_reference;
    int reference_size;

public:
    mk_test(string reference_Path, string alignment_Path, string gene_List, string input_Folder, string ouput_Path, int cuda_ID, string intermediate_Path, int ploidy, string genetic_Code, string start_Codons, string stop_Codons, string mode, string ORF_mode);

    void ingress();

    void prepration();
    // CHROM mode
    void reference_Prep(vector<pair<int, int>> TEMP_file_index);
    vector<pair<int, int>> alignment_Prep();

    // GENE mode
    void reference_Prep();
    vector<pair<int, int>> alignment_Prep(string Gene_alignment_Path, string &full_Reference, int start_Co, string temp_index_Folder);

    void codon_Alignment_print(vector<string> file_List, int start_Codon, int end_Codon, string &temp_index_Folder, string &intermediate_Reference, string &file_Name);

    void process_Genetic_code();

    void process_MK();
    void process_ORF(vector<pair<int, string>> &collect_Segregrating_site_POS, int &real_segregrating_Sites, int codon_Start, int codon_Stop, string codon_Index_File_name, int &tot_Dn, int &tot_Ds, int &tot_Pn, int &tot_Ps, float &NI);
    void ORF_search(vector<int> start_ORFs, int &found, int &ORF_start, int &ORF_stop, int gene_End);

    void print_Code(vector<string> &Code_split);

    vector<pair<int, int>> index_alignment_Folder();

    set<string> get_Log();
    void log_Write(string line);

    vector<string> compound_Interpolation_folder(vector<pair<int, int>> folder_Index, int start_Co, int end_Co);
    void backward_Search(promise<vector<int>> &backward_Found, int pos, vector<pair<int, int>> folder_Index, int start_Co, int end_Co);
    void forward_Search(promise<vector<int>> &forward_Found, int pos, vector<pair<int, int>> folder_Index, int start_Co, int end_Co);
};