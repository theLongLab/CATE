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
    /**
     * * This is the administrative class for Prometheus.
     * Prometheus is CATE's high performance mode.
     * It is available at present for the first 3 neutrality tests only and therefore, the four neutrality test functions.
     * Namely: Tajima's D, Fay and Wu and Fu and Li statistics are equipped with Prometheus.
     **/

    /**
     * ! We use the mutex form of multithreading for Prometheus.
     * ! It provided an ease of coding since Prometheus' multithreading is confined to within the class and functions.
     * Most variables used in Prometheus are global variables for the class.
     **/
    shared_mutex g_mutex;
    mutex gi_mutex;

    /**
     * @param tot_Blocks defines number of GPU blocks that are available
     **/
    int tot_Blocks;
    /**
     * @param tot_ThreadsperBlock defines number of threads that are available per GPU block.
     **/
    int tot_ThreadsperBlock;

    /**
     * ! Parameters unique to Prometheus.
     * @param CPU_cores defines the maximum number of CPU cores that can be used at any given time.
     * @param SNPs_per_Run defines the maximum number of SNPs that can be processed by the GPU at any given time.
     * @param number_of_genes_Window defines the maximum number of query regions that can be processed at a time.
     **/
    int CPU_cores;
    int SNPs_per_Run;
    int number_of_genes_Window = 0;

    // string intermediate_File;
    /**
     * @param output_File defines the path of the file to which the outputs will be written to.
     * * It is used for Window modes resume function.
     **/
    string output_File;

    // file check

    /**
     * ! If the next batch of query regions match the previous, then the SNPs from the previous read will be used.
     * @param prev_file_List is used to track the segment files used to process a specific batch of query regions.
     * @param pre_MA, @param pre_Theta_partials, @param pre_ne, @param pre_ns: GPU processed SNP values will be carried forward.
     * @param same_Files acts as a Boolean variable and will indicate the GPU to process the new list of SNPs or whether to use the existing list.
     **/
    vector<string> prev_file_List;
    int *pre_MA, *pre_Theta_partials, *pre_ne, *pre_ns;
    string same_Files = "NO";

    // WINDOW parameters
    /**
     * @param calc_Mode acts as a boolean variable to determine between GENE(FILE) mode or WINDOW mode.
     * @param sliding_Mode acts as a boolean variable and will indicate if Sliding Window mode is activated or not.
     * @param window_Size defines the base pair (bp) size of the query window for analysis. Used with WINDOW mode.
     * @param step_Size defines the number of bps to skip to reach the next query window. Used with WINDOW mode.
     **/
    string calc_Mode = "FILE";
    string sliding_Mode = "NO";
    int window_Size;
    int step_Size;

    /**
     * Prerequisite values for each neutrality test. Will be calculated in the parent class and passed to the Prometheus class.
     **/
    long int combinations;
    int N;
    float N_float;
    // TAJIMA
    float an, e1, e2;
    // FU LI
    float vd, ud, vd_star, ud_star, uf, vf, uf_star, vf_star;
    // FAY WU
    float bn, bn_plus1;

    /**
     * @param folder_Index indexed file structure data.
     * @param Multi_read boolean variable determines whether SSD based multiple reads of files should be carried out.
     **/
    vector<pair<string, string>> folder_Index;
    string Multi_read;

    // clear these after every set
    /**
     * @param gene_Size keeps track of the number of query regions collected.
     **/
    int gene_Size = 0;

    /**
     * @param all_start_Co keeps track of all the start coordinates in the collected query regions.
     * @param all_end_Co keeps track of all the end coordinates in the collected query regions.
     **/
    vector<int> all_start_Co;
    vector<int> all_end_Co;

    /**
     * ! Variables are part of the multithreaded CIS algorithm
     * @param catch_Point records all catch point file segments for each query regions.
     * @param catch_forward_index records all file segments forward in space for each query regions.
     * @param catch_back_index records all file segments backward in space for each query regions.
     * @param all_Files_index records all file segments required for all query regions, ensuring no redundancy in file segments.
     **/
    vector<int> catch_Point;
    // vector<string> catch_files_index;
    vector<vector<int>> catch_forward_index;
    vector<vector<int>> catch_back_index;
    set<int> all_Files_index;

    // set<string> all_Files;
    /**
     * @param all_Lines records all lines read from the file segments.
     **/
    vector<string> all_Lines;

    /**
     * @param tot_Segs keeps track of the total number of segregating sites being processed.
     * @param start_stop used to keep track of the start and stop seg numbers for each thread when concatenating the seg data for GPU use.
     * @param position_index_Segs keeps track of the Position of each SNP and there location in the overall arrays.
     **/
    int tot_Segs = 0;
    vector<pair<int, int>> start_stop;
    vector<pair<int, int>> position_index_Segs;

    /**
     * @param seg_catch_points_ALL is sued to ensure that a catch point was found for a particular query region.
     *
     * @param seg_catch_index_ALL records the location of the catch point index from the CBS (Compound Binary Search).
     * @param seg_backward_index_ALL records the location of the indexes from the CBS backward in space.
     * @param seg_forward_index_ALL records the location of the indexes from the CBS forward in space.
     **/
    vector<int> seg_catch_points_ALL;

    vector<int> seg_catch_index_ALL;
    vector<vector<int>> seg_backward_index_ALL;
    vector<vector<int>> seg_forward_index_ALL;

    /**
     * @param write_Lines records the lines that will be written to the output file.
     **/
    vector<string> write_Lines;

    /**
     * @param concat_Segs stores the concatenated Seg sites for GPU processing
     * @param all_site_Index stores the index data for the concatenated Seg array.
     **/
    vector<string> concat_Segs;
    vector<int *> all_site_Index;

    // vector<int *> all_end_Index;
    // vector<char **> all_segs_Matrix;
    // vector<int> max_Track;

public:
    // TAJIMA GENE FILE
    /**
     * Tajima GENE MODE CONSTRUCTOR
     **/
    prometheus(vector<pair<string, string>> folder_Index, string Multi_read, int tot_Blocks, int tot_ThreadsperBlock, long int combinations, float an, float e1, float e2, int N, int CPU_cores, int SNPs_per_Run, int number_of_genes);
    // TAJIMA WINDOW
    /**
     * Tajima WINDOW MODE CONSTRUCTOR
     **/
    prometheus(string output_File, int window_Size, int step_Size, vector<pair<string, string>> folder_Index, string Multi_read, int tot_Blocks, int tot_ThreadsperBlock, long int combinations, float an, float e1, float e2, int N, int CPU_cores, int SNPs_per_Run, int number_of_genes);

    // FU_LI
    /**
     * FU_LI GENE MODE CONSTRUCTOR
     **/
    prometheus(vector<pair<string, string>> folder_Index, string Multi_read, int tot_Blocks, int tot_ThreadsperBlock, int CPU_cores, int SNPs_per_Run, int number_of_genes, int N, long int combinations, float an, float vd, float ud, float vd_star, float ud_star, float uf, float vf, float uf_star, float vf_star);
    // FU LI WINDOW
    /**
     * FU LI WINDOW MODE CONSTRUCTOR
     **/
    prometheus(string output_File, int window_Size, int step_Size, vector<pair<string, string>> folder_Index, string Multi_read, int tot_Blocks, int tot_ThreadsperBlock, int CPU_cores, int SNPs_per_Run, int number_of_genes, int N, long int combinations, float an, float vd, float ud, float vd_star, float ud_star, float uf, float vf, float uf_star, float vf_star);

    // FAY_WU
    /**
     * FAY WU GENE MODE CONSTRUCTOR
     **/
    prometheus(vector<pair<string, string>> folder_Index, string Multi_read, int tot_Blocks, int tot_ThreadsperBlock, int CPU_cores, int SNPs_per_Run, int number_of_genes, int N, long int combinations, float an, float bn, float bn_plus1);
    // FAY_WU WINDOW
    /**
     * FAY WU WINDOW MODE CONSTRUCTOR
     **/
    prometheus(string output_File, int window_Size, int step_Size, vector<pair<string, string>> folder_Index, string Multi_read, int tot_Blocks, int tot_ThreadsperBlock, int CPU_cores, int SNPs_per_Run, int number_of_genes, int N, long int combinations, float an, float bn, float bn_plus1);

    // NEUTRALITY ALL
    /**
     * NEUTRALITY GENE MODE CONSTRUCTOR
     **/
    prometheus(vector<pair<string, string>> folder_Index, string Multi_read, int tot_Blocks, int tot_ThreadsperBlock, int CPU_cores, int SNPs_per_Run, int number_of_genes, int N, long int combinations, float an, float e1, float e2, float vd, float ud, float vd_star, float ud_star, float uf, float vf, float uf_star, float vf_star, float bn, float bn_plus1);
    // NEUTRALITY WINDOW
    /**
     * NEUTRALITY WINDOW MODE CONSTRUCTOR
     **/
    prometheus(string output_File, int window_Size, int step_Size, vector<pair<string, string>> folder_Index, string Multi_read, int tot_Blocks, int tot_ThreadsperBlock, int CPU_cores, int SNPs_per_Run, int number_of_genes, int N, long int combinations, float an, float e1, float e2, float vd, float ud, float vd_star, float ud_star, float uf, float vf, float uf_star, float vf_star, float bn, float bn_plus1);

    // vector<string> collection_Engine(vector<string> &gene_Collect);

    /**
     * Responsible for the processing of the collected query regions in gene mode.
     **/
    vector<string> collection_Engine(vector<string> &gene_Collect, string test_Type);
    /**
     * Responsible for processing Window analysis.
     **/
    void process_Window(string test);
    /**
     * Responsible for processing Sliding Window analysis.
     **/
    void process_C_sliding_Window(string test);

    /**
     * In Sliding Window mode gets the position of the respective SNP and moves forward from that point.
     **/
    void get_POS_VALID(int start_Co, int gene_ID);

    /**
     * Extracts the relevant query region information as well as the catch point for the CIS search.
     **/
    void get_Gene_info_and_catch(string gene_Combo, int gene_ID);
    /**
     * CIS forward in space.
     **/
    void forward_Search(int pos, int start_Co, int end_Co, int gene_ID);
    /**
     * CIS backward in space.
     **/
    void backward_Search(int pos, int start_Co, int end_Co, int gene_ID);
    /**
     * Compiles the three CIS functions results into one.
     **/
    void compile(int gene_ID);

    // vector<string> compound_interpolationSearch(int &start_Co, int &end_Co);

    // Multi_read
    /**
     * Reads the file segments in parallel.
     * * For SSDs.
     **/
    void file_Reader_multi(string files);
    // Single_read
    /**
     * Reads the file segments in serial.
     * * For HDDs.
     **/
    void file_Reader_single();
    // void set_values(int gene_ID, string gene_Name, string chromosome, string start_Co, string end_Co);

    /**
     * All four functions below are administrative functions that process the segregating sites and calculates the relevant tests statistics.
     **/
    void process_Tajima();
    void process_Fu_Li();
    void process_Fay_Wu();
    void process_Neutrality();

    /**
     * Concatenates the Segregating sites for processing by each GPU round.
     **/
    void seg_Concat(int round_ID, int start_Seg, int stop_Seg);
    // void seg_Concat_New(int round_ID, int start_Seg, int stop_Seg);

    /**
     * Part of the Compound Binary Search (CBS) used to find the relevant segregating sites by position to process eah query region.
     **/
    void seg_Search_catch_point(int gene_ID);
    void seg_Search_forward(int gene_ID);
    void seg_Search_backward(int gene_ID);

    /**
     * Extracts the POsition and organizes the processed segregating sites by position and their position in the overall data store.
     **/
    void seg_Indexer(int start_Seg, int stop_Seg, char *full_Char, int *VALID_or_NOT, int *pos_start_Index, int *pos_end_Index, int start);

    /**
     * Following four functions conducts the final calculations of the test statistic.
     **/
    void calc_Tajima_Segs(int gene_ID, int *MA_Count);
    void calc_Fu_Li_Segs(int gene_ID, int *MA_Count, int *ne, int *ns);
    void calc_Fay_Wu_Segs(int gene_ID, int *MA_Count, int *Theta_partials);
    void calc_Neutrality_Segs(int gene_ID, int *MA_Count, int *ne, int *ns, int *Theta_partials);

    /**
     * Clears all relevant variables for the next round of Prometheus.
     * Acts like a destructor but for each round of Prometheus execution.
     **/
    void erase();

    void intialize(vector<pair<int, int>> &coordinates);

    /**
     * These functions is used in conjunction with the constructors to set the common private variables.
     **/
    void set_Values(vector<pair<string, string>> folder_Index, int tot_Blocks, int tot_ThreadsperBlock, string Multi_read, int CPU_cores, int SNPs_per_Run, int number_of_genes);
    void set_Values_Window(string output_File, int window_Size, int step_Size, int number_of_genes_Window);
};