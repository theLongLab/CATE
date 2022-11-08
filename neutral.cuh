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

class neutral
{
    /**
     * Execution of all three Neutrality tests in unison.
     * Execution of GENE mode and WINDOW mode.
     * Execution of NORMAL mode and PROMETHEUS mode.
     **/
private:
    /**
     * @param calc_Mode acts as a boolean variable to determine between GENE(FILE) mode or WINDOW mode.
     **/
    string calc_Mode = "FILE";
    /**
     * @param window_Size defines the base pair (bp) size of the query window for analysis. Used with WINDOW mode.
     **/
    int window_Size = 0;
    /**
     * @param step_Size defines the number of bps to skip to reach the next query window. Used with WINDOW mode.
     **/
    int step_Size = 0;

    /**
     * @param gene_List defines the CATE compatible GENE list file. Used with GENE mode.
     **/
    string gene_List;

    /**
     * @param input_Folder defines the path of the CATE indexed VCF folder.
     **/
    string input_Folder;
    /**
     * @param ouput_Path defines the path of the directory to which the outputs will be written to.
     **/
    string output_Path;
    /**
     * @param intermediate_Path defines the intermediate folder path.
     **/
    string intermediate_Path;

    /**
     * @param ploidy defines number of sets of chromosomes per individual organism.
     **/
    int ploidy;

    /**
     * @param tot_Blocks defines number of GPU blocks that are available
     **/
    int tot_Blocks;
    /**
     * @param tot_ThreadsperBlock defines number of threads that are available per GPU block.
     **/
    int tot_ThreadsperBlock;

    /**
     * @param prometheus_Activate acts as a boolean variable and is used to activate Prometheus.
     **/
    string prometheus_Activate = "NO";
    /**
     * @param Multi_read acts as a boolean variable and is used to enable SSD based multiple reading of multiple files at once.
     **/
    string Multi_read = "NO";
    /**
     * @param number_of_genes defines the number of query regions Prometheus should handle at once.
     **/
    int number_of_genes = 0;
    /**
     * @param CPU_cores defines the number of CPU cores that are available at Prometheus's disposal.
     **/
    int CPU_cores = 0;
    /**
     * @param SNPs_per_Run defines the number of SNP sites the GPU will process at a time.
     **/
    int SNPs_per_Run = 0;

public:
    /**
     * Fay and Wu's normalized H and E has 4 constructor Functions.
     * They are for the 4 separate modal combinations that Neutrality tests can be configured to. Namely:
     * 1. NORMAL - GENE MODE
     * 2. NORMAL - WINDOW MODE
     * 3. PROMETHEUS - GENE MODE
     * 4. PROMETHEUS - WINDOW MODE
     * All serve the common function of assigning passed variables to the classes' private variable.
     **/

    /**
     * NORMAL - GENE MODE constructor
     **/
    neutral(string gene_List, string input_Folder, string output_Path, int cuda_ID, string intermediate_Path, int ploidy);
    /**
     * NORMAL - WINDOW MODE constructor
     **/
    neutral(string calc_Mode, int window_Size, int step_Size, string input_Folder, string output_Path, int cuda_ID, int ploidy);
    /**
     * PROMETHEUS - GENE MODE constructor
     **/
    neutral(string gene_List, string input_Folder, string output_Path, int cuda_ID, string intermediate_Path, int ploidy, string prometheus_Activate, string Multi_read, int number_of_genes, int CPU_cores, int SNPs_per_Run);
    /**
     * PROMETHEUS - WINDOW MODE constructor
     **/
    neutral(string calc_Mode, int window_Size, int step_Size, string input_Folder, string output_Path, int cuda_ID, int ploidy, string prometheus_Activate, string Multi_read, int number_of_genes, int CPU_cores, int SNPs_per_Run);

    /**
     * This function is used in conjunction with the constructors to set the common private variables.
     **/
    void set_Values(string gene_List, string input_Folder, string ouput_Path, int cuda_ID, string intermediate_Path, int ploidy);

    /**
     * Execution function.
     **/
    void ingress();

    /**
     * This is the normal WINDOW FUNCTION. Where there is a fixed step size.
     **/
    void window(string output_File, float an, float e1, float e2, float vd, float ud, float vd_star, float ud_star, float uf, float vf, float uf_star, float vf_star, float bn, float bn_plus1, float N_float, long int combinations, vector<pair<string, string>> &folder_Index);
    /**
     * This is the sliding WINDOW FUNCTION. Where there is a not a fixed step size. But the WINDOW slides from one SNP to the next in the VCF file.
     * This is activated by setting the step_Size to 0.
     **/
    void window_Sliding(string output_File, float an, float e1, float e2, float vd, float ud, float vd_star, float ud_star, float uf, float vf, float uf_star, float vf_star, float bn, float bn_plus1, float N_float, long int combinations, vector<pair<string, string>> &folder_Index);

    /**
     * These 3 functions are used to calculate the prerequisites using the sample size required for calculating the Neutrality tests.
     **/
    void get_Prerequisites(int N_tot, float &an, float &e1, float &e2, float &vd, float &ud, float &vd_star, float &ud_star, float &uf, float &vf, float &uf_star, float &vf_star, float &bn, float &bn_plus1);
    void pre_Tajima(float N_tot, float an, float bn, promise<float> &e1_set, promise<float> &e2_set);
    void pre_Fu_li(float N_float_tot, float N_float, float an, float bn, promise<float> &vd_set, promise<float> &ud_set, promise<float> &vd_star_set, promise<float> &ud_star_set, promise<float> &uf_set, promise<float> &vf_set, promise<float> &uf_star_set, promise<float> &vf_star_set);

    /**
     * Administrative function responsible for calculating the values needed for the Neutrality tests.
     * Function directly calls upon the GPU function.
     **/
    void process_Segs(vector<string> &total_Segregrating_sites, float N, int &real_segregrating_Sites, float &tot_pairwise_Differences, int &singletons_ne, int &singletons_ns, int &Total_iEi, float &theta_L, int tot_Blocks, int tot_ThreadsperBlock);

    /**
     * These functions Calculate the three Neutrality tests values using the processed data from the segregating sites.
     * The calculate_Neutrality is the administrative function. Calls the other three functions in parallel.
     **/
    void calculate_Neutrality(float N_float, float &pi, int &segregating_Sites, float &an, float &bn, float &e1, float &e2, int &singletons_ne, int &singletons_ns, float &vd, float &ud, float &vd_star, float &ud_star, float &vf, float &uf, float &vf_star, float &uf_star, float &theta_L, float &bn_plus1, string &Tajima_D, string &Fu_Li_D, string &Fu_Li_D_star, string &Fu_Li_F, string &Fu_Li_F_star, string &Fay_Wu_H, string &Fay_Wu_E);
    void Tajimas_D_thread(promise<string> &Tajimas_D_value, float pi, int segregating_Sites, float an, float e1, float e2);
    void Fu_li_thread(promise<string> &Fu_Li_D, promise<string> &Fu_Li_D_star, promise<string> &Fu_Li_F, promise<string> &Fu_Li_F_star, float N_float, float pi, int segregating_Sites, float an, int singletons_ne, int singletons_ns, float vd, float ud, float vd_star, float ud_star, float vf, float uf, float vf_star, float uf_star);
    void Fay_wu_thread(promise<string> &Fay_Wu_H, promise<string> &Fay_Wu_E, float N_float, float pi, int segregating_Sites, float an, float bn, float bn_plus1, float theta_L);
};