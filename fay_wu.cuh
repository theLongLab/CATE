#pragma once
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <cstring>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

using namespace std;

class fay_wu
{

private:
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
    fay_wu(string gene_List, string input_Folder, string ouput_Path, int cuda_ID, string intermediate_Path, int ploidy);
    fay_wu(string gene_List, string input_Folder, string ouput_Path, int cuda_ID, string intermediate_Path, int ploidy, string prometheus_Activate, string Multi_read, int number_of_genes, int CPU_cores, int SNPs_per_Run);
    void set_Values(string gene_List, string input_Folder, string ouput_Path, int cuda_ID, string intermediate_Path, int ploidy);

    void ingress();
    void calc_Pre(float &an, float &bn, float &bn_plus1, int N_tot);
    float calc_theta_L(vector<string> &total_Segregrating_sites, float N_tot, int &num_segregrating_Sites, int &Total_iEi, float &tot_pairwise_Differences);
};