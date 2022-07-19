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

class fu_li
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
    // Normal GENE mode
    fu_li(string gene_List, string input_Folder, string ouput_Path, int cuda_ID, string intermediate_Path, int ploidy);
    // PROMETHEUS GENE mode
    fu_li(string gene_List, string input_Folder, string ouput_Path, int cuda_ID, string intermediate_Path, int ploidy, string prometheus_Activate, string Multi_read, int number_of_genes, int CPU_cores, int SNPs_per_Run);
    // WINDOW mode PROMETHEUS
    fu_li(string calc_Mode, int window_Size, int step_Size, string input_Folder, string ouput_Path, int cuda_ID, int ploidy, string prometheus_Activate, string Multi_read, int number_of_genes, int CPU_cores, int SNPs_per_Run);

    void set_Values(string gene_List, string input_Folder, string ouput_Path, int cuda_ID, string intermediate_Path, int ploidy);

    void ingress();
    void calc_Pre(int N_tot, float &an, float &vd, float &ud, float &vd_star, float &ud_star, float &uf, float &vf, float &uf_star, float &vf_star);
    int outgroup_Singleton(vector<string> &info, vector<string> &positions);
};