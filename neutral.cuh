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
private:
    string calc_Mode = "FILE";
    int window_Size = 0;
    int step_Size = 0;

    string gene_List;
    string input_Folder;
    string output_Path;
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
    neutral(string gene_List, string input_Folder, string output_Path, int cuda_ID, string intermediate_Path, int ploidy);
    // PROMETHEUS GENE mode
    neutral(string gene_List, string input_Folder, string output_Path, int cuda_ID, string intermediate_Path, int ploidy, string prometheus_Activate, string Multi_read, int number_of_genes, int CPU_cores, int SNPs_per_Run);
    // WINDOW mode PROMETHEUS
    neutral(string calc_Mode, int window_Size, int step_Size, string input_Folder, string output_Path, int cuda_ID, int ploidy, string prometheus_Activate, string Multi_read, int number_of_genes, int CPU_cores, int SNPs_per_Run);

    void set_Values(string gene_List, string input_Folder, string ouput_Path, int cuda_ID, string intermediate_Path, int ploidy);

    void ingress();

    void get_Prerequisites(int N_tot, float &an, float &e1, float &e2, float &vd, float &ud, float &vd_star, float &ud_star, float &uf, float &vf, float &uf_star, float &vf_star, float &bn, float &bn_plus1);
    void pre_Tajima(float N_tot, float an, float bn, promise<float> &e1_set, promise<float> &e2_set);
    void pre_Fu_li(float N_float_tot, float N_float, float an, float bn, promise<float> &vd_set, promise<float> &ud_set, promise<float> &vd_star_set, promise<float> &ud_star_set, promise<float> &uf_set, promise<float> &vf_set, promise<float> &uf_star_set, promise<float> &vf_star_set);

    void process_Segs(vector<string> &total_Segregrating_sites, float N, int &real_segregrating_Sites, float &tot_pairwise_Differences, int &singletons_ne, int &singletons_ns, int &Total_iEi, float &theta_L, int tot_Blocks, int tot_ThreadsperBlock);

    void calculate_Neutrality(float N_float, float &pi, int &segregating_Sites, float &an, float &bn, float &e1, float &e2, int &singletons_ne, int &singletons_ns, float &vd, float &ud, float &vd_star, float &ud_star, float &vf, float &uf, float &vf_star, float &uf_star, float &theta_L, float &bn_plus1, string &Tajima_D, string &Fu_Li_D, string &Fu_Li_D_star, string &Fu_Li_F, string &Fu_Li_F_star, string &Fay_Wu_H, string &Fay_Wu_E);
    void Tajimas_D_thread(promise<string> &Tajimas_D_value, float pi, int segregating_Sites, float an, float e1, float e2);
    void Fu_li_thread(promise<string> &Fu_Li_D, promise<string> &Fu_Li_D_star, promise<string> &Fu_Li_F, promise<string> &Fu_Li_F_star, float N_float, float pi, int segregating_Sites, float an, int singletons_ne, int singletons_ns, float vd, float ud, float vd_star, float ud_star, float vf, float uf, float vf_star, float uf_star);
    void Fay_wu_thread(promise<string> &Fay_Wu_H, promise<string> &Fay_Wu_E, float N_float, float pi, int segregating_Sites, float an, float bn, float bn_plus1, float theta_L);
};