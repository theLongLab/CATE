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

public:
    fay_wu(string gene_List, string input_Folder, string ouput_Path, int cuda_ID, string intermediate_Path, int ploidy);

    void ingress();
    void calc_Pre(float &an, float &bn, float &bn_plus1, int N_tot);
    float calc_theta_L(vector<string> &total_Segregrating_sites, float N_tot, int &num_segregrating_Sites,int &Total_iEi,float &tot_pairwise_Differences);
};