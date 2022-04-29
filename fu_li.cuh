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
    string gene_List;
    string input_Folder;
    string ouput_Path;
    string intermediate_Path;

    int ploidy;

    int tot_Blocks;
    int tot_ThreadsperBlock;

public:
    fu_li(string gene_List, string input_Folder, string ouput_Path, int cuda_ID, string intermediate_Path, int ploidy);

    void ingress();
    void calc_Pre(int N_tot, float &an, float &vd, float &ud, float &vd_star, float &ud_star, float &uf, float &vf, float &uf_star, float &vf_star);
    int outgroup_Singleton(vector<string> &info, vector<string> &positions);

};