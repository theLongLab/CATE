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
// #include <mutex>
// #include <shared_mutex>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

using namespace std;

class hap_extract
{
    /**
     * Extracts the unique haplotypes in a user specified region.
     * Options are available for extracting only haplotype information, including sequence extraction,
     * if needed the complete list of sequences from the population can be extracted.
     **/

private:
    /**
     * @param input_Folder defines the path of the CATE indexed VCF folder.
     **/
    string input_Folder;

    /**
     * @param output_Path defines the path of the directory to which the outputs will be written to.
     **/
    string output_Path;

    /**
     * @param intermediate_Path defines the intermediate folder path.
     **/
    string intermediate_Path;

    /**
     * @param gene_List defines the CATE compatible GENE list file.
     **/
    string gene_List;

    /**
     * @param reference_Path defines the path of the reference FASTA file.
     * This has to be the same reference file from which the parent VCF file was created with.
     **/
    string reference_File;

    /**
     * @param pop_Out acts as a boolean variable to determine whether the entire population of FASTA files should be reconstructed.
     * If "NO" then only the haplotype sequences will be generated.
     **/
    string pop_Out = "NO";

    /**
     * @param ploidy defines number of sets of chromosomes per individual organism.
     **/
    int ploidy;
    /**
     * @param N defines number of total sequences being present per SNP.
     **/
    int N;

    /**
     * @param tot_Blocks defines number of GPU blocks that are available
     **/
    int tot_Blocks;
    /**
     * @param tot_ThreadsperBlock defines number of threads that are available per GPU block.
     **/
    int tot_ThreadsperBlock;

    /**
     * @param cuda_reference used to capture and store the reference genome sequence in GPU memory.
     **/
    char *cuda_reference;
    // int reference_size;

public:
    /**
     * Constructor Function assigns passed variables to the classes' private variable.
     **/
    hap_extract(string gene_List, string input_Folder, string output_Path, int cuda_ID, string intermediate_Path, int ploidy, string reference_File, string pop_Out);
    /**
     * This function is used in conjunction with the class constructor to set the common private variables.
     **/
    void set_Values(string gene_List, string input_Folder, string output_Path, int cuda_ID, string intermediate_Path, int ploidy);

    /**
     * Execution function.
     **/
    void ingress();

    /**
     * Administrative function responsible for haplotype reconstruction.
     * Function directly calls upon the GPU function.
     **/
    void hap_extraction(vector<string> &write_Lines, vector<string> &write_Sequences, vector<string> &total_Segregrating_sites, vector<pair<int, int>> &pos_INDEX, string gene_Name, string chr, int start_Pos, int end_Pos);
};