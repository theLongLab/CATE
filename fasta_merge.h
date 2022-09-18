#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <filesystem>
#include <list>
#include <set>
#include <sstream>
#include <vector>
#include <algorithm>

using namespace std;

class fasta_merge
{
    /**
     * Can split a merged FASTA file into its individual sequences.
     **/

private:
    /**
     * @param fasta_Folder defines the FASTA folder with the multiple sequences.
     **/
    string fasta_Folder;

    /**
     * @param output_FASTA defines the output folder to which the resultant merged FASTA sequence will be written to.
     **/
    string output_FASTA;

public:
    /**
     * Constructor Function assigns passed variables to the classes' private variable.
     **/
    fasta_merge(string fasta_Folder, string output_FASTA);

    /**
     * Execution function.
     **/
    void ingress();

    /**
     * Collects all the FASTA compatible files present in the fasta_Folder folder.
     **/
    vector<string> index_FASTA_folder();
};