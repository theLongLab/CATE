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

using namespace std;

class fasta_splitter
{
    /**
     * Splits a merged FASTA file into its individual sequences.
     * It can also be used to extract a specific sequence by its sequence ID.
     **/

private:
    /**
     * @param fasta_File defines the merged FASTA file with the multiple sequences.
     **/
    string fasta_File;

    /**
     * @param output_Folder defines the output folder to which all sequences will be written to in the FASTA format.
     **/
    string output_Folder;

    /**
     * @param fasta_Name the sequence ID of the sequence to be extracted.
     * If it is given as "all" then all sequences will be extracted.
     **/
    string fasta_Name;

public:
    /**
     * Constructor Function assigns passed variables to the classes' private variable.
     **/
    fasta_splitter(string fasta_File, string output_Path, string fasta_Name);

    /**
     * This function is used to extract all the sequences from the merged file.
     **/
    void split_all();

    /**
     * This function is used to extract a specific sequence from the merged file.
     **/
    void split_select();

    /**
     * Execution function.
     **/
    void ingress();
};