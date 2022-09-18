#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <filesystem>
#include <list>
#include <set>
#include <sstream>

using namespace std;

class gene_extract
{
    /**
     * Extracts gene sequences using from a reference sequence file.
     **/

private:
    /**
     * @param gene_List defines the CATE compatible GENE list file. It provides the base pair coordinates of the genes to be extracted.
     **/
    string gene_List;

    /**
     * @param reference_Path defines the path of the reference FASTA file.
     **/
    string reference_Path;

    /**
     * @param output_FASTA defines the output folder path.
     **/
    string output_Path;

    /**
     * @param intermediate_Path defines the intermediate folder path.
     **/
    string intermediate_Path;

public:
    /**
     * Constructor Function assigns passed variables to the classes' private variable.
     **/
    gene_extract(string gene_List, string reference_Path, string output_Path, string intermediate_Path);

    /**
     * Execution function.
     * Extracts the genes in FASTA formats.
     **/
    void ingress();
};