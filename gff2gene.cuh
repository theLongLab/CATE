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

class gff2gene
{
    /**
     * Collects and writes all GENE's in a *.GFF file to a CATE format gene list file.
     **/
    
private:
    /**
     * @param input_File defines the path of the gff file.
     **/
    string input_File;

    /**
     * @param output_Path defines the path of the outut folder to which the gene file will be written to.
     **/
    string output_Path;

public:
    /**
     * Constructor Function assigns passed variables to the classes' private variable.
     **/
    gff2gene(string input_File, string output_Path);

    /**
     * Execution function.
     * Generates CATE's gene file from a GFF file.
     **/
    void ingress();
};