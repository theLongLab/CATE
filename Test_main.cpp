#include <iostream>
#include <filesystem>
#include <unistd.h>
#include "test.cuh"
#include "cudaDevices.cuh"
#include "tajima.cuh"
#include "fu_li.cuh"
#include "fay_wu.cuh"
#include "ehh.cuh"
#include "vcf_splitter.h"
#include "mk_test.cuh"
#include "fst.cuh"
#include "fasta_splitter.h"
#include "parameter.h"
#include "fasta_merge.h"
#include <bits/stdc++.h>
#include <cuda.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

using namespace std;

void print_Cate();
void print_HELP();

int main(int argc, char *argv[])
{
     cout.precision(10);

     print_Cate();
     cout << "CATE: CUDA Accelerated Testing of Evolution" << endl;
     cout << "Evolutionary tests for large scale genome data" << endl
          << "----------------------------------------------" << endl
          << endl;

     if (argv[1] != NULL)
     {
          string function(argv[1]);
          transform(function.begin(), function.end(), function.begin(), ::tolower);

          if (function == "--help" || function == "-h")
          {
               print_HELP();
          }
          else if (function == "--godmode")
          {
               test try_something = test();
               try_something.run();
          }
          else if (function == "--cuda" || function == "-c")
          {
               cudaDevices cudaList = cudaDevices();
               cout << "All CUDA capable devices have been listed" << endl;
          }
          else
          {

               if (argv[2] != NULL)
               {
                    string parameter_File(argv[2]);

                    char tmp[256];
                    getcwd(tmp, 256);

                    parameter properties = parameter(parameter_File);

                    string output_Path = properties.where("Output path");
                    if (filesystem::exists(output_Path) == 0)
                    {
                         cout << "Creating output folder: " << output_Path << endl;

                         filesystem::create_directory(output_Path);
                    }
                    else
                    {
                         cout << "Output folder exists: " << output_Path << endl;
                    }

                    string intermediate_Path = properties.where("Intermediate path");
                    if (filesystem::exists(intermediate_Path) == 0)
                    {
                         cout << "Creating intermediate folder: " << intermediate_Path << endl;

                         filesystem::create_directory(intermediate_Path);
                    }
                    else
                    {
                         cout << "Intermediate folder exists: " << intermediate_Path << endl;
                    }

                    cout << endl;

                    if (function == "--splitvcf" || function == "-svcf")
                    {
                         // string input = properties.where("Input folder entry");

                         vcf_splitter split = vcf_splitter(tmp, properties.where("Input path"), properties.where("Population file path"), output_Path, properties.where_Int("Reference allele count"), properties.where_Int("Alternate allele count"), properties.where_Int("SNP count per file"));
                         split.index_population();
                         split.read_File();
                         cout << endl
                              << endl
                              << "VCF split has been completed." << endl;
                    }
                    else if (function == "--splitfasta" || function == "-sfasta")
                    {
                         fasta_splitter split = fasta_splitter(properties.where("Raw FASTA file"), output_Path, properties.where("Sequence"));
                         split.ingress();

                         cout << endl
                              << "FASTA split has been completed." << endl;
                    }
                    else if (function == "--mergefasta" || function == "-mfasta")
                    {
                         fasta_merge merge = fasta_merge(properties.where("FASTA files folder"), properties.where("Merge FASTA path"));
                         merge.ingress();

                         cout << endl
                              << "FASTA merge has been completed." << endl;
                    }
                    else if (function == "--tajima" || function == "-t")
                    {
                         string gene_List = properties.where("Tajima gene list");
                         if (gene_List == "universal")
                         {
                              gene_List = properties.where("Universal gene list");
                         }
                         tajima tajimasD = tajima(gene_List, properties.where("Input path"), output_Path, properties.where_Int("CUDA Device ID"), intermediate_Path, properties.where_Int("Ploidy"));
                         tajimasD.ingress();

                         cout << endl
                              << "CUDA powered Tajima's D calculator has been completed." << endl;
                    }

                    else if (function == "--fuli" || function == "-f")
                    {
                         string gene_List = properties.where("Fu and Li gene list");
                         if (gene_List == "universal")
                         {
                              gene_List = properties.where("Universal gene list");
                         }
                         fu_li fuli = fu_li(gene_List, properties.where("Input path"), output_Path, properties.where_Int("CUDA Device ID"), intermediate_Path, properties.where_Int("Ploidy"));
                         fuli.ingress();

                         cout << "CUDA powered Fu and Li's D, D*, F and F* calculator has been completed." << endl;
                    }

                    else if (function == "--faywu" || function == "-w")
                    {
                         string gene_List = properties.where("Fay and Wu gene list");
                         if (gene_List == "universal")
                         {
                              gene_List = properties.where("Universal gene list");
                         }
                         fay_wu faywu = fay_wu(gene_List, properties.where("Input path"), output_Path, properties.where_Int("CUDA Device ID"), intermediate_Path, properties.where_Int("Ploidy"));
                         faywu.ingress();

                         cout << "CUDA powered Fay and Wu's normalized H and E calculator has been completed." << endl;
                    }
                    else if (function == "--mk" || function == "-m")
                    {
                         string gene_List = properties.where("McDonald–Kreitman gene list");
                         if (gene_List == "universal")
                         {
                              gene_List = properties.where("Universal gene list");
                         }
                         mk_test mk = mk_test(properties.where("Reference genome"), properties.where("Alignment file"), gene_List, properties.where("Input path"), output_Path, properties.where_Int("CUDA Device ID"), intermediate_Path, properties.where_Int("Ploidy"), properties.where("Genetic code"), properties.where("Start codon(s)"), properties.where("Stop codon(s)"));
                         mk.ingress();

                         cout << "CUDA powered McDonald–Kreitman Neutrality Index (NI) test has been completed." << endl;
                    }
                    else if (function == "--fst" || function == "-x")
                    {
                         string gene_List = properties.where("Fst gene list");
                         if (gene_List == "universal")
                         {
                              gene_List = properties.where("Universal gene list");
                         }

                         fst fs = fst(gene_List, properties.where("Input path"), output_Path, properties.where_Int("CUDA Device ID"), intermediate_Path, properties.where_Int("Ploidy"), properties.where("Population index file path"), properties.where("Population ID"));
                         fs.ingress();

                         cout << "CUDA powered Fst (Fixation Index) calculator has been completed." << endl;
                    }
                    else if (function == "--ehh" || function == "-e")
                    {
                         string mode = properties.where("Range mode");
                         transform(mode.begin(), mode.end(), mode.begin(), ::toupper);

                         // cout << mode << endl;

                         string file_mode_Path = "NA";
                         string fixed_mode_Value = "NA";

                         string GO = "NO";

                         if (mode == "FILE")
                         {
                              file_mode_Path = properties.where("EHH FILE path");
                              if (file_mode_Path == "universal")
                              {
                                   file_mode_Path = properties.where("Universal gene list");
                              }
                              GO = "YES";
                         }
                         else if (mode == "FIXED")
                         {
                              file_mode_Path = properties.where("EHH FILE path");
                              if (file_mode_Path == "universal")
                              {
                                   file_mode_Path = properties.where("Universal gene list");
                              }
                              fixed_mode_Value = properties.where("FIXED mode");
                              GO = "YES";
                         }

                         if (GO == "YES")
                         {
                              ehh ehh_ = ehh(mode, file_mode_Path, fixed_mode_Value, properties.where("Input path"), output_Path, properties.where_Int("CUDA Device ID"), intermediate_Path, properties.where_Int("Ploidy"));
                              ehh_.ingress();

                              cout << "CUDA powered Extended Haplotype Homozygosity (EHH) calculator has been completed." << endl;
                         }
                         else
                         {
                              cout << "RANGE MODE HAS AN INCORRECT VALUE. Should be either \"FILE\" or \"FIXED\"" << endl;
                         }
                    }
                    else
                    {
                         cout << "INVALID FUNCTION PASSED AS A RUNTIME ARGUMENT" << endl
                              << "To see available functions type --help or -h." << endl;
                    }
               }
               else
               {
                    cout << "PLEASE SPECIFY A PARAMETER JSON FILE. PROGRAM REQUIRES A PARAMETER FILE TO EXECUTE.";
               }
          }
     }
     else
     {
          cout << "PLEASE ENTER A FUNCTION. ENTER -h OR --help TO LIST ALL AVAILABLE FUNCTIONS";
     }

     cout << endl
          << "Program has completed its run." << endl;

     return 0;
}

void print_Cate()
{
     cout << "   0000000    00000   0000000   0000000" << endl;
     cout << "  0000       00   00    000     00" << endl;
     cout << " 000        00     00   000     00" << endl;
     cout << "0           000000000   000     0000000" << endl;
     cout << " 000        00     00   000     00" << endl;
     cout << "  0000      00     00   000     00" << endl;
     cout << "   0000000  00     00   000     0000000" << endl;
     cout << endl;
}

void print_HELP()
{
     cout << "HELP MENU" << endl
          << "---------"
          << endl
          << endl
          << "Excecution format: \"[--function or -f] properties_file.json\"" << endl
          << endl
          << "** Available functions are (not CaSe sensitive) **" << endl
          << endl
          << "TOOLS:"
          << endl
          << "A set of simple tools used to manipulate and alter VCF and FASTA files." << endl
          << endl
          << "-svcf or --splitvcf\t: Splits the 1000 Genome's VCF file based on population and a series of pre specified parameters. Mainly used as an indexer." << endl
          << "            \t\t  Split files are placed in their respective population folders." << endl
          << "            \t\t  Split files are named as follows: \"CHROMOSOMEnumber_COUNTRY_STARTposition_ENDposition.vcf\"." << endl
          << endl
          << "-sfasta or --splitfasta\t: Split a user specified FASTA file to individual FASTA files." << endl
          << "             \t\t  Can be used to extract a singular user specified sequence as well." << endl
          << "             \t\t  Split files are placed in a user specified folder." << endl
          << "             \t\t  Each FASTA file name will be the name of the respectve sequence entry." << endl
          << endl
          << "-mfasta or --mergefasta\t: Merge all FASTA files in a user specified folder to an individual FASTA file." << endl
          << "             \t\t  Ensure that the FASTA files have the APPROPRIATE extensions: .fasta, .fna, .ffn, .faa, .frn, .fa" << endl
          << endl
          << endl
          << "NEUTRALITY TESTS:"
          << endl
          << "Core functions optimized for conducting Neutrality tests on 1000 Genome VCF files." << endl
          << endl
          << "-t or --tajima\t: Calculates the Tajima's D statistic (1989) for a predefined gene list using a (split) VCF (indexed) folder." << endl
          << "              \t  Uses a CUDA powered engine, therefore, requires a CUDA capable GPU." << endl
          << "              \t  File format is *.td (a tab deliminated text file)." << endl
          << endl
          << "-f or --fuli\t: Calculates the Fu and Li's D, D*, F and F* statistics (1993) for a predefined gene list using a (split) VCF (indexed) folder." << endl
          << "              \t  ** The D, D* and F statistics are calculated based on the original paper by Fu et al (1993)." << endl
          << "              \t  ** The F* statistic's vf* and uf* are calculated based on the corrected equations in Simonsen et al (1995)." << endl
          << "              \t  Uses a CUDA powered engine, therefore, requires a CUDA capable GPU." << endl
          << "              \t  File format is *.fl (a tab deliminated text file)." << endl
          << endl
          << "-w or --faywu\t: Calculates the Fay and Wu's normalized H and E statistics (2006) for a predefined gene list using a (split) VCF (indexed) folder." << endl
          << "              \t  Uses a CUDA powered engine, therefore, requires a CUDA capable GPU." << endl
          << "              \t  File format is *.fw (a tab deliminated text file)." << endl
          << endl
          << "-m or --mk    \t: Calculates the McDonald–Kreitman Neutrality Index (NI) (1991) for a predefined gene list using a (split) vcf (indexed) folder." << endl
          << "              \t  The reference genome must be provided in a FASTA format file." << endl
          << "              \t  Alignment file of the reference genome to the outgroup genome must also be provided in a *.maf format file." << endl
          << "              \t  PLEASE ensure that the REFERENCE sequence is first and QUERY sequence is second in the MAF file." << endl
          << "              \t  ** TIP: Chromosome wide whole genome alignment software: GSAlign (https://github.com/hsinnan75/GSAlign)." << endl
          << "              \t  Uses a CUDA powered engine, therefore, requires a CUDA capable GPU." << endl
          << "              \t  File format is *.mc (a tab deliminated text file)." << endl
          << endl
          << "-x or --fst   \t: Calculates the Fixation Index (Fst) (1965) for a predefined gene list using a (split) vcf (indexed) folder." << endl
          << "              \t  The population index of the seqeunced samples must be provided in a tab deliminated *txt format file." << endl
          << "              \t  Uses a CUDA powered engine, therefore, requires a CUDA capable GPU." << endl
          << "              \t  File format is *.fst (a tab deliminated text file)." << endl
          << endl
          << "-e or --ehh   \t: Calculates the Extended Haplotype Homozygosity (EHH) (2002) for a predefined gene list using a (split) vcf (indexed) folder." << endl
          << "              \t  The \"MODE\" used to generate the extended haplotype region must be specified. Either \"FILE\" or \"FIXED\" mode." << endl
          << "              \t  FILE mode: In the tab deliminated gene list file a tertiary column containing the extended regions dimension's will be present." << endl
          << "              \t             Formats include \"START_position:END_position\" or +VALUE or -VALUE." << endl
          << "              \t             \"+\" Causes the START_position of the gene region to be incremented by the user specified value." << endl
          << "              \t             \"-\" Causes the START_position of the gene region to be reduced by the user specified value." << endl
          << "              \t  FIXED mode: In the paramters file's \"FIXED mode\" section specify the +VALUE or -VALUE. It will be applied to all gene regions." << endl
          << "              \t  Uses a CUDA powered engine, therefore, requires a CUDA capable GPU." << endl
          << "              \t  File format is *.ehh (a tab deliminated text file)." << endl
          << endl
          << "EXTRAS:"
          << endl
          << endl
          << "-c or --cuda\t: Lists all available CUDA capable devices on machine." << endl
          << "            \t  Use the \"GPU number\" to select the desired CUDA device in the parameter file." << endl
          << endl
          << "-h or --help\t: Accesses this help menu where the software's currently available functions are listed." << endl;
}