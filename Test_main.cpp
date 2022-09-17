#include <iostream>
#include <filesystem>
#include <algorithm>
#include <unistd.h>
#include "test.cuh"
#include "cudaDevices.cuh"
#include "tajima.cuh"
#include "fu_li.cuh"
#include "fay_wu.cuh"
#include "neutral.cuh"
#include "ehh.cuh"
#include "mk_test.cuh"
#include "fst.cuh"
#include "gff2gene.cuh"
#include "hap_extract.cuh"
#include "vcf_splitter.h"
#include "print_param.h"
#include "gene_extract.h"
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

     /**
      * RENEGADES... LET US BEGIN
      * TODO: FINISH COMMENTING!
      **/

     print_Cate();
     cout << "CATE: CUDA Accelerated Testing of Evolution" << endl;
     cout << "Evolutionary tests for large scale genome data" << endl
          << "----------------------------------------------" << endl
          << endl;

     /**
      * The title of CATE is an acronym.
      * C.A.T.E stands for CUDA Accelerated Testing of Evolution.
      * It is a software designed with the aim of utilizing the computers multiprocessing technologies of both the CPU and GPU.
      * ! At present, CATE is compatible only with CUDA enabled NVIDIA GPU's.
      **/

     /**
      * Since CATE is command line executable it depends on the one or two parameters passed via the CLI.
      * @param argv is used to get these inputs.
      * * [1] will be used to get the function.
      * * [2] if required will point to the location of the parameter file.
      **/

     /**
      * @param gene_List can be specific to each function or be set as universal.
      * Universal enables all functions to access the same gene list file.
      **/

     /**
      * @param prometheus_Activate is used to activate the HPC catered high speed mode of CATE dubbed Prometheus.
      * ! Prometheus, is only available for the Neutrality tests of Tajima, Fay and Wu and, Fu and Li.
      **/

     if (argv[1] != NULL)
     {
          string function(argv[1]);
          /**
           * * Functions are converted to lowercase formats so that they will not be case sensetive.
           **/
          transform(function.begin(), function.end(), function.begin(), ::tolower);

          if (function == "--help" || function == "-h")
          {
               /**
                * Prints the help menu to the CLI.
                **/
               print_HELP();
          }
          else if (function == "--godmode")
          {
               /**
                * CATE testing framework not required nor usable by the end user.
                **/
               test try_something = test();
               try_something.thread_test();
          }
          else if (function == "--cuda" || function == "-c")
          {
               /**
                * Prints all available CUDA devices present on the current system.
                * User can use this list to determine which CUDA device to be used via the CUDA ID.
                **/
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

                    if (function == "--printparam" || function == "-pparam")
                    {
                         /**
                          * Prints a default parameter.json file for the user.
                          **/
                         print_param print = print_param(parameter_File);
                         print.ingress();

                         cout << endl
                              << "Parameter Print has been completed." << endl;
                    }
                    else
                    {
                         /**
                          * The parameter Class is used to read the parameter *.json file.
                          **/
                         parameter properties = parameter(parameter_File);

                         /**
                          * Configure, check for the presence and, create if auxillary folders are unavailable.
                          **/

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

                         /**
                          * The execution of all supplementary and main evolutionary functions.
                          **/

                         if (function == "--splitvcf" || function == "-svcf")
                         {
                              /**
                               * * Execution of the SPLIT VCF function.
                               * The function is responsible for the creation of the indexed folder hierarchy,
                               * which is crucial for the functioning of CATE.
                               *  * Enables the VCF files to be not only split by the number of SNPS, but by population and other
                               *  * pre-defined filters.
                               **/

                              // string input = properties.where("Input folder entry");
                              vcf_splitter split = vcf_splitter(tmp, properties.where("Input path"), properties.where("Population file path"), output_Path, properties.where_Int("Reference allele count"), properties.where_Int("Alternate allele count"), properties.where_Int("SNP count per file"), properties.where_Int("Sample_ID Column number"), properties.where_Int("Population_ID Column number"));
                              split.index_population();
                              split.read_File();
                              cout << endl
                                   << endl
                                   << "VCF split has been completed." << endl;
                         }
                         else if (function == "--splitfasta" || function == "-sfasta")
                         {
                              /**
                               * * Execution of the SPLIT FASTA function.
                               * Can split a merged FASTA file into its individual sequences.
                               * It can also be used to extract a specific sequence by its sequence ID.
                               **/

                              fasta_splitter split = fasta_splitter(properties.where("Raw FASTA file"), output_Path, properties.where("Sequence"));
                              split.ingress();

                              cout << endl
                                   << "FASTA split has been completed." << endl;
                         }
                         else if (function == "--mergefasta" || function == "-mfasta")
                         {
                              /**
                               * * Execution of the MERGE FASTA function.
                               * Merges all FASTA files in the folder into a single FASTA file.
                               **/

                              fasta_merge merge = fasta_merge(properties.where("FASTA files folder"), properties.where("Merge FASTA path"));
                              merge.ingress();

                              cout << endl
                                   << "FASTA merge has been completed." << endl;
                         }
                         else if (function == "--extractgenes" || function == "-egenes")
                         {
                              /**
                               * * Execution of the EXTRACT GENE function.
                               * Extracts gene sequences using a pre determined reference sequence file.
                               * The reference FASTA file should only have a single sequence.
                               * Positions of the gene list file should align with that of the reference sequence.
                               **/

                              string gene_List = properties.where("Extract gene list");
                              if (gene_List == "universal")
                              {
                                   gene_List = properties.where("Universal gene list");
                              }
                              gene_extract ge = gene_extract(gene_List, properties.where("Reference genome ex"), output_Path, intermediate_Path);
                              ge.ingress();

                              cout << endl
                                   << "Gene extractor has been completed." << endl;
                         }
                         else if (function == "--gff2gene" || function == "-g2g")
                         {
                              /**
                               * * Execution of the GFF to GENE function.
                               * Collects and writes all GENE's in a *.GFF file to a CATE format gene list file.
                               **/

                              // gff2gene(string input_File, string output_Path);
                              gff2gene g2g = gff2gene(properties.where("GFF file"), output_Path);
                              g2g.ingress();

                              cout << endl
                                   << "Completed GFF to Gene list" << endl;
                         }
                         else if (function == "--hapfromvcf" || function == "-hapext")
                         {
                              /**
                               * * Execution of the Haplotypes from VCF function.
                               * Uses the VCF file folder and the reference genome FASTA file to reconstruct all unique haplotypes in a gene region.
                               * Can also be used to recreate the FASTA file of the VCF file population. 
                               * ! Essentially ALL FASTA sequences per individual per ploidy will be generated.
                               **/

                              string gene_List = properties.where("Hap extract gene list");
                              if (gene_List == "universal")
                              {
                                   gene_List = properties.where("Universal gene list");
                              }

                              hap_extract haplotype_Extractor = hap_extract(gene_List, properties.where("Input path"), output_Path, properties.where_Int("CUDA Device ID"), intermediate_Path, properties.where_Int("Ploidy"), properties.where("Reference genome hap"), properties.where("Population out"));
                              haplotype_Extractor.ingress();

                              cout << endl
                                   << "Completed CUDA powered Haplotype extractor" << endl;
                         }
                         else if (function == "--tajima" || function == "-t")
                         {
                              string prometheus_Activate = properties.where("Prometheus activate");
                              transform(prometheus_Activate.begin(), prometheus_Activate.end(), prometheus_Activate.begin(), ::toupper);

                              string calc_Mode = properties.where("Calculation mode");
                              transform(calc_Mode.begin(), calc_Mode.end(), calc_Mode.begin(), ::toupper);

                              string gene_List;
                              if (calc_Mode == "FILE")
                              {
                                   gene_List = properties.where("Tajima gene list");
                                   if (gene_List == "universal")
                                   {
                                        gene_List = properties.where("Universal gene list");
                                   }
                              }

                              if (prometheus_Activate != "YES")
                              {
                                   if (calc_Mode == "FILE")
                                   {
                                        tajima tajimasD = tajima(gene_List, properties.where("Input path"), output_Path, properties.where_Int("CUDA Device ID"), intermediate_Path, properties.where_Int("Ploidy"));
                                        tajimasD.ingress();
                                   }
                                   else
                                   {
                                        tajima tajimasD_Window = tajima(calc_Mode, properties.where_Int("Window size"), properties.where_Int("Step size"), properties.where("Input path"), output_Path, properties.where_Int("CUDA Device ID"), properties.where_Int("Ploidy"));
                                        tajimasD_Window.ingress();
                                   }
                              }
                              else
                              {
                                   if (calc_Mode == "FILE")
                                   {
                                        tajima tajimasD_Prometheus = tajima(gene_List, properties.where("Input path"), output_Path, properties.where_Int("CUDA Device ID"), intermediate_Path, properties.where_Int("Ploidy"), prometheus_Activate, properties.where("Multi read"), properties.where_Int("Number of genes"), properties.where_Int("CPU cores"), properties.where_Int("SNPs per time"));
                                        tajimasD_Prometheus.ingress();
                                   }
                                   else
                                   {
                                        // WINDOW MODE
                                        tajima tajimasD_Prometheus_Window = tajima(calc_Mode, properties.where_Int("Window size"), properties.where_Int("Step size"), properties.where("Input path"), output_Path, properties.where_Int("CUDA Device ID"), properties.where_Int("Ploidy"), prometheus_Activate, properties.where("Multi read"), properties.where_Int("Number of genes"), properties.where_Int("CPU cores"), properties.where_Int("SNPs per time"));
                                        tajimasD_Prometheus_Window.ingress();
                                   }
                              }

                              cout << endl
                                   << "CUDA powered Tajima's D calculator has been completed." << endl;
                         }
                         else if (function == "--fuli" || function == "-f")
                         {
                              string prometheus_Activate = properties.where("Prometheus activate");
                              transform(prometheus_Activate.begin(), prometheus_Activate.end(), prometheus_Activate.begin(), ::toupper);

                              string calc_Mode = properties.where("Calculation mode");
                              transform(calc_Mode.begin(), calc_Mode.end(), calc_Mode.begin(), ::toupper);

                              string gene_List;
                              if (calc_Mode == "FILE")
                              {
                                   gene_List = properties.where("Fu and Li gene list");
                                   if (gene_List == "universal")
                                   {
                                        gene_List = properties.where("Universal gene list");
                                   }
                              }

                              if (prometheus_Activate != "YES")
                              {
                                   if (calc_Mode == "FILE")
                                   {
                                        fu_li fuli = fu_li(gene_List, properties.where("Input path"), output_Path, properties.where_Int("CUDA Device ID"), intermediate_Path, properties.where_Int("Ploidy"));
                                        fuli.ingress();
                                   }
                                   else
                                   {
                                        fu_li fuli_Window = fu_li(calc_Mode, properties.where_Int("Window size"), properties.where_Int("Step size"), properties.where("Input path"), output_Path, properties.where_Int("CUDA Device ID"), properties.where_Int("Ploidy"));
                                        fuli_Window.ingress();
                                   }
                              }
                              else
                              {
                                   if (calc_Mode == "FILE")
                                   {
                                        fu_li fuli_Prometheus = fu_li(gene_List, properties.where("Input path"), output_Path, properties.where_Int("CUDA Device ID"), intermediate_Path, properties.where_Int("Ploidy"), prometheus_Activate, properties.where("Multi read"), properties.where_Int("Number of genes"), properties.where_Int("CPU cores"), properties.where_Int("SNPs per time"));
                                        fuli_Prometheus.ingress();
                                   }
                                   else
                                   {
                                        // WINDOW MODE
                                        fu_li fuli_Prometheus_Window = fu_li(calc_Mode, properties.where_Int("Window size"), properties.where_Int("Step size"), properties.where("Input path"), output_Path, properties.where_Int("CUDA Device ID"), properties.where_Int("Ploidy"), prometheus_Activate, properties.where("Multi read"), properties.where_Int("Number of genes"), properties.where_Int("CPU cores"), properties.where_Int("SNPs per time"));
                                        fuli_Prometheus_Window.ingress();
                                   }
                              }

                              cout << "CUDA powered Fu and Li's D, D*, F and F* calculator has been completed." << endl;
                         }

                         else if (function == "--faywu" || function == "-w")
                         {
                              string prometheus_Activate = properties.where("Prometheus activate");
                              transform(prometheus_Activate.begin(), prometheus_Activate.end(), prometheus_Activate.begin(), ::toupper);

                              string calc_Mode = properties.where("Calculation mode");
                              transform(calc_Mode.begin(), calc_Mode.end(), calc_Mode.begin(), ::toupper);

                              string gene_List;
                              if (calc_Mode == "FILE")
                              {
                                   gene_List = properties.where("Fay and Wu gene list");
                                   if (gene_List == "universal")
                                   {
                                        gene_List = properties.where("Universal gene list");
                                   }
                              }

                              if (prometheus_Activate != "YES")
                              {
                                   if (calc_Mode == "FILE")
                                   {
                                        fay_wu faywu = fay_wu(gene_List, properties.where("Input path"), output_Path, properties.where_Int("CUDA Device ID"), intermediate_Path, properties.where_Int("Ploidy"));
                                        faywu.ingress();
                                   }
                                   else
                                   {
                                        // WINDOW MODE
                                        fay_wu faywu_Window = fay_wu(calc_Mode, properties.where_Int("Window size"), properties.where_Int("Step size"), properties.where("Input path"), output_Path, properties.where_Int("CUDA Device ID"), properties.where_Int("Ploidy"));
                                        faywu_Window.ingress();
                                   }
                              }
                              else
                              {
                                   if (calc_Mode == "FILE")
                                   {
                                        fay_wu faywu_Prometheus = fay_wu(gene_List, properties.where("Input path"), output_Path, properties.where_Int("CUDA Device ID"), intermediate_Path, properties.where_Int("Ploidy"), prometheus_Activate, properties.where("Multi read"), properties.where_Int("Number of genes"), properties.where_Int("CPU cores"), properties.where_Int("SNPs per time"));
                                        faywu_Prometheus.ingress();
                                   }
                                   else
                                   {
                                        // WINDOW MODE
                                        fay_wu faywu_Prometheus_Window = fay_wu(calc_Mode, properties.where_Int("Window size"), properties.where_Int("Step size"), properties.where("Input path"), output_Path, properties.where_Int("CUDA Device ID"), properties.where_Int("Ploidy"), prometheus_Activate, properties.where("Multi read"), properties.where_Int("Number of genes"), properties.where_Int("CPU cores"), properties.where_Int("SNPs per time"));
                                        faywu_Prometheus_Window.ingress();
                                   }
                              }

                              cout << "CUDA powered Fay and Wu's normalized H and E calculator has been completed." << endl;
                         }
                         else if (function == "--neutrality" || function == "-n")
                         {
                              string prometheus_Activate = properties.where("Prometheus activate");
                              transform(prometheus_Activate.begin(), prometheus_Activate.end(), prometheus_Activate.begin(), ::toupper);

                              string calc_Mode = properties.where("Calculation mode");
                              transform(calc_Mode.begin(), calc_Mode.end(), calc_Mode.begin(), ::toupper);

                              string gene_List;
                              if (calc_Mode == "FILE")
                              {
                                   gene_List = properties.where("Neutrality gene list");
                                   if (gene_List == "universal")
                                   {
                                        gene_List = properties.where("Universal gene list");
                                   }
                              }

                              if (prometheus_Activate != "YES")
                              {
                                   if (calc_Mode == "FILE")
                                   {
                                        neutral neutrality = neutral(gene_List, properties.where("Input path"), output_Path, properties.where_Int("CUDA Device ID"), intermediate_Path, properties.where_Int("Ploidy"));
                                        neutrality.ingress();
                                   }
                                   else
                                   {
                                        // WINDOW MODE
                                        neutral neutrality_Window = neutral(calc_Mode, properties.where_Int("Window size"), properties.where_Int("Step size"), properties.where("Input path"), output_Path, properties.where_Int("CUDA Device ID"), properties.where_Int("Ploidy"));
                                        neutrality_Window.ingress();
                                   }
                              }
                              else
                              {
                                   if (calc_Mode == "FILE")
                                   {
                                        neutral neutrality_Prometheus = neutral(gene_List, properties.where("Input path"), output_Path, properties.where_Int("CUDA Device ID"), intermediate_Path, properties.where_Int("Ploidy"), prometheus_Activate, properties.where("Multi read"), properties.where_Int("Number of genes"), properties.where_Int("CPU cores"), properties.where_Int("SNPs per time"));
                                        neutrality_Prometheus.ingress();
                                   }
                                   else
                                   {
                                        // WINDOW MODE
                                        neutral neutrality_Prometheus = neutral(calc_Mode, properties.where_Int("Window size"), properties.where_Int("Step size"), properties.where("Input path"), output_Path, properties.where_Int("CUDA Device ID"), properties.where_Int("Ploidy"), prometheus_Activate, properties.where("Multi read"), properties.where_Int("Number of genes"), properties.where_Int("CPU cores"), properties.where_Int("SNPs per time"));
                                        neutrality_Prometheus.ingress();
                                   }
                              }

                              cout << "CUDA powered complete neutrality test calculator has completed" << endl;
                         }
                         else if (function == "--mk" || function == "-m")
                         {
                              string gene_List = properties.where("McDonald–Kreitman gene list");
                              if (gene_List == "universal")
                              {
                                   gene_List = properties.where("Universal gene list");
                              }
                              mk_test mk = mk_test(properties.where("Reference genome mk"), properties.where("Alignment file"), gene_List, properties.where("Input path"), output_Path, properties.where_Int("CUDA Device ID"), intermediate_Path, properties.where_Int("Ploidy"), properties.where("Genetic code"), properties.where("Start codon(s)"), properties.where("Stop codon(s)"), properties.where("Alignment mode"), properties.where("ORF known"));
                              mk.ingress();

                              cout << "CUDA powered McDonald–Kreitman Neutrality Index (NI) test has been completed." << endl;
                         }
                         else if (function == "--fst" || function == "-x")
                         {
                              string calc_Mode = properties.where("Calculation mode");
                              transform(calc_Mode.begin(), calc_Mode.end(), calc_Mode.begin(), ::toupper);

                              if (calc_Mode == "FILE")
                              {
                                   string gene_List = properties.where("Fst gene list");
                                   if (gene_List == "universal")
                                   {
                                        gene_List = properties.where("Universal gene list");
                                   }

                                   fst fs = fst(gene_List, properties.where("Input path"), output_Path, properties.where_Int("CUDA Device ID"), intermediate_Path, properties.where_Int("Ploidy"), properties.where("Population index file path"), properties.where("Population ID"));
                                   fs.ingress();
                              }
                              else
                              {
                                   // fst(string calc_Mode, int window_Size, int step_Size, string gene_List, string input_Folder, string output_Path, int cuda_ID, int ploidy, string pop_Index_path, string pop_List);
                                   fst fst_Window = fst(calc_Mode, properties.where_Int("Window size"), properties.where_Int("Step size"), properties.where("Input path"), output_Path, properties.where_Int("CUDA Device ID"), properties.where_Int("Ploidy"), properties.where("Population index file path"), properties.where("Population ID"));
                                   fst_Window.ingress();
                              }

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
     cout << "   /\\\\\\\\\\\\\\\\\\     /\\\\\\\\\\\\\\\\\\     /\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\  /\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\n"
          << " /\\\\\\////////    /\\\\\\\\\\\\\\\\\\\\\\\\\\  \\///////\\\\\\/////  \\/\\\\\\///////////\n"
          << "/\\\\\\/            /\\\\\\/////////\\\\\\       \\/\\\\\\       \\/\\\\\\\n"
          << "/\\\\\\             \\/\\\\\\       \\/\\\\\\       \\/\\\\\\       \\/\\\\\\\\\\\\\\\\\\\\\\\n"
          << "\\/\\\\\\             \\/\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\       \\/\\\\\\       \\/\\\\\\///////\n"
          << " \\//\\\\\\            \\/\\\\\\/////////\\\\\\       \\/\\\\\\       \\/\\\\\\\n"
          << "   \\///\\\\\\          \\/\\\\\\       \\/\\\\\\       \\/\\\\\\       \\/\\\\\\\n"
          << "      \\////\\\\\\\\\\\\\\\\\\ \\/\\\\\\       \\/\\\\\\       \\/\\\\\\       \\/\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\n"
          << "          \\/////////  \\///        \\///        \\///        \\///////////////\n\n";

     // cout << "   0000000    00000   0000000   0000000" << endl;
     // cout << "  0000       00   00    000     00" << endl;
     // cout << " 000        00     00   000     00" << endl;
     // cout << "0           000000000   000     0000000" << endl;
     // cout << " 000        00     00   000     00" << endl;
     // cout << "  0000      00     00   000     00" << endl;
     // cout << "   0000000  00     00   000     0000000" << endl;
     // cout << endl;
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
          << "PROCESS MODES:"
          << endl
          << "A high performance and fully customizable mode called PROMETHEUS is available." << endl
          << endl
          << "PROMETHEUS is available for the three neutrality tests (Tajima's D, Fay and Wu tests and Fu and Li tests)." << endl
          << "PROMETHEUS is designed for power users on (High Performance Computing) HPC systems." << endl
          << "PROMETHEUS is activated via the parameters file. All other protocols of test execution remains the same." << endl
          << "PROMETHEUS uses a CUDA powered engine, therefore, requires a CUDA capable GPU." << endl
          << endl
          << "PROMETHEUS can be configured by the following five parameters:" << endl
          << "1. Prometheus activate: \"YES\" or \"NO\" paramters used to turn the mode ON or OFF." << endl
          << "2. CPU cores          : Controls the maximum number of cores that can be used at a time." << endl
          << "3. SNPs per time      : Controls the max number of SNPs that will be processed on the GPU at a time." << endl
          << "4. Number of genes    : Controls the number gene combinations that will be processed at a time." << endl
          << "5. Multi read         : \"YES\" or \"NO\" paramters used to control the ability to read multiple files at once." << endl
          << endl
          << "CALCULATION MODES:" << endl
          << "CATE can perform tests for pre-defined gene regions or the classical sliding window mechanism." << endl
          << endl
          << "Parameters for calculation mode is as follows:" << endl
          << "Calculation mode: Can be either \"WINDOW\" or \"FILE\"." << endl
          << endl
          << "If the calculation mode is \"WINDOW\" then the following two parameters need to be configured:" << endl
          << "1. Window size: Base pair size of the window or range of the combination." << endl
          << "2. Step size  : The base pair amount by which the next window's start will be incemented." << endl
          << "NOTE: If \"Step size\" is set to \"0\" then CATE will shift to a continous sliding window mode." << endl
          << endl
          << "If the calculation mode is \"FILE\" then the following two parameters need to be configured:" << endl
          << "1. Universal gene list: Configure the location for the tab deliminated gene list file for all tests." << endl
          << "2. * gene list        : Specify the location of the per test file or set it as \"universal\" to access the universal list." << endl
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
          << "-egenes or --extractgenes : Reads the gene list file to extract the gene sequences from the reference genome." << endl
          << "                            FASTA format reference genome must be specified." << endl
          << "                            All gene sequences will be generated into seperate FASTA files." << endl
          << endl
          << "--gff2gene or -g2g      : Creates the gene list file in a *.txt format from the input GFF3 file." << endl
          << "                          Note that only regions annotated as genes will be extracted." << endl
          << endl
          << "--hapfromvcf or -hapext : Extracts haplotypes and their sequences for a predefined gene list from a (split) VCF (indexed) folder provided the reference sequence." << endl
          << "                          The reference genome must be provided in a FASTA file." << endl
          << "                          The system will automatically identify each haplotype present." << endl
          << "                          In addition to the summary output each haplotype present for each gene will be generated in a seperate FASTA file." << endl
          << "                          IF \"Population out\" is set to \"YES\" then the entire population's FASTA configuration will be generated as well." << endl
          << "                          Uses a CUDA powered engine, therefore, requires a CUDA capable GPU." << endl
          << "                          File format is *.hsum (a tab deliminated text file)." << endl
          << endl
          << "-pparam or --printparam : Prints a sample layout of the parameter file to the specified location." << endl
          << "                          State the path with the name of the parameter file after the \"-pparam\" function." << endl
          << endl
          << "EVOLUTION TESTS:"
          << endl
          << "Core functions optimized for conducting Evolution tests on VCF files." << endl
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
          << "-n or --neutrality: Calculates the above three Neutrality tests (Tajima's Fu and Li's and Fay and Wu's) at once." << endl
          << "                    Uses a CUDA powered engine, therefore, requires a CUDA capable GPU." << endl
          << "                    File format is *.nt (a tab deliminated text file)." << endl
          << endl
          << "-m or --mk    \t: Calculates the McDonald–Kreitman Neutrality Index (NI) (1991) for a predefined gene list using a (split) vcf (indexed) folder." << endl
          << "              \t  The reference genome must be provided in a FASTA format file." << endl
          << "              \t  Two \"MODES\" exist. \"CHROM\" mode and \"GENE\" mode. Either must be specified" << endl
          << "              \t  CHROM mode: Conducts the test on an alignment does across the entire chromosomes" << endl
          << "              \t              Alignment file of the reference genome to the outgroup genome must also be provided in a *.maf format file." << endl
          << "              \t              PLEASE ensure that the REFERENCE sequence is first and OUTGROUP sequence is second in the MAF file." << endl
          << "              \t              ** TIP: Chromosome wide whole genome alignment software: GSAlign (https://github.com/hsinnan75/GSAlign)." << endl
          << "              \t  GENE mode : Conducts the test on alignments per gene between the reference gene and the outgroup gene." << endl
          << "              \t              Each gene's alignment file location must be provided as a third column in tab deliminated the gene list file." << endl
          << "              \t              Alignments must be provided in the blastn *.txt format." << endl
          << "              \t              PLEASE ensure that the REFERENCE gene sequence is the QUERY and OUTGROUP gene sequence is the SUBJECT." << endl
          << "              \t              NCBI's online blastn (https://blast.ncbi.nlm.nih.gov/Blast.cgi)," << endl
          << "              \t              or command line BLAST+ (https://anaconda.org/bioconda/blast) can be used." << endl
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