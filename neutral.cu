#include "neutral.cuh"
#include "functions.cuh"
#include "prometheus.cuh"

neutral::neutral(string gene_List, string input_Folder, string output_Path, int cuda_ID, string intermediate_Path, int ploidy)
{
     /**
      * * Constructor Function
      * NORMAL - GENE MODE constructor
      **/

     cout << "Initiating CUDA powered complete neutrality test calculator" << endl
          << "The following 3 tests will be calculated: " << endl
          << "1. Tajima's D" << endl
          << "2. Fu and Li's D, D*, F and F*" << endl
          << "3. Fay and Wu's normalized H and E" << endl
          << endl;

     set_Values(gene_List, input_Folder, output_Path, cuda_ID, intermediate_Path, ploidy);
}

neutral::neutral(string gene_List, string input_Folder, string output_Path, int cuda_ID, string intermediate_Path, int ploidy, string prometheus_Activate, string Multi_read, int number_of_genes, int CPU_cores, int SNPs_per_Run)
{
     /**
      * * Constructor Function
      * PROMETHEUS - GENE MODE constructor
      **/

     // Prometheus Neutral
     cout << "Initiating CUDA powered complete neutrality test calculator on PROMETHEUS" << endl
          << "The following 3 tests will be calculated: " << endl
          << "1. Tajima's D" << endl
          << "2. Fu and Li's D, D*, F and F*" << endl
          << "3. Fay and Wu's normalized H and E" << endl
          << endl;
     set_Values(gene_List, input_Folder, output_Path, cuda_ID, intermediate_Path, ploidy);

     this->prometheus_Activate = "YES";
     this->CPU_cores = CPU_cores;
     this->SNPs_per_Run = SNPs_per_Run;
     transform(Multi_read.begin(), Multi_read.end(), Multi_read.begin(), ::toupper);
     this->Multi_read = Multi_read;
     this->number_of_genes = number_of_genes;
}

neutral::neutral(string calc_Mode, int window_Size, int step_Size, string input_Folder, string output_Path, int cuda_ID, int ploidy, string prometheus_Activate, string Multi_read, int number_of_genes, int CPU_cores, int SNPs_per_Run)
{
     /**
      * * Constructor Function
      * PROMETHEUS - WINDOW MODE constructor
      **/

     // PROMETHEUS WINDOW MODE
     cout << "Initiating CUDA powered complete neutrality test calculator on PROMETHEUS" << endl
          << "The following 3 tests will be calculated: " << endl
          << "1. Tajima's D" << endl
          << "2. Fu and Li's D, D*, F and F*" << endl
          << "3. Fay and Wu's normalized H and E" << endl
          << endl;

     this->calc_Mode = "WINDOW";
     this->window_Size = window_Size;
     this->step_Size = step_Size;

     set_Values("", input_Folder, output_Path, cuda_ID, "", ploidy);

     this->prometheus_Activate = "YES";
     this->CPU_cores = CPU_cores;
     this->SNPs_per_Run = SNPs_per_Run;
     transform(Multi_read.begin(), Multi_read.end(), Multi_read.begin(), ::toupper);
     this->Multi_read = Multi_read;
     this->number_of_genes = number_of_genes;
}

neutral::neutral(string calc_Mode, int window_Size, int step_Size, string input_Folder, string output_Path, int cuda_ID, int ploidy)
{
     /**
      * * Constructor Function
      * NORMAL - WINDOW MODE constructor
      **/

     cout << "Initiating CUDA powered complete neutrality test calculator" << endl
          << "The following 3 tests will be calculated: " << endl
          << "1. Tajima's D" << endl
          << "2. Fu and Li's D, D*, F and F*" << endl
          << "3. Fay and Wu's normalized H and E" << endl
          << endl;

     this->calc_Mode = "WINDOW";
     this->window_Size = window_Size;
     this->step_Size = step_Size;

     set_Values("", input_Folder, output_Path, cuda_ID, "", ploidy);
}

void neutral::set_Values(string gene_List, string input_Folder, string output_Path, int cuda_ID, string intermediate_Path, int ploidy)
{
     /**
      * This function is used in conjunction with the constructor to set the common private variables.
      * Notifies the user if it is WINDOW mode or GENE (FILE) mode.
      * If WINDOW user is also notified if it is sliding window or normal step wise window mode.
      * Here the first call to the selected CUDA device occurs.
      **/

     if (this->calc_Mode != "WINDOW")
     {
          cout << "Calculation mode: FILE" << endl;
          this->gene_List = gene_List;
          cout << "Gene list file path\t: " << gene_List << endl;
     }
     else
     {
          cout << "Calculation mode: WINDOW" << endl;
          cout << "Window size: " << this->window_Size << endl;
          if (step_Size != 0)
          {
               cout << "Step size: " << this->step_Size << endl;
          }
          else
          {
               cout << "Sliding Window mode" << endl;
          }
          cout << endl;
     }
     cout << endl;
     this->input_Folder = input_Folder;
     this->output_Path = output_Path;
     this->intermediate_Path = intermediate_Path;
     this->ploidy = ploidy;

     cudaSetDevice(cuda_ID);
     cout << "Properties of selected CUDA GPU:" << endl;
     cudaDeviceProp prop;
     cudaGetDeviceProperties(&prop, cuda_ID);
     cout << "GPU number\t: " << cuda_ID << endl;
     cout << "GPU name\t: " << prop.name << endl;
     size_t l_free = 0;
     size_t l_Total = 0;
     cudaError_t error_id = cudaMemGetInfo(&l_free, &l_Total);
     cout << "GPU memory (GB)\t: " << l_Total / (1000 * 1000 * 1000) << endl;
     cout << "GPU number of multiprocessor(s)\t: " << prop.multiProcessorCount << endl;
     cout << "GPU block(s) per multiprocessor\t: " << prop.maxBlocksPerMultiProcessor << endl;
     this->tot_Blocks = prop.maxBlocksPerMultiProcessor;
     this->tot_ThreadsperBlock = prop.maxThreadsPerBlock;
     cout << "GPU thread(s) per block\t: " << tot_ThreadsperBlock << endl
          << endl;
}

void neutral::ingress()
{
     /**
      * Execution function.
      **/

     /**
      * Call the "functions" class. Bespoke functions commonly used by CATE.
      **/
     functions function = functions();

     /**
      * CATE indexed VCF folder is analyzed to extract the available super populations.
      * @param countries vector captures the available super populations.
      * Each population is processed separately.
      **/
     vector<string> countries = function.get_Countries(this->input_Folder);
     cout << countries.size() << " populations were found: ";
     for (int count = 0; count < countries.size(); count++)
     {
          string folder = countries[count];
          cout << folder.substr(folder.find_last_of("/") + 1, folder.length());
          if (count < countries.size() - 1)
          {
               cout << ", ";
          }
     }
     cout << endl
          << endl;

     for (string country : countries)
     {
          /**
           * To reiterate each population is processed separately.
           **/
          cout << "Processing country\t: " << country.substr(country.find_last_of("/") + 1, country.length()) << endl
               << endl;

          /**
           * @param folder_Index vector captures the sorted and indexed VCF file list from the query population folder.
           **/
          vector<pair<string, string>> folder_Index = function.index_Folder(country);
          cout << "Completed indexing folder\t: " << country << endl;

          cout << endl;
          /**
           * The first VCF file is read to obtain information of the sample size.
           * @param samples captures the sample size of the population under study.
           **/
          int samples = function.getN_Split(folder_Index[0].second);
          cout << "Number of samples in " << country.substr(country.find_last_of("/") + 1, country.length()) << " population\t: " << samples << endl;

          /**
           * @param N defines number of total sequences being present per SNP.
           **/
          int N = samples * ploidy;
          float N_float = (float)N;
          cout << "Number of sequences in " << country.substr(country.find_last_of("/") + 1, country.length()) << " population [ " << samples << " x " << ploidy << " ] (N)\t: " << N << endl;

          /**
           * @param combinations defines number of total pairwise combinations being present.
           **/
          long int combinations = function.combos_N(N);
          cout << "Pairwise combinations\t: " << combinations << endl;
          cout << endl;

          /**
           * Pre-requisite values needed for determination of all three neutrality tests' statistics.
           **/
          float an, e1, e2, vd, ud, vd_star, ud_star, uf, vf, uf_star, vf_star, bn, bn_plus1 = 0;
          get_Prerequisites(N, an, e1, e2, vd, ud, vd_star, ud_star, uf, vf, uf_star, vf_star, bn, bn_plus1);

          /**
           * @param test is used by Prometheus, to tell it which test is being processed.
           * T   = Tajima
           * FU  = Fu and Li
           * FA  = Fay and Wu
           * * N = All 3 Neutrality tests
           **/
          string test = "N";

          /**
           * Ensures which mode is being run. GENE (FILE) mode or WINDOW mode.
           **/
          if (this->calc_Mode != "FILE")
          {
               /**
                * * WINDOW mode configuration:
                **/

               /**
                * Output file is created for the population in the output folder for WINDOW mode.
                * @param output_File stores the output file's location.
                * The file name is a combination of the the country, window size and step size. Sliding window files will have a step size of 0.
                **/
               string output_File = output_Path + "/" +
                                    country.substr(country.find_last_of("/") + 1, country.length()) + "_" +
                                    to_string(window_Size) + "_" + to_string(step_Size) +
                                    ".nt";

               /**
                * Ensures if PROMETHEUS is being activated.
                **/
               if (prometheus_Activate == "YES")
               {
                    /**
                     * If Prometheus is being ACTIVATED then it is initialised accordingly.
                     **/
                    prometheus pro_Neutrality_Window = prometheus(output_File, window_Size, step_Size, folder_Index, Multi_read, tot_Blocks, tot_ThreadsperBlock, CPU_cores, SNPs_per_Run, number_of_genes, N, combinations, an, e1, e2, vd, ud, vd_star, ud_star, uf, vf, uf_star, vf_star, bn, bn_plus1);
                    /**
                     * Ensures if it is NORMAL window or SLIDING window mode.
                     * If step_Size is = 0 then it is sliding window mode.
                     **/
                    if (step_Size != 0)
                    {
                         pro_Neutrality_Window.process_Window(test);
                    }
                    else
                    {
                         /**
                          * Initiates processing of neutrality tests on PROMETHEUS on sliding window mode.
                          **/
                         pro_Neutrality_Window.process_C_sliding_Window(test);
                    }
               }
               else
               {
                    // PROMETHEUS OFF WINDOW MODE
                    /**
                     * If Prometheus is NOT being activated the window calls be done accordingly.
                     **/
                    if (step_Size != 0)
                    {
                         /**
                          * Initiates processing of Fu and Li on step wise window mode.
                          **/
                         window(output_File, an, e1, e2, vd, ud, vd_star, ud_star, uf, vf, uf_star, vf_star, bn, bn_plus1, N_float, combinations, folder_Index);
                    }
                    else
                    {
                         /**
                          * Initiates processing of Fu and Li on sliding window mode.
                          **/
                         window_Sliding(output_File, an, e1, e2, vd, ud, vd_star, ud_star, uf, vf, uf_star, vf_star, bn, bn_plus1, N_float, combinations, folder_Index);
                    }
               }
          }
          else
          {

               /**
                * * GENE (FILE) mode configuration:
                **/

               cout << "Prerequisites:" << endl
                    << endl;

               cout << "Common prerequisites:" << endl;
               cout << "an: " << an << "\t"
                    << "bn: " << bn << endl
                    << endl;

               cout << "Tajima's D prerequisites:" << endl;
               cout << "e1\t: " << e1 << "\t"
                    << "e2: " << e2 << endl
                    << endl;

               cout << "Fu and Li's prerequisites:" << endl;
               cout << "ud\t: " << ud << "\t"
                    << "vd\t: " << vd << endl;

               cout << "ud*\t: " << ud_star << "\t"
                    << "vd*\t: " << vd_star << endl;

               cout << endl;

               cout << "uf\t: " << uf << "\t"
                    << "vf\t: " << vf << endl;

               cout << "uf*\t: " << uf_star << "\t"
                    << "vf*\t: " << vf_star << endl
                    << endl;

               cout << "Fay and Wu Prerequisites:" << endl;
               cout << "bn_plus_1: " << bn_plus1 << endl
                    << endl;

               fstream gene_File;
               gene_File.open(gene_List, ios::in);
               cout << "Processing gene list:" << endl;

               /**
                * Output file is created for the population in the output folder for FILE mode.
                * @param output_File stores the output file's location.
                * The file name is a combination of the the country, and gene file name.
                **/
               string output_File = output_Path + "/" +
                                    country.substr(country.find_last_of("/") + 1, country.length()) + "_" +
                                    filesystem::path(gene_List).stem().string() +
                                    ".nt";

               /**
                * Log file created in the intermediate folder for the population.
                * @param intermediate_File stores the log file's location.
                * ! This helps with the resume function. Automatically resumes from the last completely processed gene in the event of a program crash.
                **/
               string intermediate_File = intermediate_Path + "/" +
                                          country.substr(country.find_last_of("/") + 1, country.length()) + "_" +
                                          filesystem::path(gene_List).stem().string() +
                                          ".log_nt";
               cout << endl;
               cout << "Writing to file\t: " << output_File << endl;
               cout << endl;

               /**
                * Initiate the reading of the gene file.
                **/
               if (gene_File.is_open())
               {
                    /**
                     * @param gene_Combo used to capture and extract info of each gene combination.
                     **/
                    string gene_Combo;

                    /**
                     * If the output file is absent this run will be considered as a brand new run of this query and,
                     * the output file and the intermediate log file will be created.
                     **/
                    if (filesystem::exists(output_File) == 0)
                    {
                         function.createFile(output_File, "Gene_name\tCoordinates\tPi\tS\tne\tns\tTotal_iEi\tTajimas_D\tD\tD_star\tF\tF_star\tFay_Wu_Normalized_H\tFay_Wu_Normalized_E");
                         function.createFile(intermediate_File);
                    }
                    else
                    {
                         /**
                          * If the intermediate log file present then the resume process will initiated.
                          * This is a unintelligent resume. Essentially it matches the each read line written with the lines read from the gene file.
                          * The break will occur as soon as their is a mismatch.
                          * To counter any errors it is advised to have a new gene file name or a new intermediate folder per new run.
                          **/
                         fstream intermediate;
                         intermediate.open(intermediate_File, ios::in);

                         /**
                          * @param get_finished comparison variable. Used o compare the intermediate file data with that of the gene file.
                          **/
                         string get_finished;
                         while (getline(intermediate, get_finished))
                         {
                              getline(gene_File, gene_Combo);
                              if (gene_Combo != get_finished)
                              {
                                   break;
                              }
                         }
                         intermediate.close();
                    }

                    fstream output;
                    fstream intermediate;
                    output.open(output_File, ios::app);
                    intermediate.open(intermediate_File, ios::app);

                    // PROMETHEUS HERE
                    /**
                     * Ensures if PROMETHEUS is being activated.
                     **/
                    if (prometheus_Activate == "YES")
                    {
                         cout << "Initializing Prometheus:" << endl
                              << endl;

                         /**
                          * If Prometheus is being ACTIVATED then it is initialised accordingly.
                          **/
                         prometheus pro_Neutrality = prometheus(folder_Index, Multi_read, tot_Blocks, tot_ThreadsperBlock, CPU_cores, SNPs_per_Run, number_of_genes, N, combinations, an, e1, e2, vd, ud, vd_star, ud_star, uf, vf, uf_star, vf_star, bn, bn_plus1);

                         /**
                          * @param gene_Collect vector is used to collect the batch of query regions to be processed by Prometheus at once.
                          **/
                         vector<string> gene_Collect;

                         while (getline(gene_File, gene_Combo))
                         {
                              gene_Collect.push_back(gene_Combo);
                              /**
                               * Ensures that the number of collected query regions match the user set limit to be processed at a time.
                               **/
                              if (gene_Collect.size() == number_of_genes)
                              {
                                   cout << "Prometheus batch intialized" << endl;
                                   cout << "From: " << gene_Collect[0] << endl;
                                   cout << "To  : " << gene_Collect[gene_Collect.size() - 1] << endl
                                        << endl;

                                   /**
                                    * LAUNCH Prometheus to process the collected query batch.
                                    * @param write_Lines vector collects the lines that should be written to the output file.
                                    */
                                   // launch prometheus
                                   vector<string> write_Lines = pro_Neutrality.collection_Engine(gene_Collect, test);
                                   // print
                                   cout << "System is writing Neutrality tests results" << endl;
                                   /**
                                    * Outputs are written and logs are made.
                                    **/
                                   for (size_t i = 0; i < write_Lines.size(); i++)
                                   {
                                        output << write_Lines[i] << "\n";
                                        intermediate << gene_Combo << "\n";
                                   }
                                   // clear prometheus
                                   output.flush();
                                   intermediate.flush();
                                   pro_Neutrality.erase();
                                   gene_Collect.clear();
                                   cout << endl;
                              }
                         }
                         if (gene_Collect.size() != 0)
                         {
                              /**
                               * If so then Prometheus is executed to process these regions.
                               **/
                              // RUN PROMETHEUS for remaining
                              // launch prometheus
                              cout << "Prometheus batch intialized" << endl;
                              cout << "From: " << gene_Collect[0] << endl;
                              cout << "To  : " << gene_Collect[gene_Collect.size() - 1] << endl
                                   << endl;

                              /**
                               * LAUNCH Prometheus to process the collected query batch.
                               * @param write_Lines vector collects the lines that should be written to the output file.
                               */
                              vector<string> write_Lines = pro_Neutrality.collection_Engine(gene_Collect, test);
                              // print
                              cout << "System is writing Neutrality tests results" << endl;
                              for (size_t i = 0; i < write_Lines.size(); i++)
                              {
                                   if (write_Lines[i] != "")
                                   {
                                        output << write_Lines[i] << "\n";
                                        intermediate << gene_Combo << "\n";
                                   }
                              }
                              cout << endl;
                         }

                         output.flush();
                         intermediate.flush();
                         pro_Neutrality.erase();
                         gene_Collect.clear();
                    }
                    else
                    {
                         /**
                          * If Prometheus is NOT activated each query gene region in the gene file is handled individually.
                          * This will be suitable for low powered systems and normal users.
                          * Because there will be no excessive use of resources nor any requirement to have extensive knowledge of your system.
                          **/
                         while (getline(gene_File, gene_Combo))
                         {
                              /**
                               * @param split_Data vector captures split function's outputs on the genes information.
                               **/
                              vector<string> split_Data;
                              function.split(split_Data, gene_Combo, '\t');

                              /**
                               * @param gene_Name captures the gene's name.
                               **/
                              string gene_Name = split_Data[0];
                              cout << "Gene name\t: " << gene_Name << endl;

                              /**
                               * @param coordinates vector captures split function's outputs on gene coordinates.
                               * [0] = chromosome
                               * [1] = start position
                               * [2] = end position
                               **/
                              vector<string> coordinates;
                              function.split(coordinates, split_Data[1], ':');

                              /**
                               * @param start_Co captures query gene's start position as an integer.
                               **/
                              /**
                               * @param end_Co captures query gene's end position as an integer.
                               **/
                              int start_Co = stoi(coordinates[1]);
                              int end_Co = stoi(coordinates[2]);
                              cout << "Coordinates\t: Chromosome: " << coordinates[0] << " Start: " << start_Co << " End: " << end_Co << endl;

                              /**
                               * @param tot_pairwise_Differences Fu and Li also requires the tot_pairwise_Differences in the query region to determine the average number of pairwise differences in the region.
                               * @param segregating_Sites  Fu and Li requires the total number of segregating sites/ SNPS in the query region.
                               * @param singletons_ns accounts for all singleton mutations (mutations present only once in the population) in the region.
                               * @param singletons_ne accounts for all singleton mutations that are different from the AA.
                               * @param theta_L  The mean number of mutations accumulated in each gene since the most recent common ancestor.
                               **/
                              float tot_pairwise_Differences = 0;
                              float theta_L = 0.00;
                              int segregating_Sites = 0;
                              int singletons_ns = 0;
                              int singletons_ne = 0;
                              int Total_iEi = 0;

                              vector<string> file_List;
                              cout << endl;
                              cout << "System is retrieving file(s)" << endl;
                              if (folder_Index.size() > 1)
                              {
                                   file_List = function.compound_interpolationSearch(folder_Index, start_Co, end_Co);
                              }
                              else
                              {
                                   file_List.push_back(folder_Index[0].second);
                              }
                              cout << "System has retrieved all file(s)" << endl;
                              cout << "System is collecting segregating site(s)" << endl;
                              vector<string> collect_Segregrating_sites;

                              for (string files : file_List)
                              {
                                   fstream file;
                                   file.open(files, ios::in);
                                   if (file.is_open())
                                   {
                                        string line;
                                        getline(file, line);
                                        while (getline(file, line))
                                        {
                                             vector<string> positions;
                                             function.split_getPos_ONLY(positions, line, '\t');
                                             int pos = stoi(positions[1]);

                                             if (pos >= start_Co && pos <= end_Co)
                                             {
                                                  collect_Segregrating_sites.push_back(line);
                                             }
                                             else if (pos > end_Co)
                                             {
                                                  break;
                                             }
                                        }
                                        file.close();
                                   }
                              }

                              // CUDA combined function
                              /**
                               * Calls the function to process the Segregating sites and obtain the necessary values to process all three neutrality tests.
                               */
                              process_Segs(collect_Segregrating_sites, N_float, segregating_Sites, tot_pairwise_Differences, singletons_ne, singletons_ns, Total_iEi, theta_L, tot_Blocks, tot_ThreadsperBlock);

                              cout << endl;
                              cout << "Total segregating sites (S)\t: " << segregating_Sites << endl;
                              cout << endl;

                              float pi = 0;

                              string Tajima_D, Fu_Li_D, Fu_Li_D_star, Fu_Li_F, Fu_Li_F_star, Fay_Wu_H, Fay_Wu_E;

                              if (segregating_Sites != 0)
                              {
                                   pi = (float)tot_pairwise_Differences / combinations;
                                   cout << "Average pairwise polymorphisms (pi)\t: " << pi << endl
                                        << endl;

                                   // float theta_squared, = 0;
                                   calculate_Neutrality(N_float, pi, segregating_Sites,
                                                        an, bn, e1, e2,
                                                        singletons_ne, singletons_ns, vd, ud, vd_star, ud_star, vf, uf, vf_star, uf_star,
                                                        theta_L, bn_plus1,
                                                        Tajima_D, Fu_Li_D, Fu_Li_D_star, Fu_Li_F, Fu_Li_F_star, Fay_Wu_H, Fay_Wu_E);

                                   cout << "Tajima's D\t: " << Tajima_D << endl
                                        << endl;

                                   // cout << "Total ns singletons\t: " << singletons_ns << endl;
                                   // cout << "Total ne singletons\t: " << singletons_ne << endl
                                   //      << endl;

                                   cout << "Fu and Li's D\t: " << Fu_Li_D << endl;
                                   cout << "Fu and Li's D*\t: " << Fu_Li_D_star << endl;
                                   cout << "Fu and Li's F\t: " << Fu_Li_F << endl;
                                   cout << "Fu and Li's F*\t: " << Fu_Li_F_star << endl
                                        << endl;

                                   // cout << "Theta_squared\t: " << theta_squared << endl;
                                   // cout << "Theta_L\t: " << theta_L << endl
                                   //      << endl;

                                   cout << "Fay and Wu's normalized H\t: " << Fay_Wu_H << endl;
                                   cout << "Fay and Wu's normalized E\t: " << Fay_Wu_E << endl;
                              }
                              else
                              {
                                   // cout << endl;
                                   cout << "Neutrality tests: Not Available" << endl;
                                   Tajima_D = "NA";
                                   Fu_Li_D = "NA";
                                   Fu_Li_D_star = "NA";
                                   Fu_Li_F = "NA";
                                   Fu_Li_F_star = "NA";
                                   Fay_Wu_H = "NA";
                                   Fay_Wu_E = "NA";
                              }

                              cout << endl;

                              //"Gene_name\tCoordinates\tPi\tS\tne\tns\tTotal_iEi\tTajimas_D\tD\tD_star\tF\tF_star\tFay_Wu_Normalized_H\tFay_Wu_Normalized_E"
                              output << gene_Name << "\t"
                                     << coordinates[0] << ":" << to_string(start_Co) << ":" << to_string(end_Co)

                                     << "\t" << to_string(pi)
                                     << "\t" << to_string(segregating_Sites)

                                     << "\t" << to_string(singletons_ne)
                                     << "\t" << to_string(singletons_ns)

                                     << "\t" << to_string(Total_iEi)

                                     << "\t" << Tajima_D

                                     << "\t" << Fu_Li_D
                                     << "\t" << Fu_Li_D_star
                                     << "\t" << Fu_Li_F
                                     << "\t" << Fu_Li_F_star

                                     << "\t" << Fay_Wu_H
                                     << "\t" << Fay_Wu_E
                                     << "\n";

                              intermediate << gene_Combo << "\n";
                              output.flush();
                              intermediate.flush();

                              // REMOVE AFTER TEST
                              // break;
                         }
                    }
                    output.close();
                    intermediate.close();
                    gene_File.close();
               }
               // REMOVE AFTER TEST
               // break;
          }
     }
}

void neutral::window_Sliding(string output_File, float an, float e1, float e2, float vd, float ud, float vd_star, float ud_star, float uf, float vf, float uf_star, float vf_star, float bn, float bn_plus1, float N_float, long int combinations, vector<pair<string, string>> &folder_Index)
{
     /**
      * NORMAL MODE SLIDING WINDOW FUNCTION
      **/

     /**
      * Call the "functions" class. Bespoke functions commonly used by CATE.
      **/
     functions function = functions();
     cout << "Writing to file\t: " << output_File << endl;
     cout << endl;

     /**
      * WINDOW functions have their own bespoke resume function that does not need an intermediate log file.
      * @param file_Count_Start is used to keep track of the files that have already been processed.
      * @param line_Num is used to keep track of the number of lines in that file that have already been processed.
      * ! This helps with the resume function. Automatically resumes from the last completely processed gene in the event of a program crash.
      **/
     int file_Count_Start = 0;
     int line_Num = 0;

     /**
      * If the output file is absent this run will be considered as a brand new run of this query and,
      * the output file and the intermediate log file will be created.
      **/
     if (filesystem::exists(output_File) == 0)
     {
          /**
           * Window outputs have NO gene name column.
           **/
          function.createFile(output_File, "Coordinates\tPi\tS\tne\tns\tTotal_iEi\tTajimas_D\tD\tD_star\tF\tF_star\tFay_Wu_Normalized_H\tFay_Wu_Normalized_E");
     }
     else
     {
          /**
           * If the output file is already present then the resume process will initiated.
           * This is a unintelligent resume. Essentially it matches the each read line written with the lines read from the gene file.
           * The break will occur as soon as their is a mismatch.
           * To counter any errors it is advised to have a new gene file name or a new intermediate folder per new run.
           * @param found acts as a boolean variable. found = 0 if the lines need to be skipped and will equal 1 when the resume position is found.
           **/
          int found = 0;

          /**
           * Open the output file to be begin finding the resume point.
           **/
          fstream output_Check;
          output_Check.open(output_File, ios::in);
          if (output_Check.is_open())
          {
               /**
                * @param line_Check is used to get the line from the output file to be compared.
                * First line is skipped cause it is a header line containing column names.
                **/
               string line_Check;
               getline(output_Check, line_Check); // skip first header line

               /**
                * We go through the files in the folder hierarchy one by one till we find the resume point.
                **/
               for (int file_Count = 0; file_Count < folder_Index.size(); file_Count++)
               {
                    /**
                     * @param file_Path gets the path of the query file being checked.
                     * @param line_Current gets the line number currently being checked.
                     **/
                    string file_Path = folder_Index[file_Count].second;
                    fstream file;
                    file.open(file_Path, ios::in);
                    int line_Current = 0;

                    if (file.is_open())
                    {
                         string line;
                         /**
                          * The first line of each VCF is skipped as it is the header line.
                          **/
                         getline(file, line); // skip first header line
                         while (getline(file, line))
                         {
                              line_Current++;
                              /**
                               * Checks if the line being queried is a valid seg site.
                               * If so it is checked if it has been already processed.
                               **/
                              int VALID = function.get_Valid(line);
                              if (VALID != -1)
                              {
                                   getline(output_Check, line_Check);
                                   string trim = line_Check.substr(0, line_Check.find('\t'));

                                   vector<string> positions;
                                   function.split_getPos_ONLY(positions, line, '\t');
                                   string pos = positions[1] + ":" + to_string((stoi(positions[1]) + window_Size));

                                   /**
                                    * Ensures the query line does not match that of the output
                                    **/
                                   if (pos != trim)
                                   {
                                        found = 1;
                                        file_Count_Start = file_Count;
                                        line_Num = line_Current;
                                        break;
                                   }
                              }
                         }
                         file.close();
                    }
                    /**
                     * If found is 1 that means the resume location has been found and the loop is broken.
                     **/
                    if (found == 1)
                    {
                         break;
                    }
               }
               output_Check.close();
          }
     }

     fstream output;
     output.open(output_File, ios::app);

     /**
      * @param line_Current is used to skip over the lines that have already been processed.
      **/
     int line_Current = 0;

     for (int file_Count = file_Count_Start; file_Count < folder_Index.size(); file_Count++)
     {
          /**
           * @param file_Path gets the path of the file being processed.
           **/
          string file_Path = folder_Index[file_Count].second;
          fstream file_Main;
          file_Main.open(file_Path, ios::in);

          if (file_Main.is_open())
          {
               string line_Main;
               /**
                * The first line of each VCF is skipped as it is the header line.
                **/
               getline(file_Main, line_Main); // skip first header line
               while (getline(file_Main, line_Main))
               {
                    /**
                     * Skips over lines that have already been processed.
                     **/
                    if (line_Current < line_Num)
                    {
                         line_Current++;
                    }
                    else
                    {
                         // check VALID
                         /**
                          * Checks if the line being queried is a valid seg site.
                          * If so it is processed.
                          * @param VALID captures the position of the query site if it is valid, else it returns -1.
                          **/
                         int VALID = function.get_Valid(line_Main);
                         // cout << line_Main << endl;
                         if (VALID != -1)
                         {
                              /**
                               * @param start_Co captures the start position as an integer.
                               * @param end_Co captures the end position as an integer.
                               **/
                              int start_Co = VALID;
                              int end_Co = start_Co + window_Size;

                              cout << "Coordinates\t: Start: " << start_Co << " End: " << end_Co << endl;

                              float tot_pairwise_Differences = 0;
                              float theta_L = 0.00;
                              int segregating_Sites = 0;
                              int singletons_ns = 0;
                              int singletons_ne = 0;
                              int Total_iEi = 0;

                              /**
                               * @param file_List vector is used to store the list of VCF files (found via CATES CIS algorithm) that satisfy the query region.
                               **/
                              vector<string> file_List;
                              cout << endl;
                              cout << "System is retrieving file(s)" << endl;
                              if (folder_Index.size() > 1)
                              {
                                   /**
                                    * IF only one file is present in the index folder that file will be used as is.
                                    **/
                                   file_List = function.compound_interpolationSearch(folder_Index, start_Co, end_Co);
                              }
                              else
                              {
                                   file_List.push_back(folder_Index[0].second);
                              }
                              cout << "System has retrieved all file(s)" << endl;
                              cout << "System is collecting segregating site(s)" << endl;
                              /**
                               * The SNPs (Segregating sites) that fall within the query region are collected from the VCF's.
                               * @param collect_Segregrating_sites vector stores the collected SNPs.
                               **/
                              vector<string> collect_Segregrating_sites;

                              /**
                               * Once the required files are found they are read sequentially to get the required SNP data for processing.
                               **/
                              for (string files : file_List)
                              {
                                   fstream file;
                                   file.open(files, ios::in);
                                   if (file.is_open())
                                   {
                                        string line;
                                        /**
                                         * The first line of each VCF is skipped as it is the header line.
                                         **/
                                        getline(file, line);
                                        while (getline(file, line))
                                        {
                                             vector<string> positions;
                                             /**
                                              * @param positions vector is used to capture the SNP data upto the position column (Column 2 (non zero count)).
                                              **/
                                             function.split_getPos_ONLY(positions, line, '\t');
                                             int pos = stoi(positions[1]);
                                             /**
                                              * Ensures that the query SNP's position satisfies the query region's range.
                                              **/
                                             if (pos >= start_Co && pos <= end_Co)
                                             {
                                                  /**
                                                   * If the SNP is between the range of the query region, it is collected.
                                                   * Information from the SNP is extracted via the GPU.
                                                   **/
                                                  collect_Segregrating_sites.push_back(line);
                                             }
                                             else if (pos > end_Co)
                                             {
                                                  /**
                                                   * If the read files query SNP exceeds the query regions range then the read loop is broken.
                                                   * This is because VCF's by nature, are sorted by position.
                                                   **/
                                                  break;
                                             }
                                        }
                                        file.close();
                                   }
                              }

                              // CUDA combined function
                              /**
                               * Calls the function to process the Segregating sites and obtain the necessary values to process all three neutrality tests.
                               */
                              process_Segs(collect_Segregrating_sites, N_float, segregating_Sites, tot_pairwise_Differences, singletons_ne, singletons_ns, Total_iEi, theta_L, tot_Blocks, tot_ThreadsperBlock);

                              cout << endl;
                              cout << "Total segregating sites (S)\t: " << segregating_Sites << endl;
                              cout << endl;

                              float pi = 0;

                              string Tajima_D, Fu_Li_D, Fu_Li_D_star, Fu_Li_F, Fu_Li_F_star, Fay_Wu_H, Fay_Wu_E;

                              if (segregating_Sites != 0)
                              {
                                   pi = (float)tot_pairwise_Differences / combinations;
                                   cout << "Average pairwise polymorphisms (pi)\t: " << pi << endl
                                        << endl;

                                   // float theta_squared, = 0;
                                   calculate_Neutrality(N_float, pi, segregating_Sites,
                                                        an, bn, e1, e2,
                                                        singletons_ne, singletons_ns, vd, ud, vd_star, ud_star, vf, uf, vf_star, uf_star,
                                                        theta_L, bn_plus1,
                                                        Tajima_D, Fu_Li_D, Fu_Li_D_star, Fu_Li_F, Fu_Li_F_star, Fay_Wu_H, Fay_Wu_E);

                                   cout << "Tajima's D\t: " << Tajima_D << endl
                                        << endl;

                                   // cout << "Total ns singletons\t: " << singletons_ns << endl;
                                   // cout << "Total ne singletons\t: " << singletons_ne << endl
                                   //      << endl;

                                   cout << "Fu and Li's D\t: " << Fu_Li_D << endl;
                                   cout << "Fu and Li's D*\t: " << Fu_Li_D_star << endl;
                                   cout << "Fu and Li's F\t: " << Fu_Li_F << endl;
                                   cout << "Fu and Li's F*\t: " << Fu_Li_F_star << endl
                                        << endl;

                                   // cout << "Theta_squared\t: " << theta_squared << endl;
                                   // cout << "Theta_L\t: " << theta_L << endl
                                   //      << endl;

                                   cout << "Fay and Wu's normalized H\t: " << Fay_Wu_H << endl;
                                   cout << "Fay and Wu's normalized E\t: " << Fay_Wu_E << endl;
                              }
                              else
                              {
                                   // cout << endl;
                                   cout << "Neutrality tests: Not Available" << endl;
                                   Tajima_D = "NA";
                                   Fu_Li_D = "NA";
                                   Fu_Li_D_star = "NA";
                                   Fu_Li_F = "NA";
                                   Fu_Li_F_star = "NA";
                                   Fay_Wu_H = "NA";
                                   Fay_Wu_E = "NA";
                              }

                              cout << endl;

                              //"Gene_name\tCoordinates\tPi\tS\tne\tns\tTotal_iEi\tTajimas_D\tD\tD_star\tF\tF_star\tFay_Wu_Normalized_H\tFay_Wu_Normalized_E"
                              output << to_string(start_Co) << ":" << to_string(end_Co)

                                     << "\t" << to_string(pi)
                                     << "\t" << to_string(segregating_Sites)

                                     << "\t" << to_string(singletons_ne)
                                     << "\t" << to_string(singletons_ns)

                                     << "\t" << to_string(Total_iEi)

                                     << "\t" << Tajima_D

                                     << "\t" << Fu_Li_D
                                     << "\t" << Fu_Li_D_star
                                     << "\t" << Fu_Li_F
                                     << "\t" << Fu_Li_F_star

                                     << "\t" << Fay_Wu_H
                                     << "\t" << Fay_Wu_E
                                     << "\n";

                              output.flush();
                         }
                    }
               }
               file_Main.close();
          }
     }
     output.close();
}

void neutral::window(string output_File, float an, float e1, float e2, float vd, float ud, float vd_star, float ud_star, float uf, float vf, float uf_star, float vf_star, float bn, float bn_plus1, float N_float, long int combinations, vector<pair<string, string>> &folder_Index)
{
     functions function = functions();

     int start_Value = stoi(folder_Index[0].first.substr(0, folder_Index[0].first.find('_')));
     int end_Value = stoi(folder_Index[folder_Index.size() - 1].first.substr(folder_Index[folder_Index.size() - 1].first.find('_') + 1));

     // cout << start_Value << endl;
     // cout << end_Value << endl;

     int start_Co = 0;
     int end_Co = start_Co + window_Size;

     while (start_Value > end_Co)
     {
          start_Co = start_Co + step_Size;
          end_Co = start_Co + window_Size;
     }

     cout << "Writing to file\t: " << output_File << endl;
     cout << endl;

     if (filesystem::exists(output_File) == 0)
     {
          function.createFile(output_File, "Coordinates\tPi\tS\tne\tns\tTotal_iEi\tTajimas_D\tD\tD_star\tF\tF_star\tFay_Wu_Normalized_H\tFay_Wu_Normalized_E");
     }
     else
     {
          // RESUME FUNCTION
          int caught = 0;

          // skipper
          fstream output;
          output.open(output_File, ios::in);

          string output_Line;

          while (start_Co <= end_Value)
          {

               // skip header
               getline(output, output_Line);

               while (getline(output, output_Line))
               {
                    string trim = output_Line.substr(0, output_Line.find('\t'));
                    string check = to_string(start_Co) + ":" + to_string(end_Co);
                    if (trim != check)
                    {
                         caught = 1;
                         break;
                    }
               }

               if (caught == 1)
               {
                    break;
               }

               start_Co = start_Co + step_Size;
               end_Co = start_Co + window_Size;
          }
          output.close();
     }

     fstream output;
     output.open(output_File, ios::app);

     while (start_Co <= end_Value)
     {
          cout << "Coordinates\t: Start: " << start_Co << " End: " << end_Co << endl;

          float tot_pairwise_Differences = 0;
          float theta_L = 0.00;
          int segregating_Sites = 0;
          int singletons_ns = 0;
          int singletons_ne = 0;
          int Total_iEi = 0;

          vector<string> file_List;
          cout << endl;
          cout << "System is retrieving file(s)" << endl;
          if (folder_Index.size() > 1)
          {
               file_List = function.compound_interpolationSearch(folder_Index, start_Co, end_Co);
          }
          else
          {
               file_List.push_back(folder_Index[0].second);
          }
          cout << "System has retrieved all file(s)" << endl;
          cout << "System is collecting segregrating site(s)" << endl;
          vector<string> collect_Segregrating_sites;

          for (string files : file_List)
          {
               fstream file;
               file.open(files, ios::in);
               if (file.is_open())
               {
                    string line;
                    getline(file, line);
                    while (getline(file, line))
                    {
                         vector<string> positions;
                         function.split_getPos_ONLY(positions, line, '\t');
                         int pos = stoi(positions[1]);

                         if (pos >= start_Co && pos <= end_Co)
                         {
                              collect_Segregrating_sites.push_back(line);
                         }
                         else if (pos > end_Co)
                         {
                              break;
                         }
                    }
                    file.close();
               }
          }

          // CUDA combined function
          process_Segs(collect_Segregrating_sites, N_float, segregating_Sites, tot_pairwise_Differences, singletons_ne, singletons_ns, Total_iEi, theta_L, tot_Blocks, tot_ThreadsperBlock);

          cout << endl;
          cout << "Total segregating sites (S)\t: " << segregating_Sites << endl;
          cout << endl;

          float pi = 0;

          string Tajima_D, Fu_Li_D, Fu_Li_D_star, Fu_Li_F, Fu_Li_F_star, Fay_Wu_H, Fay_Wu_E;

          if (segregating_Sites != 0)
          {
               pi = (float)tot_pairwise_Differences / combinations;
               cout << "Average pairwise polymorphisms (pi)\t: " << pi << endl
                    << endl;

               // float theta_squared, = 0;
               calculate_Neutrality(N_float, pi, segregating_Sites,
                                    an, bn, e1, e2,
                                    singletons_ne, singletons_ns, vd, ud, vd_star, ud_star, vf, uf, vf_star, uf_star,
                                    theta_L, bn_plus1,
                                    Tajima_D, Fu_Li_D, Fu_Li_D_star, Fu_Li_F, Fu_Li_F_star, Fay_Wu_H, Fay_Wu_E);

               cout << "Tajima's D\t: " << Tajima_D << endl
                    << endl;

               // cout << "Total ns singletons\t: " << singletons_ns << endl;
               // cout << "Total ne singletons\t: " << singletons_ne << endl
               //      << endl;

               cout << "Fu and Li's D\t: " << Fu_Li_D << endl;
               cout << "Fu and Li's D*\t: " << Fu_Li_D_star << endl;
               cout << "Fu and Li's F\t: " << Fu_Li_F << endl;
               cout << "Fu and Li's F*\t: " << Fu_Li_F_star << endl
                    << endl;

               // cout << "Theta_squared\t: " << theta_squared << endl;
               // cout << "Theta_L\t: " << theta_L << endl
               //      << endl;

               cout << "Fay and Wu's normalized H\t: " << Fay_Wu_H << endl;
               cout << "Fay and Wu's normalized E\t: " << Fay_Wu_E << endl;
          }
          else
          {
               // cout << endl;
               cout << "Neutrality tests: Not Available" << endl;
               Tajima_D = "NA";
               Fu_Li_D = "NA";
               Fu_Li_D_star = "NA";
               Fu_Li_F = "NA";
               Fu_Li_F_star = "NA";
               Fay_Wu_H = "NA";
               Fay_Wu_E = "NA";
          }

          cout << endl;

          //"Gene_name\tCoordinates\tPi\tS\tne\tns\tTotal_iEi\tTajimas_D\tD\tD_star\tF\tF_star\tFay_Wu_Normalized_H\tFay_Wu_Normalized_E"
          output << to_string(start_Co) << ":" << to_string(end_Co)

                 << "\t" << to_string(pi)
                 << "\t" << to_string(segregating_Sites)

                 << "\t" << to_string(singletons_ne)
                 << "\t" << to_string(singletons_ns)

                 << "\t" << to_string(Total_iEi)

                 << "\t" << Tajima_D

                 << "\t" << Fu_Li_D
                 << "\t" << Fu_Li_D_star
                 << "\t" << Fu_Li_F
                 << "\t" << Fu_Li_F_star

                 << "\t" << Fay_Wu_H
                 << "\t" << Fay_Wu_E
                 << "\n";

          output.flush();

          start_Co = start_Co + step_Size;
          end_Co = start_Co + window_Size;
     }

     output.close();
}

void neutral::calculate_Neutrality(float N_float, float &pi, int &segregating_Sites, float &an, float &bn, float &e1, float &e2, int &singletons_ne, int &singletons_ns, float &vd, float &ud, float &vd_star, float &ud_star, float &vf, float &uf, float &vf_star, float &uf_star, float &theta_L, float &bn_plus1, string &Tajima_D, string &Fu_Li_D, string &Fu_Li_D_star, string &Fu_Li_F, string &Fu_Li_F_star, string &Fay_Wu_H, string &Fay_Wu_E)
{
     promise<string> promise_Tajima;
     future<string> future_Tajima = promise_Tajima.get_future();

     promise<string> promise_D;
     future<string> future_D = promise_D.get_future();
     promise<string> promise_D_star;
     future<string> future_D_star = promise_D_star.get_future();
     promise<string> promise_F;
     future<string> future_F = promise_F.get_future();
     promise<string> promise_F_star;
     future<string> future_F_star = promise_F_star.get_future();

     promise<string> promise_H;
     future<string> future_H = promise_H.get_future();
     promise<string> promise_E;
     future<string> future_E = promise_E.get_future();

     // Tajimas_D_thread(promise<string> &Tajimas_D_value, float pi, int segregating_Sites, float an, float e1, float e2);
     thread tajima_thread{&neutral::Tajimas_D_thread, this, ref(promise_Tajima), pi, segregating_Sites, an, e1, e2};
     // Fu_li_thread(promise<string> Fu_Li_D, promise<string> Fu_Li_D_star, promise<string> Fu_Li_F, promise<string> Fu_Li_F_star, float N_float, float pi, int segregating_Sites, float an, int singletons_ne, int singletons_ns, float vd, float ud, float vd_star, float ud_star, float vf, float uf, float vf_star, float uf_star)
     thread fu_li_thread{&neutral::Fu_li_thread, this, ref(promise_D), ref(promise_D_star), ref(promise_F), ref(promise_F_star), N_float, pi, segregating_Sites, an, singletons_ne, singletons_ns, vd, ud, vd_star, ud_star, vf, uf, vf_star, uf_star};
     // Fay_wu_thread(promise<string> Fay_Wu_H, promise<string> Fay_Wu_E, float N_float, float pi, int segregating_Sites, float an, float bn, float bn_plus1, float theta_L)
     thread fay_wu_thread{&neutral::Fay_wu_thread, this, ref(promise_H), ref(promise_E), N_float, pi, segregating_Sites, an, bn, bn_plus1, theta_L};

     Tajima_D = future_Tajima.get();

     Fu_Li_D = future_D.get();
     Fu_Li_D_star = future_D_star.get();
     Fu_Li_F = future_F.get();
     Fu_Li_F_star = future_F_star.get();

     Fay_Wu_H = future_H.get();
     Fay_Wu_E = future_E.get();

     tajima_thread.join();
     fu_li_thread.join();
     fay_wu_thread.join();
}

void neutral::Fay_wu_thread(promise<string> &Fay_Wu_H, promise<string> &Fay_Wu_E, float N_float, float pi, int segregating_Sites, float an, float bn, float bn_plus1, float theta_L)
{
     float S = (float)segregating_Sites;
     float theta_squared = (float)(S * (S - 1)) / (pow(an, 2) + bn);
     float theta_W = (float)S / an;

     float VAR_pi_MINUS_theta_L = (float)(((N_float - 2.0) / (6.0 * (N_float - 1.0))) * theta_W) + ((((18.0 * pow(N_float, 2) * ((3.0 * N_float) + 2.0) * bn_plus1) - ((88.0 * pow(N_float, 3)) + (9.0 * pow(N_float, 2)) - (13.0 * N_float) + 6.0)) / (9.0 * N_float * pow(N_float - 1, 2))) * theta_squared);
     float VAR_theta_L_MINUS_theta_W = (float)(((N_float / (2.0 * (N_float - 1.0))) - (1.0 / an)) * theta_W) + (((bn / (pow(an, 2))) + (2.0 * pow((N_float / (N_float - 1.0)), 2) * bn) - ((2.0 * ((N_float * bn) - N_float + 1.0)) / ((N_float - 1.0) * an)) - (((3.0 * N_float) + 1) / (N_float - 1.0))) * theta_squared);

     float H = ((float)(pi - theta_L)) / (sqrt(VAR_pi_MINUS_theta_L));
     float E = ((float)(theta_L - theta_W)) / (sqrt(VAR_theta_L_MINUS_theta_W));

     Fay_Wu_H.set_value(to_string(H));
     Fay_Wu_E.set_value(to_string(E));
}

void neutral::Tajimas_D_thread(promise<string> &Tajimas_D_value, float pi, int segregating_Sites, float an, float e1, float e2)
{
     float D_Tajima = (float)(pi - (segregating_Sites / an)) / sqrt(((e1 * segregating_Sites) + (e2 * segregating_Sites * (segregating_Sites - 1))));

     Tajimas_D_value.set_value(to_string(D_Tajima));
}

void neutral::Fu_li_thread(promise<string> &Fu_Li_D, promise<string> &Fu_Li_D_star, promise<string> &Fu_Li_F, promise<string> &Fu_Li_F_star, float N_float, float pi, int segregating_Sites, float an, int singletons_ne, int singletons_ns, float vd, float ud, float vd_star, float ud_star, float vf, float uf, float vf_star, float uf_star)
{
     float D = (float)(segregating_Sites - (an * singletons_ne)) / sqrt(((ud * segregating_Sites) + (vd * (pow(segregating_Sites, 2)))));
     float D_star = (float)(((N_float / (N_float - 1)) * segregating_Sites) - (an * singletons_ns)) / sqrt(((ud_star * segregating_Sites) + (vd_star * (pow(segregating_Sites, 2.0)))));
     float F = (float)(pi - singletons_ne) / sqrt(((uf * segregating_Sites) + (vf * (pow(segregating_Sites, 2)))));
     float F_star = (float)(pi - (((N_float - 1) / N_float) * singletons_ns)) / sqrt(((uf_star * segregating_Sites) + (vf_star * (pow(segregating_Sites, 2.0)))));

     Fu_Li_D.set_value(to_string(D));
     Fu_Li_D_star.set_value(to_string(D_star));
     Fu_Li_F.set_value(to_string(F));
     Fu_Li_F_star.set_value(to_string(F_star));
}

__global__ void cuda_process_Segs(char *sites, int *index, int num_Segregrating_sites, int *theta_Partials, int *VALID_or_NOT, int *MA_count, int *ne, int *ns)
{
     int tid = threadIdx.x + blockIdx.x * blockDim.x;

     while (tid < num_Segregrating_sites)
     {
          int column = 0;
          int site_Start = index[tid];
          int site_End = index[tid + 1];

          int i = site_Start;

          char REF = 'N';
          while (column < 3)
          {
               if (sites[i] == '\t')
               {
                    column++;
               }
               i++;
          }

          if (sites[i] >= 97)
          {
               REF = sites[i] - 32;
          }
          else
          {
               REF = sites[i];
          }

          char ALT = 'N';
          while (column < 4)
          {
               if (sites[i] == '\t')
               {
                    column++;
               }
               i++;
          }

          if (sites[i] >= 97)
          {
               ALT = sites[i] - 32;
          }
          else
          {
               ALT = sites[i];
          }

          // printf("REF: %c", REF);
          // printf("\n");
          // printf("ALT: %c", ALT);
          // printf("\n");

          while (column < 7)
          {
               if (sites[i] == '\t')
               {
                    column++;
               }
               i++;
          }

          char AA = 'N';

          while (column < 8)
          {
               if ((sites[i] == 'A' && sites[i + 1] == 'A' && sites[i + 2] == '=') || (sites[i] == 'a' && sites[i + 1] == 'a' && sites[i + 2] == '='))
               {
                    char CHECK = 'N';

                    if (sites[i + 3] >= 97)
                    {
                         CHECK = sites[i + 3] - 32;
                    }
                    else
                    {
                         CHECK = sites[i + 3];
                    }

                    if (CHECK == 'A' || CHECK == 'T' || CHECK == 'G' || CHECK == 'C')
                    {
                         AA = CHECK;
                    }
               }

               if (sites[i] == '\t')
               {
                    column++;
               }
               i++;
          }

          if (AA == 'N')
          {
               AA = REF;
          }

          while (column < 9)
          {
               if (sites[i] == '\t')
               {
                    column++;
               }
               i++;
          }

          // printf("Column 1: %c\n", sites[i]);

          int ALT_count = 0;
          int REF_count = 0;
          while (i < site_End)
          {
               if (sites[i] == '1')
               {
                    ALT_count = ALT_count + 1;
               }
               else if (sites[i] == '0')
               {
                    REF_count = REF_count + 1;
               }

               i++;
          }

          if (ALT_count == 0 || REF_count == 0)
          {
               VALID_or_NOT[tid] = 0;
               MA_count[tid] = 0;
               theta_Partials[tid] = 0;
          }
          else
          {
               VALID_or_NOT[tid] = 1;
               char MA = 'N';

               if (ALT_count < REF_count)
               {
                    MA = ALT;
                    MA_count[tid] = ALT_count;
               }
               else
               {
                    MA = REF;
                    MA_count[tid] = REF_count;
               }

               if (MA_count[tid] == 1)
               {
                    ns[tid] = 1;

                    if (MA == AA)
                    {
                         ne[tid] = 0;
                    }
                    else
                    {
                         ne[tid] = 1;
                    }
               }
               else
               {
                    ne[tid] = 0;
                    ns[tid] = 0;
               }

               int theta_Partial = 0;
               if (MA != AA)
               {
                    theta_Partial = MA_count[tid];
               }

               theta_Partials[tid] = theta_Partial;
          }

          tid += blockDim.x * gridDim.x;
     }
}

void neutral::process_Segs(vector<string> &total_Segregrating_sites, float N, int &real_segregrating_Sites, float &tot_pairwise_Differences, int &singletons_ne, int &singletons_ns, int &Total_iEi, float &theta_L, int tot_Blocks, int tot_ThreadsperBlock)
{
     cout << "System is processing and filtering segregrating site(s)" << endl;

     int num_segregrating_Sites = total_Segregrating_sites.size();
     theta_L = 0.00;
     string Seg_sites = "";
     int site_Index[num_segregrating_Sites + 1];
     site_Index[0] = 0;

     for (size_t i = 0; i < num_segregrating_Sites; i++)
     {
          Seg_sites.append(total_Segregrating_sites[i]);
          site_Index[i + 1] = site_Index[i] + total_Segregrating_sites[i].size();
     }

     char *full_Char;
     full_Char = (char *)malloc((Seg_sites.size() + 1) * sizeof(char));
     strcpy(full_Char, Seg_sites.c_str());

     total_Segregrating_sites.clear();

     char *cuda_full_Char;
     cudaMallocManaged(&cuda_full_Char, (Seg_sites.size() + 1) * sizeof(char));
     int *cuda_site_Index;
     cudaMallocManaged(&cuda_site_Index, (num_segregrating_Sites + 1) * sizeof(int));

     int *cuda_VALID_or_NOT, *VALID_or_NOT;
     cudaMallocManaged(&cuda_VALID_or_NOT, num_segregrating_Sites * sizeof(int));
     VALID_or_NOT = (int *)malloc(num_segregrating_Sites * sizeof(int));

     int *cuda_MA_Count, *MA_Count;
     cudaMallocManaged(&cuda_MA_Count, num_segregrating_Sites * sizeof(int));
     MA_Count = (int *)malloc(num_segregrating_Sites * sizeof(int));

     int *cuda_Theta_partials, *Theta_partials;
     cudaMallocManaged(&cuda_Theta_partials, num_segregrating_Sites * sizeof(int));
     Theta_partials = (int *)malloc(num_segregrating_Sites * sizeof(int));

     int *ne_CUDA, *ne, *ns_CUDA, *ns;
     cudaMallocManaged(&ne_CUDA, num_segregrating_Sites * sizeof(int));
     cudaMallocManaged(&ns_CUDA, num_segregrating_Sites * sizeof(int));
     ne = (int *)malloc(num_segregrating_Sites * sizeof(int));
     ns = (int *)malloc(num_segregrating_Sites * sizeof(int));

     cudaMemcpy(cuda_full_Char, full_Char, (Seg_sites.size() + 1) * sizeof(char), cudaMemcpyHostToDevice);
     cudaMemcpy(cuda_site_Index, site_Index, (num_segregrating_Sites + 1) * sizeof(int), cudaMemcpyHostToDevice);

     // cuda_process_Segs(char *sites, int *index, int num_Segregrating_sites, int *theta_Partials, int *VALID_or_NOT, int *MA_count, int *ne, int *ns)
     cuda_process_Segs<<<tot_Blocks, tot_ThreadsperBlock>>>(cuda_full_Char, cuda_site_Index, num_segregrating_Sites, cuda_Theta_partials, cuda_VALID_or_NOT, cuda_MA_Count, ne_CUDA, ns_CUDA);
     cudaDeviceSynchronize();

     cudaMemcpy(VALID_or_NOT, cuda_VALID_or_NOT, num_segregrating_Sites * sizeof(int), cudaMemcpyDeviceToHost);
     cudaMemcpy(MA_Count, cuda_MA_Count, num_segregrating_Sites * sizeof(int), cudaMemcpyDeviceToHost);
     cudaMemcpy(ne, ne_CUDA, num_segregrating_Sites * sizeof(int), cudaMemcpyDeviceToHost);
     cudaMemcpy(ns, ns_CUDA, num_segregrating_Sites * sizeof(int), cudaMemcpyDeviceToHost);
     cudaMemcpy(Theta_partials, cuda_Theta_partials, num_segregrating_Sites * sizeof(int), cudaMemcpyDeviceToHost);

     free(full_Char);
     cudaFree(cuda_full_Char);
     cudaFree(cuda_site_Index);
     cudaFree(cuda_MA_Count);
     cudaFree(cuda_VALID_or_NOT);
     cudaFree(cuda_Theta_partials);
     cudaFree(ns_CUDA);
     cudaFree(ne_CUDA);

     real_segregrating_Sites = 0;
     tot_pairwise_Differences = 0;
     singletons_ne = 0;
     singletons_ns = 0;
     int total_iTheta = 0;

     for (size_t i = 0; i < num_segregrating_Sites; i++)
     {
          if (VALID_or_NOT[i] == 1)
          {
               real_segregrating_Sites = real_segregrating_Sites + 1;
               float MAF = (float)MA_Count[i] / N;
               tot_pairwise_Differences = tot_pairwise_Differences + (MAF * (1 - MAF) * pow(N, 2));
               total_iTheta = total_iTheta + Theta_partials[i];
               singletons_ne = singletons_ne + ne[i];
               singletons_ns = singletons_ns + ns[i];
          }
     }

     Total_iEi = total_iTheta;
     theta_L = (float)total_iTheta / (float)(N - 1);

     free(Theta_partials);
     free(MA_Count);
     free(VALID_or_NOT);
     free(ns);
     free(ne);

     cout << "System has completed processing the segregrating site(s)" << endl;
}

__global__ void cuda_pre_Calculation(int N, float *a1_CUDA, float *a2_CUDA)
{
     int tid = threadIdx.x + blockIdx.x * blockDim.x;

     while (tid < N)
     {
          a1_CUDA[tid] = (float)1 / (tid + 1);
          a2_CUDA[tid] = (float)1 / ((tid + 1) * (tid + 1));
          tid += blockDim.x * gridDim.x;
     }
}

void neutral::get_Prerequisites(int N_tot, float &an, float &e1, float &e2, float &vd, float &ud, float &vd_star, float &ud_star, float &uf, float &vf, float &uf_star, float &vf_star, float &bn, float &bn_plus1)
{
     int N = N_tot - 1;
     float *a1_CUDA, *a2_CUDA;
     float *a1_partial, *a2_partial;

     a1_partial = (float *)malloc(N * sizeof(float));
     a2_partial = (float *)malloc(N * sizeof(float));

     cudaMallocManaged(&a1_CUDA, N * sizeof(int));
     cudaMallocManaged(&a2_CUDA, N * sizeof(int));

     cuda_pre_Calculation<<<tot_Blocks, tot_ThreadsperBlock>>>(N, a1_CUDA, a2_CUDA);
     cudaDeviceSynchronize();

     cudaMemcpy(a1_partial, a1_CUDA, N * sizeof(float), cudaMemcpyDeviceToHost);
     cudaMemcpy(a2_partial, a2_CUDA, N * sizeof(float), cudaMemcpyDeviceToHost);

     cudaFree(a1_CUDA);
     cudaFree(a2_CUDA);

     an = 0;
     bn = 0;
     for (size_t i = 0; i < N; i++)
     {
          an += a1_partial[i];
          bn += a2_partial[i];
     }

     free(a1_partial);
     free(a2_partial);

     float N_float_tot = (float)N_tot;
     float N_float = (float)N;

     promise<float> promise_e1;
     promise<float> promise_e2;
     future<float> future_e1 = promise_e1.get_future();
     future<float> future_e2 = promise_e2.get_future();

     promise<float> promise_vd;
     promise<float> promise_ud;
     promise<float> promise_vd_star;
     promise<float> promise_ud_star;
     promise<float> promise_uf;
     promise<float> promise_vf;
     promise<float> promise_uf_star;
     promise<float> promise_vf_star;

     future<float> future_vd = promise_vd.get_future();
     future<float> future_ud = promise_ud.get_future();
     future<float> future_vd_star = promise_vd_star.get_future();
     future<float> future_ud_star = promise_ud_star.get_future();
     future<float> future_uf = promise_uf.get_future();
     future<float> future_vf = promise_vf.get_future();
     future<float> future_uf_star = promise_uf_star.get_future();
     future<float> future_vf_star = promise_vf_star.get_future();

     thread tajima_thread{&neutral::pre_Tajima, this, N_float_tot, an, bn, ref(promise_e1), ref(promise_e2)};
     thread fu_li_thread{&neutral::pre_Fu_li, this, N_float_tot, N_float, an, bn, ref(promise_vd), ref(promise_ud), ref(promise_vd_star), ref(promise_ud_star), ref(promise_uf), ref(promise_vf), ref(promise_uf_star), ref(promise_vf_star)};

     e1 = future_e1.get();
     e2 = future_e2.get();

     vd = future_vd.get();
     ud = future_ud.get();
     vd_star = future_vd_star.get();
     ud_star = future_ud_star.get();
     uf = future_uf.get();
     vf = future_vf.get();
     uf_star = future_uf_star.get();
     vf_star = future_vf_star.get();

     tajima_thread.join();
     fu_li_thread.join();

     bn_plus1 = bn + (1.0 / pow(N_float_tot, 2));
}

void neutral::pre_Tajima(float N_tot, float an, float bn, promise<float> &e1_set, promise<float> &e2_set)
{
     float b1 = (float)(N_tot + 1) / (3 * (N_tot - 1));
     float b2 = (float)(2 * (pow(N_tot, 2.0) + N_tot + 3)) / ((9 * N_tot) * (N_tot - 1));

     float c1 = b1 - (1 / an);
     float c2 = b2 - ((N_tot + 2) / (an * N_tot)) + (bn / pow(an, 2.0));

     float e1 = c1 / an;
     float e2 = c2 / (pow(an, 2.0) + bn);

     e1_set.set_value(e1);
     e2_set.set_value(e2);
}

void neutral::pre_Fu_li(float N_float_tot, float N_float, float an, float bn, promise<float> &vd_set, promise<float> &ud_set, promise<float> &vd_star_set, promise<float> &ud_star_set, promise<float> &uf_set, promise<float> &vf_set, promise<float> &uf_star_set, promise<float> &vf_star_set)
{
     // float N_float = N_float_tot - 1;

     float an_plus_1 = (float)an + (1 / N_float_tot);

     float cn = (float)2 * (((N_float_tot * an) - (2 * N_float)) / (N_float * (N_float_tot - 2)));
     float an_square = (float)pow(an, 2.0);

     float vd = (float)1 + ((an_square / (bn + an_square)) * (cn - (N_float_tot + 1) / N_float));
     float ud = (float)an - 1 - vd;

     float dn = (float)cn + ((N_float_tot - 2) / pow(N_float, 2.0)) + ((2 / N_float) * (1.5 - (((2 * an_plus_1) - 3) / (N_float_tot - 2)) - (1 / N_float_tot)));

     float vd_star = (float)((pow((N_float_tot / N_float), 2.0) * bn) + (pow(an, 2.0) * dn) - (2 * ((N_float_tot * an * (an + 1)) / (pow(N_float, 2.0))))) / (pow(an, 2.0) + bn);
     float ud_star = (float)N_float_tot / N_float * (an - (N_float_tot / N_float)) - vd_star;

     float vf = (float)(cn + ((2 * (pow(N_float_tot, 2) + N_float_tot + 3)) / (9 * N_float_tot * N_float)) - (2 / N_float)) / (pow(an, 2) + bn);
     float uf = (float)((1 + ((N_float_tot + 1) / (3 * N_float)) - ((4 * ((N_float_tot + 1) / pow(N_float, 2))) * (an_plus_1 - ((2 * N_float_tot) / (N_float_tot + 1))))) / an) - vf;

     // corrected equations Simonsen et al 1995
     float vf_star = (float)((((2 * pow(N_float_tot, 3.0)) + (110 * pow(N_float_tot, 2.0)) - (225 * N_float_tot) + 153) / (9 * pow(N_float_tot, 2.0) * N_float)) + ((2 * N_float * an) / (pow(N_float_tot, 2.0))) - ((8 * bn) / N_float_tot)) / (pow(an, 2.0) + bn);
     float uf_star = (float)((((4 * pow(N_float_tot, 2.0)) + (19 * N_float_tot) + 3 - (12 * (N_float_tot + 1) * an_plus_1)) / (3 * N_float_tot * (N_float))) / an) - vf_star;

     vd_set.set_value(vd);
     ud_set.set_value(ud);
     vd_star_set.set_value(vd_star);
     ud_star_set.set_value(ud_star);
     uf_set.set_value(uf);
     vf_set.set_value(vf);
     uf_star_set.set_value(uf_star);
     vf_star_set.set_value(vf_star);
}