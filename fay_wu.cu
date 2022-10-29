#include "fay_wu.cuh"
#include "functions.cuh"
#include "prometheus.cuh"

fay_wu::fay_wu(string gene_List, string input_Folder, string ouput_Path, int cuda_ID, string intermediate_Path, int ploidy)
{
     /**
      * * Constructor Function
      * NORMAL - GENE MODE constructor
      **/

     cout << "Initiating CUDA powered Fay and Wu's normalized H and E calculator" << endl
          << endl;

     set_Values(gene_List, input_Folder, ouput_Path, cuda_ID, intermediate_Path, ploidy);
     // this->gene_List = gene_List;
     // cout << "Gene list file path\t: " << gene_List << endl
     //      << endl;
     // this->input_Folder = input_Folder;
     // this->ouput_Path = ouput_Path;
     // this->intermediate_Path = intermediate_Path;
     // this->ploidy = ploidy;

     // cudaSetDevice(cuda_ID);
     // cout << "Properties of selected CUDA GPU:" << endl;
     // cudaDeviceProp prop;
     // cudaGetDeviceProperties(&prop, cuda_ID);
     // cout << "GPU number\t: " << cuda_ID << endl;
     // cout << "GPU name\t: " << prop.name << endl;
     // size_t l_free = 0;
     // size_t l_Total = 0;
     // cudaError_t error_id = cudaMemGetInfo(&l_free, &l_Total);
     // cout << "GPU memory (GB)\t: " << l_Total / (1000 * 1000 * 1000) << endl;
     // cout << "GPU number of multiprocessor(s)\t: " << prop.multiProcessorCount << endl;
     // cout << "GPU block(s) per multiprocessor\t: " << prop.maxBlocksPerMultiProcessor << endl;
     // this->tot_Blocks = prop.maxBlocksPerMultiProcessor;
     // this->tot_ThreadsperBlock = prop.maxThreadsPerBlock;
     // cout << "GPU thread(s) per block\t: " << tot_ThreadsperBlock << endl
     //      << endl;
}

fay_wu::fay_wu(string gene_List, string input_Folder, string ouput_Path, int cuda_ID, string intermediate_Path, int ploidy, string prometheus_Activate, string Multi_read, int number_of_genes, int CPU_cores, int SNPs_per_Run)
{
     /**
      * * Constructor Function
      * PROMETHEUS - GENE MODE constructor
      **/

     // PROMETHEUS CONSTRUCTOR
     cout << "Initiating CUDA powered Fay and Wu's normalized H and E calculator on PROMETHEUS" << endl
          << endl;

     set_Values(gene_List, input_Folder, ouput_Path, cuda_ID, intermediate_Path, ploidy);

     this->prometheus_Activate = "YES";
     this->CPU_cores = CPU_cores;
     this->SNPs_per_Run = SNPs_per_Run;
     transform(Multi_read.begin(), Multi_read.end(), Multi_read.begin(), ::toupper);
     this->Multi_read = Multi_read;
     this->number_of_genes = number_of_genes;
}

fay_wu::fay_wu(string calc_Mode, int window_Size, int step_Size, string input_Folder, string ouput_Path, int cuda_ID, int ploidy, string prometheus_Activate, string Multi_read, int number_of_genes, int CPU_cores, int SNPs_per_Run)
{
     /**
      * * Constructor Function
      * PROMETHEUS - WINDOW MODE constructor
      **/

     // PROMETHEUS WINDOW MODE
     cout << "Initiating CUDA powered Fay and Wu's normalized H and E calculator on PROMETHEUS" << endl
          << endl;

     this->calc_Mode = "WINDOW";
     this->window_Size = window_Size;
     this->step_Size = step_Size;

     set_Values("", input_Folder, ouput_Path, cuda_ID, "", ploidy);

     this->prometheus_Activate = "YES";
     this->CPU_cores = CPU_cores;
     this->SNPs_per_Run = SNPs_per_Run;
     transform(Multi_read.begin(), Multi_read.end(), Multi_read.begin(), ::toupper);
     this->Multi_read = Multi_read;
     this->number_of_genes = number_of_genes;
}
fay_wu::fay_wu(string calc_Mode, int window_Size, int step_Size, string input_Folder, string ouput_Path, int cuda_ID, int ploidy)
{
     /**
      * * Constructor Function
      * NORMAL - WINDOW MODE constructor
      **/

     // NORMAL WINDOW CONSTRUCTOR

     cout << "Initiating CUDA powered Fay and Wu's normalized H and E calculator" << endl
          << endl;

     this->calc_Mode = "WINDOW";
     this->window_Size = window_Size;
     this->step_Size = step_Size;

     set_Values("", input_Folder, ouput_Path, cuda_ID, "", ploidy);
}

void fay_wu::set_Values(string gene_List, string input_Folder, string ouput_Path, int cuda_ID, string intermediate_Path, int ploidy)
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
     this->ouput_Path = ouput_Path;
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

void fay_wu::ingress()
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
          // first: start_stop second: filename
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
           * Pre-requisite values needed for determination of Fay and Wu.
           **/
          float an, bn, bn_plus1;
          calc_Pre(an, bn, bn_plus1, N);

          /**
           * @param test is used by Prometheus, to tell it which test is being processed.
           * * T  = Tajima
           * FU   = Fu and Li
           * FA   = Fay and Wu
           * N    = All 3 Neutrality tests
           **/
          string test = "FA";

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
               string output_File = ouput_Path + "/" +
                                    country.substr(country.find_last_of("/") + 1, country.length()) + "_" +
                                    to_string(window_Size) + "_" + to_string(step_Size) +
                                    ".fw";

               /**
                * Ensures if PROMETHEUS is being activated.
                **/
               if (prometheus_Activate == "YES")
               {
                    /**
                     * If Prometheus is being ACTIVATED then it is initialised accordingly.
                     **/
                    prometheus pro_Fay_Wu_Window = prometheus(output_File, window_Size, step_Size, folder_Index, Multi_read, tot_Blocks, tot_ThreadsperBlock, CPU_cores, SNPs_per_Run, number_of_genes, N, combinations, an, bn, bn_plus1);
                    /**
                     * Ensures if it is NORMAL window or SLIDING window mode.
                     * If step_Size is = 0 then it is sliding window mode.
                     **/
                    if (step_Size != 0)
                    {
                         /**
                          * Initiates processing of Tajima on PROMETHEUS on step wise window mode.
                          **/
                         pro_Fay_Wu_Window.process_Window(test);
                    }
                    else
                    {
                         /**
                          * Initiates processing of Tajima on PROMETHEUS on sliding window mode.
                          **/
                         pro_Fay_Wu_Window.process_C_sliding_Window(test);
                    }
               }
               else
               {
                    /**
                     * If Prometheus is NOT being activated the window calls be done accordingly.
                     **/
                    // Prometheus OFF Window Mode
                    if (step_Size != 0)
                    {
                         /**
                          * Initiates processing of Tajima on step wise window mode.
                          **/
                         window(output_File, an, bn, bn_plus1, N_float, combinations, folder_Index);
                    }
                    else
                    {
                         /**
                          * Initiates processing of Tajima on sliding window mode.
                          **/
                         window_Sliding(output_File, an, bn, bn_plus1, N_float, combinations, folder_Index);
                    }
               }
          }
          else
          {
               /**
                * * GENE (FILE) mode configuration:
                **/
               fstream gene_File;
               gene_File.open(gene_List, ios::in);
               cout << "Processing gene list:" << endl;
               /**
                * Output file is created for the population in the output folder for FILE mode.
                * @param output_File stores the output file's location.
                * The file name is a combination of the the country, and gene file name.
                **/
               string output_File = ouput_Path + "/" +
                                    country.substr(country.find_last_of("/") + 1, country.length()) + "_" +
                                    filesystem::path(gene_List).stem().string() +
                                    ".fw";
               /**
                * Log file created in the intermediate folder for the population.
                * @param intermediate_File stores the log file's location.
                * ! This helps with the resume function. Automatically resumes from the last completely processed gene in the event of a program crash.
                **/
               string intermediate_File = intermediate_Path + "/" +
                                          country.substr(country.find_last_of("/") + 1, country.length()) + "_" +
                                          filesystem::path(gene_List).stem().string() +
                                          ".log_fw";
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
                         function.createFile(output_File, "Gene_name\tCoordinates\tPi\tS\tTotal_iEi\tFay_Wu_Normalized_H\tFay_Wu_Normalized_E");
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

                    /**
                     * Ensures if PROMETHEUS is being activated.
                     **/
                    // PROMETHEUS HERE
                    if (prometheus_Activate == "YES")
                    {
                         cout << "Initializing Prometheus:" << endl
                              << endl;

                         /**
                          * If Prometheus is being ACTIVATED then it is initialised accordingly.
                          **/
                         prometheus pro_Fay_Wu = prometheus(folder_Index, Multi_read, tot_Blocks, tot_ThreadsperBlock, CPU_cores, SNPs_per_Run, number_of_genes, N, combinations, an, bn, bn_plus1);

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
                                   vector<string> write_Lines = pro_Fay_Wu.collection_Engine(gene_Collect, test);
                                   // print
                                   cout << "System is writing Fay and Wu results" << endl;
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
                                   pro_Fay_Wu.erase();
                                   gene_Collect.clear();
                                   cout << endl;
                              }
                         }
                         /**
                          * Ensures that there are no left over collected regions after finishing reading the gene file.
                          **/
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
                              vector<string> write_Lines = pro_Fay_Wu.collection_Engine(gene_Collect, test);
                              // print
                              cout << "System is writing Fay and Wu results" << endl;
                              /**
                               * Outputs are written and logs are made.
                               **/
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
                         pro_Fay_Wu.erase();
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
                              int start_Co = stoi(coordinates[1]);
                              /**
                               * @param end_Co captures query gene's end position as an integer.
                               **/
                              int end_Co = stoi(coordinates[2]);
                              cout << "Coordinates\t: Chromosome: " << coordinates[0] << " Start: " << start_Co << " End: " << end_Co << endl;

                              /**
                               * @param tot_pairwise_Differences Tajima's D also requires the tot_pairwise_Differences in the query region to determine the average number of pairwise differences in the region.
                               * @param segregating_Sites Tajima's D requires the total number of segregating sites/ SNPS in the query region.
                               **/
                              float tot_pairwise_Differences = 0;

                              /**
                               * The SNPs (Segregating sites) that fall within the query region are collected from the VCF's.
                               * @param collect_Segregrating_sites vector stores the collected SNPs.
                               **/
                              vector<string> collect_Segregrating_sites;

                              /**
                               * @param file_List vector is used to store the list of VCF files (found via CATES CIS algorithm) that satisfy the query region.
                               **/
                              vector<string> file_List;
                              cout << endl;
                              cout << "System is retrieving file(s)" << endl;
                              if (folder_Index.size() > 1)
                              {
                                   file_List = function.compound_interpolationSearch(folder_Index, start_Co, end_Co);
                              }
                              else
                              {
                                   /**
                                    * IF only one file is present in the index folder that file will be used as is.
                                    **/
                                   file_List.push_back(folder_Index[0].second);
                              }
                              cout << "System has retrieved all file(s)" << endl;
                              cout << endl;

                              cout << "System is collecting segregrating site(s)" << endl;

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
                                        getline(file, line); // skip first header line
                                        while (getline(file, line))
                                        {
                                             /**
                                              * @param positions vector is used to capture the SNP data upto the position column (Column 2 (non zero count)).
                                              **/
                                             vector<string> positions;
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

                                                  // string check_0 = country.substr(country.find_last_of("/") + 1, country.length()) + "_AF=0";
                                                  // string seg_Check = "GO";
                                                  // vector<string> info;
                                                  // function.split(info, positions[7], ";");
                                                  // for (string AF_check : info)
                                                  // {
                                                  //      if (AF_check == check_0)
                                                  //      {
                                                  //           seg_Check = "NO";
                                                  //           break;
                                                  //      }
                                                  // }
                                                  // if (seg_Check == "GO")
                                                  // {
                                                  //      string check_AF_country = country.substr(country.find_last_of("/") + 1, country.length()) + "_AF";
                                                  //      float MAF_float = 0.0000;
                                                  //      // collect_Segregrating_sites.push_back(line);
                                                  //      for (string AF_check : info)
                                                  //      {
                                                  //           vector<string> split_info;
                                                  //           function.split(split_info, AF_check, "=");
                                                  //           if (split_info[0] == check_AF_country)
                                                  //           {
                                                  //                MAF_float = stof(split_info[1]);
                                                  //                if (MAF_float > 0.5)
                                                  //                {
                                                  //                     MAF_float = 1 - MAF_float;
                                                  //                }
                                                  //                break;
                                                  //           }
                                                  //      }
                                                  //      tot_pairwise_Differences = tot_pairwise_Differences + (MAF_float * (1 - MAF_float) * pow(N_float, 2));
                                                  // }
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

                              int num_segregrating_Sites;
                              string Fay_Wu_H, Fay_Wu_E;
                              float pi = 0.0;
                              int Total_iEi = 0;

                              /**
                               * Calls the function to process the segregating sites and calculate theta L requried for calculating the Fay and Wu values.
                               **/
                              float theta_L = calc_theta_L(collect_Segregrating_sites, N_float, num_segregrating_Sites, Total_iEi, tot_pairwise_Differences);

                              cout << "Total segregating sites (S)\t: " << num_segregrating_Sites << endl;
                              cout << endl;
                              if (num_segregrating_Sites != 0)
                              {
                                   float S = (float)num_segregrating_Sites;
                                   float theta_squared = (float)(S * (S - 1)) / (pow(an, 2) + bn);
                                   cout << "Theta_squared\t: " << theta_squared << endl;
                                   cout << "Theta_L\t: " << theta_L << endl;
                                   float theta_W = (float)S / an;
                                   cout << "Theta_W\t: " << theta_W << endl;
                                   pi = (float)tot_pairwise_Differences / combinations;
                                   cout << "Average pairwise polymorphisms (pi)\t: " << pi << endl;
                                   cout << endl;

                                   float VAR_pi_MINUS_theta_L = (float)(((N_float - 2.0) / (6.0 * (N_float - 1.0))) * theta_W) + ((((18.0 * pow(N_float, 2) * ((3.0 * N_float) + 2.0) * bn_plus1) - ((88.0 * pow(N_float, 3)) + (9.0 * pow(N_float, 2)) - (13.0 * N_float) + 6.0)) / (9.0 * N_float * pow(N_float - 1, 2))) * theta_squared);
                                   // cout << "VAR_pi_MINUS_theta_L: " << VAR_pi_MINUS_theta_L << endl;
                                   float VAR_theta_L_MINUS_theta_W = (float)(((N_float / (2.0 * (N_float - 1.0))) - (1.0 / an)) * theta_W) + (((bn / (pow(an, 2))) + (2.0 * pow((N_float / (N_float - 1.0)), 2) * bn) - ((2.0 * ((N_float * bn) - N_float + 1.0)) / ((N_float - 1.0) * an)) - (((3.0 * N_float) + 1) / (N_float - 1.0))) * theta_squared);
                                   // cout << "VAR_theta_L_MINUS_theta_W: " << VAR_theta_L_MINUS_theta_W << endl;

                                   float H = (float)(pi - theta_L) / (sqrt(VAR_pi_MINUS_theta_L));
                                   Fay_Wu_H = to_string(H);
                                   cout << "Fay and Wu's normalized H\t: " << Fay_Wu_H << endl;
                                   float E = (float)(theta_L - theta_W) / (sqrt(VAR_theta_L_MINUS_theta_W));
                                   Fay_Wu_E = to_string(E);
                                   cout << "Fay and Wu's normalized E\t: " << Fay_Wu_E << endl;
                              }
                              else
                              {
                                   cout << "Fay and Wu's H and E\t: "
                                        << "Not Available" << endl;
                                   Fay_Wu_H = "NA";
                                   Fay_Wu_E = "NA";
                              }

                              cout << endl;

                              // Gene_name\tCoordinates\tPi\tS\tTotal_iEi\tFay_Wu_Normalized_H\tFay_Wu_Normalized_E
                              output << gene_Name << "\t"
                                     << coordinates[0] << ":" << to_string(start_Co) << ":" << to_string(end_Co)
                                     << "\t" << to_string(pi)
                                     << "\t" << to_string(num_segregrating_Sites)
                                     << "\t" << to_string(Total_iEi)

                                     << "\t" << Fay_Wu_H
                                     << "\t" << Fay_Wu_E << "\n";

                              intermediate << gene_Combo << "\n";
                              output.flush();
                              intermediate.flush();
                         }
                    }
                    output.close();
                    intermediate.close();
                    gene_File.close();
               }
          }
     }
}

void fay_wu::window_Sliding(string output_File, float an, float bn, float bn_plus1, float N_float, long int combinations, vector<pair<string, string>> &folder_Index)
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
          function.createFile(output_File, "Coordinates\tPi\tS\tTotal_iEi\tFay_Wu_Normalized_H\tFay_Wu_Normalized_E");
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
                                        /**
                                         * If they do not match,
                                         * found is set to 1 indicating the resume position has been found and,
                                         * the loop is broken.
                                         **/
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

                              /**
                               * @param tot_pairwise_Differences Fay and Wu also requires the tot_pairwise_Differences in the query region to determine the average number of pairwise differences in the region.
                               **/
                              float tot_pairwise_Differences = 0;

                              /**
                               * @param file_List vector is used to store the list of VCF files (found via CATES CIS algorithm) that satisfy the query region.
                               **/
                              vector<string> file_List;
                              cout << endl;
                              cout << "System is retrieving file(s)" << endl;
                              if (folder_Index.size() > 1)
                              {
                                   file_List = function.compound_interpolationSearch(folder_Index, start_Co, end_Co);
                              }
                              else
                              {
                                   /**
                                    * IF only one file is present in the index folder that file will be used as is.
                                    **/
                                   file_List.push_back(folder_Index[0].second);
                              }
                              cout << "System has retrieved all file(s)" << endl;

                              cout << "System is collecting segregrating site(s)" << endl;
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
                                   // cout << files << endl;
                                   fstream file;
                                   file.open(files, ios::in);
                                   if (file.is_open())
                                   {
                                        string line;
                                        /**
                                         * The first line of each VCF is skipped as it is the header line.
                                         **/
                                        getline(file, line); // skip first header line
                                        while (getline(file, line))
                                        {
                                             /**
                                              * @param positions vector is used to capture the SNP data upto the position column (Column 2 (non zero count)).
                                              **/
                                             vector<string> positions;
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

                              int num_segregrating_Sites;
                              string Fay_Wu_H, Fay_Wu_E;
                              float pi = 0.0;
                              int Total_iEi = 0;

                              /**
                               * Calls the function to process the segregating sites and calculate theta L requried for calculating the Fay and Wu values.
                               **/
                              float theta_L = calc_theta_L(collect_Segregrating_sites, N_float, num_segregrating_Sites, Total_iEi, tot_pairwise_Differences);

                              cout << "Total segregating sites (S)\t: " << num_segregrating_Sites << endl;
                              cout << endl;

                              if (num_segregrating_Sites != 0)
                              {
                                   float S = (float)num_segregrating_Sites;
                                   float theta_squared = (float)(S * (S - 1)) / (pow(an, 2) + bn);
                                   cout << "Theta_squared\t: " << theta_squared << endl;
                                   cout << "Theta_L\t: " << theta_L << endl;
                                   float theta_W = (float)S / an;
                                   cout << "Theta_W\t: " << theta_W << endl;
                                   pi = (float)tot_pairwise_Differences / combinations;
                                   cout << "Average pairwise polymorphisms (pi)\t: " << pi << endl;
                                   cout << endl;

                                   float VAR_pi_MINUS_theta_L = (float)(((N_float - 2.0) / (6.0 * (N_float - 1.0))) * theta_W) + ((((18.0 * pow(N_float, 2) * ((3.0 * N_float) + 2.0) * bn_plus1) - ((88.0 * pow(N_float, 3)) + (9.0 * pow(N_float, 2)) - (13.0 * N_float) + 6.0)) / (9.0 * N_float * pow(N_float - 1, 2))) * theta_squared);
                                   // cout << "VAR_pi_MINUS_theta_L: " << VAR_pi_MINUS_theta_L << endl;
                                   float VAR_theta_L_MINUS_theta_W = (float)(((N_float / (2.0 * (N_float - 1.0))) - (1.0 / an)) * theta_W) + (((bn / (pow(an, 2))) + (2.0 * pow((N_float / (N_float - 1.0)), 2) * bn) - ((2.0 * ((N_float * bn) - N_float + 1.0)) / ((N_float - 1.0) * an)) - (((3.0 * N_float) + 1) / (N_float - 1.0))) * theta_squared);
                                   // cout << "VAR_theta_L_MINUS_theta_W: " << VAR_theta_L_MINUS_theta_W << endl;

                                   float H = (float)(pi - theta_L) / (sqrt(VAR_pi_MINUS_theta_L));
                                   Fay_Wu_H = to_string(H);
                                   cout << "Fay and Wu's normalized H\t: " << Fay_Wu_H << endl;
                                   float E = (float)(theta_L - theta_W) / (sqrt(VAR_theta_L_MINUS_theta_W));
                                   Fay_Wu_E = to_string(E);
                                   cout << "Fay and Wu's normalized E\t: " << Fay_Wu_E << endl;
                              }
                              else
                              {
                                   cout << "Fay and Wu's H and E\t: "
                                        << "Not Available" << endl;
                                   Fay_Wu_H = "NA";
                                   Fay_Wu_E = "NA";
                              }

                              cout << endl;

                              output << to_string(start_Co) << ":" << to_string(end_Co)
                                     << "\t" << to_string(pi)
                                     << "\t" << to_string(num_segregrating_Sites)
                                     << "\t" << to_string(Total_iEi)

                                     << "\t" << Fay_Wu_H
                                     << "\t" << Fay_Wu_E << "\n";

                              output.flush();

                              start_Co = start_Co + step_Size;
                              end_Co = start_Co + window_Size;
                         }
                    }
               }

               file_Main.close();
          }
     }

     output.close();
}

void fay_wu::window(string output_File, float an, float bn, float bn_plus1, float N_float, long int combinations, vector<pair<string, string>> &folder_Index)
{
     /**
      * NORMAL MODE WINDOW FUNCTION
      **/

     /**
      * Call the "functions" class. Bespoke functions commonly used by CATE.
      **/
     functions function = functions();

     /**
      * @param start_Value is used to capture the lowest POSITION available in the VCF range.
      * @param end_Value is used to capture the highest POSITION available in the VCF range.
      **/
     int start_Value = stoi(folder_Index[0].first.substr(0, folder_Index[0].first.find('_')));
     int end_Value = stoi(folder_Index[folder_Index.size() - 1].first.substr(folder_Index[folder_Index.size() - 1].first.find('_') + 1));

     // cout << start_Value << endl;
     // cout << end_Value << endl;

     /**
      * @param start_Co captures the start position as an integer.
      * @param end_Co captures the end position as an integer.
      **/
     int start_Co = 0;
     int end_Co = start_Co + window_Size;

     /**
      * We cycle through till we fall into the range of SNPs available in our VCFs.
      **/
     while (start_Value > end_Co)
     {
          start_Co = start_Co + step_Size;
          end_Co = start_Co + window_Size;
     }

     cout << "Writing to file\t: " << output_File << endl;
     cout << endl;

     /**
      * If the output file is absent this run will be considered as a brand new run of this query and,
      * the output file and the intermediate log file will be created.
      **/
     if (filesystem::exists(output_File) == 0)
     {
          function.createFile(output_File, "Coordinates\tPi\tS\tTotal_iEi\tFay_Wu_Normalized_H\tFay_Wu_Normalized_E");
     }
     else
     {
          // RESUME FUNCTION
          /**
           * If the output file is already present then the resume process will initiated.
           * This is a unintelligent resume. Essentially it matches the each read line written with the lines read from the gene file.
           * The break will occur as soon as their is a mismatch.
           * To counter any errors it is advised to have a new gene file name or a new intermediate folder per new run.
           * @param caught acts as a boolean variable. caught = 0 if the lines need to be skipped and will equal 1 when the resume position is found.
           * ! This helps with the resume function. Automatically resumes from the last completely processed gene in the event of a program crash.
           **/
          int caught = 0;

          // skipper
          /**
           * Open the output file to be begin finding the resume point.
           **/
          fstream output;
          output.open(output_File, ios::in);

          string output_Line;

          while (start_Co <= end_Value)
          {
               /**
                * First line is skipped cause it is a header line containing column names.
                **/
               // skip header
               getline(output, output_Line);

               while (getline(output, output_Line))
               {
                    string trim = output_Line.substr(0, output_Line.find('\t'));
                    string check = to_string(start_Co) + ":" + to_string(end_Co);
                    /**
                     * Ensures the query line does not match that of the output
                     **/
                    if (trim != check)
                    {
                         /**
                          * If they do not match,
                          * caught is set to 1 indicating the resume position has been found and,
                          * the loop is broken.
                          **/
                         caught = 1;
                         break;
                    }
               }
               /**
                * If caught is 1 that means the resume location has been found and the loop is broken.
                **/
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

          /**
           * @param tot_pairwise_Differences Fay and Wu also requires the tot_pairwise_Differences in the query region to determine the average number of pairwise differences in the region.
           * @param segregating_Sites Fay and Wu requires the total number of segregating sites/ SNPS in the query region.
           **/
          float tot_pairwise_Differences = 0;

          /**
           * The SNPs (Segregating sites) that fall within the query region are collected from the VCF's.
           * @param collect_Segregrating_sites vector stores the collected SNPs.
           **/
          vector<string> collect_Segregrating_sites;

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
          cout << endl;

          cout << "System is collecting segregrating site(s)" << endl;

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
                    getline(file, line); // skip first header line
                    while (getline(file, line))
                    {
                         vector<string> positions;
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

          int num_segregrating_Sites;
          string Fay_Wu_H, Fay_Wu_E;
          float pi = 0.0;
          int Total_iEi = 0;

          /**
           * Calls the function to process the segregating sites and calculate theta L requried for calculating the Fay and Wu values.
           **/
          float theta_L = calc_theta_L(collect_Segregrating_sites, N_float, num_segregrating_Sites, Total_iEi, tot_pairwise_Differences);

          cout << "Total segregating sites (S)\t: " << num_segregrating_Sites << endl;
          cout << endl;

          if (num_segregrating_Sites != 0)
          {
               float S = (float)num_segregrating_Sites;
               float theta_squared = (float)(S * (S - 1)) / (pow(an, 2) + bn);
               cout << "Theta_squared\t: " << theta_squared << endl;
               cout << "Theta_L\t: " << theta_L << endl;
               float theta_W = (float)S / an;
               cout << "Theta_W\t: " << theta_W << endl;
               pi = (float)tot_pairwise_Differences / combinations;
               cout << "Average pairwise polymorphisms (pi)\t: " << pi << endl;
               cout << endl;

               float VAR_pi_MINUS_theta_L = (float)(((N_float - 2.0) / (6.0 * (N_float - 1.0))) * theta_W) + ((((18.0 * pow(N_float, 2) * ((3.0 * N_float) + 2.0) * bn_plus1) - ((88.0 * pow(N_float, 3)) + (9.0 * pow(N_float, 2)) - (13.0 * N_float) + 6.0)) / (9.0 * N_float * pow(N_float - 1, 2))) * theta_squared);
               // cout << "VAR_pi_MINUS_theta_L: " << VAR_pi_MINUS_theta_L << endl;
               float VAR_theta_L_MINUS_theta_W = (float)(((N_float / (2.0 * (N_float - 1.0))) - (1.0 / an)) * theta_W) + (((bn / (pow(an, 2))) + (2.0 * pow((N_float / (N_float - 1.0)), 2) * bn) - ((2.0 * ((N_float * bn) - N_float + 1.0)) / ((N_float - 1.0) * an)) - (((3.0 * N_float) + 1) / (N_float - 1.0))) * theta_squared);
               // cout << "VAR_theta_L_MINUS_theta_W: " << VAR_theta_L_MINUS_theta_W << endl;

               float H = (float)(pi - theta_L) / (sqrt(VAR_pi_MINUS_theta_L));
               Fay_Wu_H = to_string(H);
               cout << "Fay and Wu's normalized H\t: " << Fay_Wu_H << endl;
               float E = (float)(theta_L - theta_W) / (sqrt(VAR_theta_L_MINUS_theta_W));
               Fay_Wu_E = to_string(E);
               cout << "Fay and Wu's normalized E\t: " << Fay_Wu_E << endl;
          }
          else
          {
               cout << "Fay and Wu's H and E\t: "
                    << "Not Available" << endl;
               Fay_Wu_H = "NA";
               Fay_Wu_E = "NA";
          }

          cout << endl;

          output << to_string(start_Co) << ":" << to_string(end_Co)
                 << "\t" << to_string(pi)
                 << "\t" << to_string(num_segregrating_Sites)
                 << "\t" << to_string(Total_iEi)

                 << "\t" << Fay_Wu_H
                 << "\t" << Fay_Wu_E << "\n";

          output.flush();

          start_Co = start_Co + step_Size;
          end_Co = start_Co + window_Size;
     }
     output.close();
}

__global__ void cuda_theta_L(char *sites, int *index, int num_Segregrating_sites, int *theta_Partials, int *VALID_or_NOT, int *MA_count)
{
     /**
      * @param tid is used to get the unique thread ID. In this instance thread ID is used to keep track of the Seg/SNP sites.
      **/
     int tid = threadIdx.x + blockIdx.x * blockDim.x;

     while (tid < num_Segregrating_sites)
     {
          /**
           * We like to think CATE strides or skates through the 1D SNP array based on the SEG site the individual thread is assigned.
           * Each thread has a boundary wall within which it will stay. Similar to a skater's rink.
           * This wall is defined by two variables, namely:
           * @param site_Start is used to define the start of the Seg site assigned to the thread.
           * @param site_End is used to define the end of the Seg site assigned to the thread.
           * The thread will move across this region extracting information as needed.
           **/

          /**
           * @param column is used to keep track of the columns being strided through.
           * Since VCF's are tab eliminated we are abe to track the columns by keeping track of the '\t' we come across.
           * VCF's by default have 9 columns. Beyond these 9 are information of the samples that have been sequenced.
           * [0] = chromosome number or index
           * [1] = position
           * [2] = snp id
           * [3] = reference allele
           * [4] = alternate allele
           * [5] = quality
           * [6] = filter
           * [7] = info
           * [8] = format
           * [9 ..] = sample information
           **/
          int column = 0;
          int site_Start = index[tid];
          int site_End = index[tid + 1];

          /**
           * @param i is used to track the navigation through the SNP information. This is our skater, an increment variable.
           **/
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

          /**
           * @param REF is used to capture the Reference allele from column 3 (column index above). It is converted to uppercase if in lowercase.
           * @param ALT is used to capture the Reference allele from column 4 (column index above). It is converted to uppercase if in lowercase.
           **/

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

          /**
           * @param AA is used to capture the Ancestral Allele (AA) information from the AA tag in the INFO column 7 (column index above).
           **/

          char AA = 'N';

          while (column < 8)
          {
               if ((sites[i] == 'A' && sites[i + 1] == 'A' && sites[i + 2] == '=') || (sites[i] == 'a' && sites[i + 1] == 'a' && sites[i + 2] == '='))
               {
                    char CHECK = 'N';

                    /**
                     * If the AA data is present the ancestral allele is validated and converted to uppercase.
                     **/

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

          // printf("AA: %c\n", AA);
          // printf("\n");

          // printf("%c", sites[i]);

          while (column < 9)
          {
               if (sites[i] == '\t')
               {
                    column++;
               }
               i++;
          }

          // printf("Column 1: %c\n", sites[i]);

          /**
           * @param ALT_count is used capture number of instances the ALTERNATE allele is present.
           * @param REF_count is used capture number of instances the REFERENCE allele is present.
           **/
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

          /**
           * If either allele has a count of zero this means that the site is not a SEG site for that population.
           **/
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

               int theta_Partial = 0;
               if (MA != AA)
               {
                    theta_Partial = MA_count[tid];
               }

               // printf("theta partial: %d\n", theta_Partial);

               theta_Partials[tid] = theta_Partial;
          }

          // char MA = 'N';
          // int MA_count = 0;
          // if (ALT_count < REF_count)
          // {
          //      MA = ALT;
          //      MA_count = ALT_count;
          // }
          // else
          // {
          //      MA = REF;
          //      MA_count = REF_count;
          // }

          // //printf("Minor allele: %c\n", MA);
          // //printf("Minor allele count: %d\n", MA_count);

          // int theta_Partial = 0;
          // if (MA != AA)
          // {
          //      theta_Partial = MA_count;
          // }

          // //printf("theta partial: %d\n", theta_Partial);

          // theta_Partials[tid] = theta_Partial;

          tid += blockDim.x * gridDim.x;
     }
}

float fay_wu::calc_theta_L(vector<string> &total_Segregrating_sites, float N_tot, int &real_segregrating_Sites, int &Total_iEi, float &tot_pairwise_Differences)
{
     /**
      * Administrative function responsible for collection of data required for Tajima's D.
      * 1. Conversion of SNP strings into char pointers for GPU accessability.
      * 2. Call GPU for extracting MAs (Minor allele) and MAF's (Minor Allele Frequencies).
      * 3. Calculates the prerequisites required for determining Theta L by accounting for the alleles that are not equal to the AA.
      **/
     cout << "System is processing and filtering segregrating site(s)" << endl;

     /**
      * @param num_segregrating_Sites is used to track the number of SNPs collected for the query region.
      * This track is vital for navigating through the data in the GPU. For the data is stored in the form of a 1D array.
      **/
     int num_segregrating_Sites = total_Segregrating_sites.size();

     float theta_L = 0.00;

     /**
      * @param Seg_sites is used to stitch the SNP data end to end, before converting it to a char array.
      **/
     string Seg_sites = "";

     /**
      * @param site_Index is used to keep track of the start and ends of each SNP's data.
      **/
     int site_Index[num_segregrating_Sites + 1];
     site_Index[0] = 0;
     // cout << "copy string" << endl;

     /**
      * Conversion of vector SNP information into a 1D array by concatenating the vector data into a single string.
      **/
     for (size_t i = 0; i < num_segregrating_Sites; i++)
     {
          Seg_sites.append(total_Segregrating_sites[i]);
          site_Index[i + 1] = site_Index[i] + total_Segregrating_sites[i].size();
     }
     // cout << "copy string done" << endl;
     // cout << "Seg size : " << Seg_sites.size() << endl;

     /**
      * Final assignment of concatented string into a 1D char array.
      **/
     char *full_Char;
     full_Char = (char *)malloc((Seg_sites.size() + 1) * sizeof(char));
     // cout << "copy char " << endl;
     strcpy(full_Char, Seg_sites.c_str());
     // cout << "copy char done" << endl;

     /**
      * RAM is released to prevent redundancy.
      **/
     total_Segregrating_sites.clear();

     /**
      * @param cuda_full_Char is used by the GPU. Is a COPY of full_Char.
      * @param cuda_site_Index is used by the GPU. Is a COPY of site_Index.
      *
      * * These 4 variables work together in all instances.
      **/
     char *cuda_full_Char;
     cudaMallocManaged(&cuda_full_Char, (Seg_sites.size() + 1) * sizeof(char));
     int *cuda_site_Index;
     cudaMallocManaged(&cuda_site_Index, (num_segregrating_Sites + 1) * sizeof(int));

     /**
      * Certain collected SNP sites though present in the VCF, specially if the index was created without a filter,
      * will have MAF of zero. Indicating that region is conserved for that super-pop and therfore not a seg site.
      * These have to be filtered out. This is done using the variables below:
      * @param cuda_VALID_or_NOT is used to determine if a site is seg site which is VALID or NOT.
      * @param VALID_or_NOT is used by the CPU. Is a COPY of cuda_VALID_or_NOT.
      **/
     int *cuda_VALID_or_NOT, *VALID_or_NOT;
     cudaMallocManaged(&cuda_VALID_or_NOT, num_segregrating_Sites * sizeof(int));
     VALID_or_NOT = (int *)malloc(num_segregrating_Sites * sizeof(int));

     /**
      * @param cuda_MA_Count is used to record the MA's count.
      * @param MA_Count is used by the CPU. Is a COPY of cuda_MA_Count.
      **/
     int *cuda_MA_Count, *MA_Count;
     cudaMallocManaged(&cuda_MA_Count, num_segregrating_Sites * sizeof(int));
     MA_Count = (int *)malloc(num_segregrating_Sites * sizeof(int));

     int *cuda_Theta_partials, *Theta_partials;
     cudaMallocManaged(&cuda_Theta_partials, num_segregrating_Sites * sizeof(int));
     Theta_partials = (int *)malloc(num_segregrating_Sites * sizeof(int));

     /**
      * Transfer of data to the GPU.
      **/
     cudaMemcpy(cuda_full_Char, full_Char, (Seg_sites.size() + 1) * sizeof(char), cudaMemcpyHostToDevice);
     cudaMemcpy(cuda_site_Index, site_Index, (num_segregrating_Sites + 1) * sizeof(int), cudaMemcpyHostToDevice);
     // cout << "GPU" << endl;

     /**
      * CALL THE GPU.
      * * GPU WILL PROCESS THE COLLECTED SEG SITES
      **/
     cuda_theta_L<<<tot_Blocks, tot_ThreadsperBlock>>>(cuda_full_Char, cuda_site_Index, num_segregrating_Sites, cuda_Theta_partials, cuda_VALID_or_NOT, cuda_MA_Count);
     cudaDeviceSynchronize();

     cudaMemcpy(VALID_or_NOT, cuda_VALID_or_NOT, num_segregrating_Sites * sizeof(int), cudaMemcpyDeviceToHost);
     cudaMemcpy(MA_Count, cuda_MA_Count, num_segregrating_Sites * sizeof(int), cudaMemcpyDeviceToHost);
     cudaMemcpy(Theta_partials, cuda_Theta_partials, num_segregrating_Sites * sizeof(int), cudaMemcpyDeviceToHost);

     free(full_Char);
     cudaFree(cuda_full_Char);
     cudaFree(cuda_site_Index);
     cudaFree(cuda_MA_Count);
     cudaFree(cuda_VALID_or_NOT);
     cudaFree(cuda_Theta_partials);

     // replace with CUDA addition
     real_segregrating_Sites = 0;
     tot_pairwise_Differences = 0;
     int total_iTheta = 0;

     for (size_t i = 0; i < num_segregrating_Sites; i++)
     {
          if (VALID_or_NOT[i] == 1)
          {
               real_segregrating_Sites = real_segregrating_Sites + 1;
               float MAF = (float)MA_Count[i] / N_tot;
               tot_pairwise_Differences = tot_pairwise_Differences + (MAF * (1 - MAF) * pow(N_tot, 2));
               total_iTheta = total_iTheta + Theta_partials[i];
          }
     }
     // cout << "GPU DONE" << endl;
     // cout << "total iTheta: " << total_iTheta << endl;

     Total_iEi = total_iTheta;
     theta_L = (float)(1 / (N_tot - 1)) * (float)total_iTheta;

     // cout << "Theta_L: " << theta_L << endl;

     free(Theta_partials);
     free(MA_Count);
     free(VALID_or_NOT);

     cout << "System has completed processing the segregrating site(s)" << endl;
     return theta_L;
}

__global__ void faywu_Calculation(int N, float *a1_CUDA, float *a2_CUDA)
{
     int tid = threadIdx.x + blockIdx.x * blockDim.x;

     while (tid < N)
     {
          a1_CUDA[tid] = (float)1 / (tid + 1);
          a2_CUDA[tid] = (float)1 / ((tid + 1) * (tid + 1));
          tid += blockDim.x * gridDim.x;
     }
}

void fay_wu::calc_Pre(float &an, float &bn, float &bn_plus1, int N_tot)
{
     /**
      * Calculates the prerequisite values required for Fay and Wu
      **/
     functions function = functions();

     int N = N_tot - 1;
     float *a1_CUDA, *a2_CUDA;
     float *a1_partial, *a2_partial;

     a1_partial = (float *)malloc(N * sizeof(float));
     a2_partial = (float *)malloc(N * sizeof(float));

     cudaMallocManaged(&a1_CUDA, N * sizeof(int));
     cudaMallocManaged(&a2_CUDA, N * sizeof(int));

     faywu_Calculation<<<tot_Blocks, tot_ThreadsperBlock>>>(N, a1_CUDA, a2_CUDA);
     cudaDeviceSynchronize();

     cudaMemcpy(a1_partial, a1_CUDA, N * sizeof(float), cudaMemcpyDeviceToHost);
     cudaMemcpy(a2_partial, a2_CUDA, N * sizeof(float), cudaMemcpyDeviceToHost);

     cudaFree(a1_CUDA);
     cudaFree(a2_CUDA);

     an = function.add(N, a1_partial);
     bn = function.add(N, a2_partial);
     float N_tot_float = (float)N_tot;
     bn_plus1 = bn + (1.0 / pow(N_tot_float, 2));

     free(a1_partial);
     free(a2_partial);

     cout << "an: " << an << "\t"
          << "bn: " << bn << "\t"
          << "bn_plus_1: " << bn_plus1 << endl;

     cout << endl;
}