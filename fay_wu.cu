#include "fay_wu.cuh"
#include "functions.cuh"
#include "prometheus.cuh"

fay_wu::fay_wu(string gene_List, string input_Folder, string ouput_Path, int cuda_ID, string intermediate_Path, int ploidy)
{
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
     functions function = functions();

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
          cout << "Processing country\t: " << country.substr(country.find_last_of("/") + 1, country.length()) << endl
               << endl;

          vector<pair<string, string>> folder_Index = function.index_Folder(country);
          cout << "Completed indexing folder\t: " << country << endl;

          cout << endl;

          int samples = function.getN_Split(folder_Index[0].second);
          cout << "Number of samples in " << country.substr(country.find_last_of("/") + 1, country.length()) << " population\t: " << samples << endl;

          int N = samples * ploidy;
          float N_float = (float)N;
          cout << "Number of sequences in " << country.substr(country.find_last_of("/") + 1, country.length()) << " population [ " << samples << " x " << ploidy << " ] (N)\t: " << N << endl;

          long int combinations = function.combos_N(N);
          cout << "Pairwise combinations\t: " << combinations << endl;

          cout << endl;

          float an, bn, bn_plus1;
          calc_Pre(an, bn, bn_plus1, N);
          string test = "FA";

          if (this->calc_Mode != "FILE")
          {
               string output_File = ouput_Path + "/" +
                                    country.substr(country.find_last_of("/") + 1, country.length()) + "_" +
                                    to_string(window_Size) + "_" + to_string(step_Size) +
                                    ".fw";
               if (prometheus_Activate == "YES")
               {
                    prometheus pro_Fay_Wu_Window = prometheus(output_File, window_Size, step_Size, folder_Index, Multi_read, tot_Blocks, tot_ThreadsperBlock, CPU_cores, SNPs_per_Run, number_of_genes, N, combinations, an, bn, bn_plus1);
                    if (step_Size != 0)
                    {
                         pro_Fay_Wu_Window.process_Window(test);
                    }
                    else
                    {
                         pro_Fay_Wu_Window.process_C_sliding_Window(test);
                    }
               }
               else
               {
                    // Prometheus OFF Window Mode
                    if (step_Size != 0)
                    {
                         window(output_File, an, bn, bn_plus1, N_float, combinations, folder_Index);
                    }
                    else
                    {
                         window_Sliding(output_File, an, bn, bn_plus1, N_float, combinations, folder_Index);
                    }
               }
          }
          else
          {
               fstream gene_File;
               gene_File.open(gene_List, ios::in);
               cout << "Processing gene list:" << endl;
               string output_File = ouput_Path + "/" +
                                    country.substr(country.find_last_of("/") + 1, country.length()) + "_" +
                                    filesystem::path(gene_List).stem().string() +
                                    ".fw";
               string intermediate_File = intermediate_Path + "/" +
                                          country.substr(country.find_last_of("/") + 1, country.length()) + "_" +
                                          filesystem::path(gene_List).stem().string() +
                                          ".log_fw";
               cout << endl;
               cout << "Writing to file\t: " << output_File << endl;
               cout << endl;

               if (gene_File.is_open())
               {
                    string gene_Combo;

                    if (filesystem::exists(output_File) == 0)
                    {
                         function.createFile(output_File, "Gene_name\tCoordinates\tPi\tS\tTotal_iEi\tFay_Wu_Normalized_H\tFay_Wu_Normalized_E");
                         function.createFile(intermediate_File);
                    }
                    else
                    {
                         fstream intermediate;
                         intermediate.open(intermediate_File, ios::in);
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
                    if (prometheus_Activate == "YES")
                    {
                         cout << "Initializing Prometheus:" << endl
                              << endl;

                         prometheus pro_Fay_Wu = prometheus(folder_Index, Multi_read, tot_Blocks, tot_ThreadsperBlock, CPU_cores, SNPs_per_Run, number_of_genes, N, combinations, an, bn, bn_plus1);
                         vector<string> gene_Collect;

                         while (getline(gene_File, gene_Combo))
                         {
                              gene_Collect.push_back(gene_Combo);
                              if (gene_Collect.size() == number_of_genes)
                              {
                                   cout << "Prometheus batch intialized" << endl;
                                   cout << "From: " << gene_Collect[0] << endl;
                                   cout << "To  : " << gene_Collect[gene_Collect.size() - 1] << endl
                                        << endl;
                                   // launch prometheus
                                   vector<string> write_Lines = pro_Fay_Wu.collection_Engine(gene_Collect, test);
                                   // print
                                   cout << "System is writing Fay and Wu results" << endl;
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
                         if (gene_Collect.size() != 0)
                         {
                              // RUN PROMETHEUS for remaining
                              // launch prometheus
                              cout << "Prometheus batch intialized" << endl;
                              cout << "From: " << gene_Collect[0] << endl;
                              cout << "To  : " << gene_Collect[gene_Collect.size() - 1] << endl
                                   << endl;

                              vector<string> write_Lines = pro_Fay_Wu.collection_Engine(gene_Collect, test);
                              // print
                              cout << "System is writing Fay and Wu results" << endl;
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
                         while (getline(gene_File, gene_Combo))
                         {
                              vector<string> split_Data;
                              function.split(split_Data, gene_Combo, '\t');
                              string gene_Name = split_Data[0];
                              cout << "Gene name\t: " << gene_Name << endl;
                              vector<string> coordinates;
                              function.split(coordinates, split_Data[1], ':');
                              int start_Co = stoi(coordinates[1]);
                              int end_Co = stoi(coordinates[2]);
                              cout << "Coordinates\t: Chromosome: " << coordinates[0] << " Start: " << start_Co << " End: " << end_Co << endl;

                              float tot_pairwise_Differences = 0;
                              vector<string> collect_Segregrating_sites;

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
                              cout << endl;

                              cout << "System is collecting segregrating site(s)" << endl;
                              for (string files : file_List)
                              {
                                   fstream file;
                                   file.open(files, ios::in);
                                   if (file.is_open())
                                   {
                                        string line;
                                        getline(file, line); // skip first header line
                                        while (getline(file, line))
                                        {
                                             vector<string> positions;
                                             function.split_getPos_ONLY(positions, line, '\t');
                                             int pos = stoi(positions[1]);

                                             if (pos >= start_Co && pos <= end_Co)
                                             {
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
     functions function = functions();
     cout << "Writing to file\t: " << output_File << endl;
     cout << endl;

     int file_Count_Start = 0;
     int line_Num = 0;

     if (filesystem::exists(output_File) == 0)
     {
          function.createFile(output_File, "Coordinates\tPi\tS\tTotal_iEi\tFay_Wu_Normalized_H\tFay_Wu_Normalized_E");
     }
     else
     {
          int found = 0;

          fstream output_Check;
          output_Check.open(output_File, ios::in);
          if (output_Check.is_open())
          {
               string line_Check;
               getline(output_Check, line_Check); // skip first header line

               for (int file_Count = 0; file_Count < folder_Index.size(); file_Count++)
               {
                    string file_Path = folder_Index[file_Count].second;
                    fstream file;
                    file.open(file_Path, ios::in);
                    int line_Current = 0;

                    if (file.is_open())
                    {
                         string line;
                         getline(file, line); // skip first header line
                         while (getline(file, line))
                         {
                              line_Current++;
                              int VALID = function.get_Valid(line);
                              if (VALID != -1)
                              {
                                   getline(output_Check, line_Check);
                                   string trim = line_Check.substr(0, line_Check.find('\t'));

                                   vector<string> positions;
                                   function.split_getPos_ONLY(positions, line, '\t');
                                   string pos = positions[1] + ":" + to_string((stoi(positions[1]) + window_Size));

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

     int line_Current = 0;

     for (int file_Count = file_Count_Start; file_Count < folder_Index.size(); file_Count++)
     {
          string file_Path = folder_Index[file_Count].second;
          fstream file_Main;
          file_Main.open(file_Path, ios::in);

          if (file_Main.is_open())
          {
               string line_Main;
               getline(file_Main, line_Main); // skip first header line
               while (getline(file_Main, line_Main))
               {
                    if (line_Current < line_Num)
                    {
                         line_Current++;
                    }
                    else
                    {
                         // check VALID
                         int VALID = function.get_Valid(line_Main);
                         // cout << line_Main << endl;
                         if (VALID != -1)
                         {
                              int start_Co = VALID;
                              int end_Co = start_Co + window_Size;

                              cout << "Coordinates\t: Start: " << start_Co << " End: " << end_Co << endl;

                              float tot_pairwise_Differences = 0;

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
                                   // cout << files << endl;
                                   fstream file;
                                   file.open(files, ios::in);
                                   if (file.is_open())
                                   {
                                        string line;
                                        getline(file, line); // skip first header line
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

                              int num_segregrating_Sites;
                              string Fay_Wu_H, Fay_Wu_E;
                              float pi = 0.0;
                              int Total_iEi = 0;

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
     functions function = functions();

     int start_Value = stoi(folder_Index[0].first.substr(0, folder_Index[0].first.find('_')));
     int end_Value = stoi(folder_Index[folder_Index.size() - 1].first.substr(folder_Index[folder_Index.size() - 1].first.find('_') + 1));

     cout << start_Value << endl;
     cout << end_Value << endl;

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
          function.createFile(output_File, "Coordinates\tPi\tS\tTotal_iEi\tFay_Wu_Normalized_H\tFay_Wu_Normalized_E");
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
          vector<string> collect_Segregrating_sites;

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
          cout << endl;

          cout << "System is collecting segregrating site(s)" << endl;

          for (string files : file_List)
          {
               fstream file;
               file.open(files, ios::in);
               if (file.is_open())
               {
                    string line;
                    getline(file, line); // skip first header line
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

          int num_segregrating_Sites;
          string Fay_Wu_H, Fay_Wu_E;
          float pi = 0.0;
          int Total_iEi = 0;

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
     cout << "System is processing and filtering segregrating site(s)" << endl;
     int num_segregrating_Sites = total_Segregrating_sites.size();
     float theta_L = 0.00;
     string Seg_sites = "";
     int site_Index[num_segregrating_Sites + 1];
     site_Index[0] = 0;
     // cout << "copy string" << endl;
     for (size_t i = 0; i < num_segregrating_Sites; i++)
     {
          Seg_sites.append(total_Segregrating_sites[i]);
          site_Index[i + 1] = site_Index[i] + total_Segregrating_sites[i].size();
     }
     // cout << "copy string done" << endl;
     // cout << "Seg size : " << Seg_sites.size() << endl;
     char *full_Char;
     full_Char = (char *)malloc((Seg_sites.size() + 1) * sizeof(char));
     // cout << "copy char " << endl;
     strcpy(full_Char, Seg_sites.c_str());
     // cout << "copy char done" << endl;
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

     cudaMemcpy(cuda_full_Char, full_Char, (Seg_sites.size() + 1) * sizeof(char), cudaMemcpyHostToDevice);
     cudaMemcpy(cuda_site_Index, site_Index, (num_segregrating_Sites + 1) * sizeof(int), cudaMemcpyHostToDevice);
     // cout << "GPU" << endl;
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