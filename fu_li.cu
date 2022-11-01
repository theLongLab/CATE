#include "fu_li.cuh"
#include "functions.cuh"
#include "prometheus.cuh"

fu_li::fu_li(string gene_List, string input_Folder, string ouput_Path, int cuda_ID, string intermediate_Path, int ploidy)
{
    /**
     * * Constructor Function
     * NORMAL - GENE MODE constructor
     **/

    cout << "Initiating CUDA powered Fu and Li's D, D*, F and F* calculator" << endl
         << endl;

    set_Values(gene_List, input_Folder, ouput_Path, cuda_ID, intermediate_Path, ploidy);
}

fu_li::fu_li(string gene_List, string input_Folder, string ouput_Path, int cuda_ID, string intermediate_Path, int ploidy, string prometheus_Activate, string Multi_read, int number_of_genes, int CPU_cores, int SNPs_per_Run)
{
    /**
     * * Constructor Function
     * PROMETHEUS - GENE MODE constructor
     **/

    // PROMETHEUS Constructor
    cout << "Initiating CUDA powered Fu and Li's D, D*, F and F* calculator on PROMETHEUS" << endl
         << endl;

    set_Values(gene_List, input_Folder, ouput_Path, cuda_ID, intermediate_Path, ploidy);

    this->prometheus_Activate = "YES";
    this->CPU_cores = CPU_cores;
    this->SNPs_per_Run = SNPs_per_Run;
    transform(Multi_read.begin(), Multi_read.end(), Multi_read.begin(), ::toupper);
    this->Multi_read = Multi_read;
    this->number_of_genes = number_of_genes;
}

fu_li::fu_li(string calc_Mode, int window_Size, int step_Size, string input_Folder, string ouput_Path, int cuda_ID, int ploidy, string prometheus_Activate, string Multi_read, int number_of_genes, int CPU_cores, int SNPs_per_Run)
{
    /**
     * * Constructor Function
     * PROMETHEUS - WINDOW MODE constructor
     **/

    // PROMETHEUS Constructor WINDOW
    cout << "Initiating CUDA powered Fu and Li's D, D*, F and F* calculator on PROMETHEUS" << endl
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

fu_li::fu_li(string calc_Mode, int window_Size, int step_Size, string input_Folder, string ouput_Path, int cuda_ID, int ploidy)
{
    /**
     * * Constructor Function
     * NORMAL - WINDOW MODE constructor
     **/

    // NORMAL WINDOW CONSTRUCTOR

    cout << "Initiating CUDA powered Fu and Li's D, D*, F and F* calculator" << endl
         << endl;

    this->calc_Mode = "WINDOW";
    this->window_Size = window_Size;
    this->step_Size = step_Size;

    set_Values("", input_Folder, ouput_Path, cuda_ID, "", ploidy);
}

void fu_li::set_Values(string gene_List, string input_Folder, string ouput_Path, int cuda_ID, string intermediate_Path, int ploidy)
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

void fu_li::ingress()
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
        // float soft_Singl = 1 / N_float;
        // string SOFT_singleton_MAF = function.roundoff(soft_Singl, 4);
        // cout << "Singleton MAF\t: " << SOFT_singleton_MAF << endl;
        cout << endl;

        // calculate prerequisites
        /**
         * Pre-requisite values needed for determination of Fu an Li statistics.
         **/

        float an, vd, ud, vd_star, ud_star, uf, vf, uf_star, vf_star;
        calc_Pre(N, an, vd, ud, vd_star, ud_star, uf, vf, uf_star, vf_star);

        /**
         * @param test is used by Prometheus, to tell it which test is being processed.
         * * T  = Tajima
         * FU   = Fu and Li
         * FA   = Fay and Wu
         * N    = All 3 Neutrality tests
         **/
        string test = "FU";

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
                                 ".fl";

            /**
             * Ensures if PROMETHEUS is being activated.
             **/
            if (prometheus_Activate == "YES")
            {
                /**
                 * If Prometheus is being ACTIVATED then it is initialised accordingly.
                 **/
                prometheus pro_Fu_Li_Window = prometheus(output_File, window_Size, step_Size, folder_Index, Multi_read, tot_Blocks, tot_ThreadsperBlock, CPU_cores, SNPs_per_Run, number_of_genes, N, combinations, an, vd, ud, vd_star, ud_star, uf, vf, uf_star, vf_star);
                /**
                 * Ensures if it is NORMAL window or SLIDING window mode.
                 * If step_Size is = 0 then it is sliding window mode.
                 **/
                if (step_Size != 0)
                {
                    pro_Fu_Li_Window.process_Window(test);
                }
                else
                {
                    /**
                     * Initiates processing of Fu and Li on PROMETHEUS on sliding window mode.
                     **/
                    pro_Fu_Li_Window.process_C_sliding_Window(test);
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
                     * Initiates processing of Fu and Li on step wise window mode.
                     **/
                    window(output_File, an, vd, ud, vd_star, ud_star, uf, vf, uf_star, vf_star, N_float, combinations, folder_Index);
                }
                else
                {
                    /**
                     * Initiates processing of Fu and Li on sliding window mode.
                     **/
                    window_Sliding(output_File, an, vd, ud, vd_star, ud_star, uf, vf, uf_star, vf_star, N_float, combinations, folder_Index);
                }
            }
        }
        else
        {
            /**
             * * GENE (FILE) mode configuration:
             **/

            /**
             * Output file is created for the population in the output folder for FILE mode.
             * @param output_File stores the output file's location.
             * The file name is a combination of the the country, and gene file name.
             **/
            fstream gene_File;
            gene_File.open(gene_List, ios::in);
            cout << "Processing gene list:" << endl;
            string output_File = ouput_Path + "/" +
                                 country.substr(country.find_last_of("/") + 1, country.length()) + "_" +
                                 filesystem::path(gene_List).stem().string() +
                                 ".fl";
            /**
             * Log file created in the intermediate folder for the population.
             * @param intermediate_File stores the log file's location.
             * ! This helps with the resume function. Automatically resumes from the last completely processed gene in the event of a program crash.
             **/
            string intermediate_File = intermediate_Path + "/" +
                                       country.substr(country.find_last_of("/") + 1, country.length()) + "_" +
                                       filesystem::path(gene_List).stem().string() +
                                       ".log_fl";
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
                    function.createFile(output_File, "Gene_name\tCoordinates\tPi\tS\tne\tns\tD\tD_star\tF\tF_star");
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
                if (prometheus_Activate == "YES")
                {
                    cout << "Initializing Prometheus:" << endl
                         << endl;

                    /**
                     * If Prometheus is being ACTIVATED then it is initialised accordingly.
                     **/
                    prometheus pro_Fu_Li = prometheus(folder_Index, Multi_read, tot_Blocks, tot_ThreadsperBlock, CPU_cores, SNPs_per_Run, number_of_genes, N, combinations, an, vd, ud, vd_star, ud_star, uf, vf, uf_star, vf_star);

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
                            vector<string> write_Lines = pro_Fu_Li.collection_Engine(gene_Collect, test);
                            // print
                            cout << "System is writing Fu and Li results" << endl;
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
                            pro_Fu_Li.erase();
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
                        cout << "Prometheus batch initalized" << endl;
                        cout << "From: " << gene_Collect[0] << endl;
                        cout << "To  : " << gene_Collect[gene_Collect.size() - 1] << endl
                             << endl;

                        /**
                         * LAUNCH Prometheus to process the collected query batch.
                         * @param write_Lines vector collects the lines that should be written to the output file.
                         */
                        vector<string> write_Lines = pro_Fu_Li.collection_Engine(gene_Collect, test);
                        // print
                        cout << "System is writing Fu and Li results" << endl;
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
                    pro_Fu_Li.erase();
                    gene_Collect.clear();

                    // cout << endl;
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
                         * @param segregating_Sites Fu and Li requires the total number of segregating sites/ SNPS in the query region.
                         **/
                        float tot_pairwise_Differences = 0;
                        int segregating_Sites = 0;
                        int singletons_ns = 0;
                        int singletons_ne = 0;

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
                                        collect_Segregrating_sites.push_back(line);
                                        // string check_0 = country.substr(country.find_last_of("/") + 1, country.length()) + "_AF=0";
                                        // string seg_Check = "GO";
                                        // vector<string> info;
                                        // function.split(info, positions[7], ";");
                                        // for (string AF_check : info)
                                        // {
                                        //     if (AF_check == check_0)
                                        //     {
                                        //         seg_Check = "NO";
                                        //         break;
                                        //     }
                                        // }
                                        // if (seg_Check == "GO")
                                        // {
                                        //     string check_AF_country = country.substr(country.find_last_of("/") + 1, country.length()) + "_AF";
                                        //     float MAF_float = 0.0000;
                                        //     segregating_Sites = segregating_Sites + 1;
                                        //     for (string AF_check : info)
                                        //     {
                                        //         vector<string> split_info;
                                        //         function.split(split_info, AF_check, "=");
                                        //         if (split_info[0] == check_AF_country)
                                        //         {
                                        //             MAF_float = stof(split_info[1]);
                                        //             string MAF = split_info[1];
                                        //             string singleton_MAF;
                                        //             string check;
                                        //             if (MAF_float > 0.5)
                                        //             {
                                        //                 MAF_float = 1 - MAF_float;
                                        //                 //singleton_MAF = function.roundoff(soft_Singl, MAF.length() - 2);
                                        //                 string new_MAF = function.roundoff(MAF_float, MAF.length() - 2);
                                        //                 check = check_AF_country + "=" + new_MAF;
                                        //             }
                                        //             else
                                        //             {
                                        //                 //MAF = split_info[1];
                                        //                 check = AF_check;
                                        //             }

                                        //             singleton_MAF = function.roundoff(soft_Singl, MAF.length() - 2);
                                        //             string check_singleton = check_AF_country + "=" + singleton_MAF;
                                        //             if (check == check_singleton)
                                        //             {
                                        //                 //cout << "Singleton MAF check\t: " << check_singleton << "\t" << AF_check << endl;
                                        //                 singletons_ns = singletons_ns + 1;
                                        //                 singletons_ne = singletons_ne + outgroup_Singleton(info, positions);
                                        //             }
                                        //             break;
                                        //         }
                                        //     }
                                        //     //int pairwise_Differences = function.calc_Pairwise(line, N, this->tot_Blocks, this->tot_ThreadsperBlock);
                                        //     tot_pairwise_Differences = tot_pairwise_Differences + (MAF_float * (1 - MAF_float) * pow(N_float, 2));
                                        //     //break;
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

                        function.process_Seg_sites_fu_li(collect_Segregrating_sites, N_float, segregating_Sites, tot_pairwise_Differences, singletons_ne, singletons_ns, this->tot_Blocks, this->tot_ThreadsperBlock);

                        cout << endl;

                        // test
                        //  segregating_Sites = 18;
                        //  singletons_ne = 9;
                        //  singletons_ns = 10;
                        //"Gene_name\tCoordinates\tPi\tS\tne\tns\tD\tD_star\tF\tF_star"
                        cout << "Total segregating sites (S)\t: " << segregating_Sites << endl;

                        /**
                         * @param Fu_Li_D is used to store the Fu and Li D value for the query region.
                         * @param Fu_Li_D_star is used to store the Fu and Li D star value for the query region.
                         * @param Fu_Li_F is used to store the Fu and Li F value for the query region.
                         * @param Fu_Li_F_star is used to store the Fu and Li F star value for the query region.
                         **/
                        string Fu_Li_D;
                        string Fu_Li_D_star;
                        string Fu_Li_F;
                        string Fu_Li_F_star;
                        float pi = 0;

                        if (segregating_Sites != 0)
                        {
                            // test
                            // N_float = 24.0;
                            float D = (float)(segregating_Sites - (an * singletons_ne)) / sqrt(((ud * segregating_Sites) + (vd * (pow(segregating_Sites, 2)))));
                            float D_star = (float)(((N_float / (N_float - 1)) * segregating_Sites) - (an * singletons_ns)) / sqrt(((ud_star * segregating_Sites) + (vd_star * (pow(segregating_Sites, 2.0)))));

                            pi = (float)tot_pairwise_Differences / combinations;
                            // test
                            // pi = 3.16;

                            cout << "Average pairwise polymorphisms (pi)\t: " << pi << endl;
                            cout << "Total ns singletons\t: " << singletons_ns << endl;
                            cout << "Total ne singletons\t: " << singletons_ne << endl;
                            cout << endl;

                            float F = (float)(pi - singletons_ne) / sqrt(((uf * segregating_Sites) + (vf * (pow(segregating_Sites, 2)))));
                            float F_star = (float)(pi - (((N_float - 1) / N_float) * singletons_ns)) / sqrt(((uf_star * segregating_Sites) + (vf_star * (pow(segregating_Sites, 2.0)))));

                            Fu_Li_D = to_string(D);
                            Fu_Li_D_star = to_string(D_star);
                            Fu_Li_F = to_string(F);
                            Fu_Li_F_star = to_string(F_star);

                            cout << "Fu and Li's D\t: " << Fu_Li_D << endl;
                            cout << "Fu and Li's D*\t: " << Fu_Li_D_star << endl;
                            cout << "Fu and Li's F\t: " << Fu_Li_F << endl;
                            cout << "Fu and Li's F*\t: " << Fu_Li_F_star << endl;
                        }
                        else
                        {
                            cout << endl;
                            Fu_Li_D = "NA";
                            Fu_Li_D_star = "NA";
                            Fu_Li_F = "NA";
                            Fu_Li_F_star = "NA";
                            cout << "Fu and Li's D, D*, F and F*\t: "
                                 << "Not Available" << endl;
                        }

                        cout << endl;
                        //"Gene_name\tCoordinates\tPi\tS\tne\tns\tD\tD_star\tF\tF_star"
                        output << gene_Name << "\t"
                               << coordinates[0] << ":" << to_string(start_Co) << ":" << to_string(end_Co)
                               << "\t" << to_string(pi)
                               << "\t" << to_string(segregating_Sites)

                               << "\t" << to_string(singletons_ne)
                               << "\t" << to_string(singletons_ns)

                               << "\t" << Fu_Li_D
                               << "\t" << Fu_Li_D_star
                               << "\t" << Fu_Li_F
                               << "\t" << Fu_Li_F_star << "\n";

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

void fu_li::window_Sliding(string output_File, float an, float vd, float ud, float vd_star, float ud_star, float uf, float vf, float uf_star, float vf_star, float N_float, long int combinations, vector<pair<string, string>> &folder_Index)
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
        function.createFile(output_File, "Coordinates\tPi\tS\tne\tns\tD\tD_star\tF\tF_star");
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
                    /**
                     * Checks if the line being queried is a valid seg site.
                     * If so it is processed.
                     * @param VALID captures the position of the query site if it is valid, else it returns -1.
                     **/
                    // check VALID
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
                         * @param tot_pairwise_Differences Fu and Li also requires the tot_pairwise_Differences in the query region to determine the average number of pairwise differences in the region.
                         * @param segregating_Sites  Fu and Li requires the total number of segregating sites/ SNPS in the query region.
                         * @param singletons_ns accounts for all singleton mutations (mutations present only once in the population) in the region.
                         * @param singletons_ne accounts for all singleton mutations that are different from the AA.
                         **/
                        float tot_pairwise_Differences = 0;
                        int segregating_Sites = 0;
                        int singletons_ns = 0;
                        int singletons_ne = 0;

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

                        /**
                         * CALLs the function to calculate the actual segregating sites in the region where (MAF != 0) and the total pairwise differences.
                         **/
                        function.process_Seg_sites_fu_li(collect_Segregrating_sites, N_float, segregating_Sites, tot_pairwise_Differences, singletons_ne, singletons_ns, this->tot_Blocks, this->tot_ThreadsperBlock);
                        cout << endl;

                        cout << "Total segregating sites (S)\t: " << segregating_Sites << endl;

                        /**
                         * @param Fu_Li_D is used to store the Fu and Li D value for the query region.
                         * @param Fu_Li_D_star is used to store the Fu and Li D star value for the query region.
                         * @param Fu_Li_F is used to store the Fu and Li F value for the query region.
                         * @param Fu_Li_F_star is used to store the Fu and Li F star value for the query region.
                         **/
                        string Fu_Li_D;
                        string Fu_Li_D_star;
                        string Fu_Li_F;
                        string Fu_Li_F_star;
                        /**
                         * @param pi is used to store the average pairwise polymorphisms
                         **/
                        float pi = 0;

                        /**
                         * Ensure that there are segregating sites in the query region.
                         **/
                        if (segregating_Sites != 0)
                        {
                            float D = (float)(segregating_Sites - (an * singletons_ne)) / sqrt(((ud * segregating_Sites) + (vd * (pow(segregating_Sites, 2)))));
                            float D_star = (float)(((N_float / (N_float - 1)) * segregating_Sites) - (an * singletons_ns)) / sqrt(((ud_star * segregating_Sites) + (vd_star * (pow(segregating_Sites, 2.0)))));

                            pi = (float)tot_pairwise_Differences / combinations;

                            cout << "Average pairwise polymorphisms (pi)\t: " << pi << endl;
                            cout << "Total ns singletons\t: " << singletons_ns << endl;
                            cout << "Total ne singletons\t: " << singletons_ne << endl;
                            cout << endl;

                            float F = (float)(pi - singletons_ne) / sqrt(((uf * segregating_Sites) + (vf * (pow(segregating_Sites, 2)))));
                            float F_star = (float)(pi - (((N_float - 1) / N_float) * singletons_ns)) / sqrt(((uf_star * segregating_Sites) + (vf_star * (pow(segregating_Sites, 2.0)))));

                            Fu_Li_D = to_string(D);
                            Fu_Li_D_star = to_string(D_star);
                            Fu_Li_F = to_string(F);
                            Fu_Li_F_star = to_string(F_star);

                            cout << "Fu and Li's D\t: " << Fu_Li_D << endl;
                            cout << "Fu and Li's D*\t: " << Fu_Li_D_star << endl;
                            cout << "Fu and Li's F\t: " << Fu_Li_F << endl;
                            cout << "Fu and Li's F*\t: " << Fu_Li_F_star << endl;
                        }
                        else
                        {
                            cout << endl;
                            Fu_Li_D = "NA";
                            Fu_Li_D_star = "NA";
                            Fu_Li_F = "NA";
                            Fu_Li_F_star = "NA";
                            cout << "Fu and Li's D, D*, F and F*\t: "
                                 << "Not Available" << endl;
                        }

                        cout << endl;

                        output << to_string(start_Co) << ":" << to_string(end_Co)
                               << "\t" << to_string(pi)
                               << "\t" << to_string(segregating_Sites)

                               << "\t" << to_string(singletons_ne)
                               << "\t" << to_string(singletons_ns)

                               << "\t" << Fu_Li_D
                               << "\t" << Fu_Li_D_star
                               << "\t" << Fu_Li_F
                               << "\t" << Fu_Li_F_star << "\n";

                        output.flush();
                    }
                }
            }
            file_Main.close();
        }
    }

    output.close();
}

void fu_li::window(string output_File, float an, float vd, float ud, float vd_star, float ud_star, float uf, float vf, float uf_star, float vf_star, float N_float, long int combinations, vector<pair<string, string>> &folder_Index)
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
     * WE cycle through till we fall into the range of SNPs available in our VCFs.
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
        function.createFile(output_File, "Coordinates\tPi\tS\tne\tns\tD\tD_star\tF\tF_star");
    }
    else
    {
        /**
         * If the output file is already present then the resume process will initiated.
         * This is a unintelligent resume. Essentially it matches the each read line written with the lines read from the gene file.
         * The break will occur as soon as their is a mismatch.
         * To counter any errors it is advised to have a new gene file name or a new intermediate folder per new run.
         * @param caught acts as a boolean variable. caught = 0 if the lines need to be skipped and will equal 1 when the resume position is found.
         * ! This helps with the resume function. Automatically resumes from the last completely processed gene in the event of a program crash.
         **/
        // RESUME FUNCTION
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

            // skip header
            /**
             * First line is skipped cause it is a header line containing column names.
             **/
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
         * @param tot_pairwise_Differences Fu and Li also requires the tot_pairwise_Differences in the query region to determine the average number of pairwise differences in the region.
         * @param segregating_Sites  Fu and Li requires the total number of segregating sites/ SNPS in the query region.
         * @param singletons_ns accounts for all singleton mutations (mutations present only once in the population) in the region.
         * @param singletons_ne accounts for all singleton mutations that are different from the AA.
         **/
        float tot_pairwise_Differences = 0;
        int segregating_Sites = 0;
        int singletons_ns = 0;
        int singletons_ne = 0;

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

        /**
         * The SNPs (Segregating sites) that fall within the query region are collected from the VCF's.
         * @param collect_Segregrating_sites vector stores the collected SNPs.
         **/
        cout << "System is collecting segregating site(s)" << endl;
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

        /**
         * CALLs the function to calculate the actual segregating sites in the region where (MAF != 0) and the total pairwise differences.
         **/
        function.process_Seg_sites_fu_li(collect_Segregrating_sites, N_float, segregating_Sites, tot_pairwise_Differences, singletons_ne, singletons_ns, this->tot_Blocks, this->tot_ThreadsperBlock);
        cout << endl;

        cout << "Total segregating sites (S)\t: " << segregating_Sites << endl;

        /**
         * @param Fu_Li_D is used to store the Fu and Li D value for the query region.
         * @param Fu_Li_D_star is used to store the Fu and Li D star value for the query region.
         * @param Fu_Li_F is used to store the Fu and Li F value for the query region.
         * @param Fu_Li_F_star is used to store the Fu and Li F star value for the query region.
         **/
        string Fu_Li_D;
        string Fu_Li_D_star;
        string Fu_Li_F;
        string Fu_Li_F_star;
        /**
         * @param pi is used to store the average pairwise polymorphisms
         **/
        float pi = 0;

        if (segregating_Sites != 0)
        {
            float D = (float)(segregating_Sites - (an * singletons_ne)) / sqrt(((ud * segregating_Sites) + (vd * (pow(segregating_Sites, 2)))));
            float D_star = (float)(((N_float / (N_float - 1)) * segregating_Sites) - (an * singletons_ns)) / sqrt(((ud_star * segregating_Sites) + (vd_star * (pow(segregating_Sites, 2.0)))));

            pi = (float)tot_pairwise_Differences / combinations;

            cout << "Average pairwise polymorphisms (pi)\t: " << pi << endl;
            cout << "Total ns singletons\t: " << singletons_ns << endl;
            cout << "Total ne singletons\t: " << singletons_ne << endl;
            cout << endl;

            float F = (float)(pi - singletons_ne) / sqrt(((uf * segregating_Sites) + (vf * (pow(segregating_Sites, 2)))));
            float F_star = (float)(pi - (((N_float - 1) / N_float) * singletons_ns)) / sqrt(((uf_star * segregating_Sites) + (vf_star * (pow(segregating_Sites, 2.0)))));

            Fu_Li_D = to_string(D);
            Fu_Li_D_star = to_string(D_star);
            Fu_Li_F = to_string(F);
            Fu_Li_F_star = to_string(F_star);

            cout << "Fu and Li's D\t: " << Fu_Li_D << endl;
            cout << "Fu and Li's D*\t: " << Fu_Li_D_star << endl;
            cout << "Fu and Li's F\t: " << Fu_Li_F << endl;
            cout << "Fu and Li's F*\t: " << Fu_Li_F_star << endl;
        }
        else
        {
            cout << endl;
            Fu_Li_D = "NA";
            Fu_Li_D_star = "NA";
            Fu_Li_F = "NA";
            Fu_Li_F_star = "NA";
            cout << "Fu and Li's D, D*, F and F*\t: "
                 << "Not Available" << endl;
        }

        cout << endl;

        output << to_string(start_Co) << ":" << to_string(end_Co)
               << "\t" << to_string(pi)
               << "\t" << to_string(segregating_Sites)

               << "\t" << to_string(singletons_ne)
               << "\t" << to_string(singletons_ns)

               << "\t" << Fu_Li_D
               << "\t" << Fu_Li_D_star
               << "\t" << Fu_Li_F
               << "\t" << Fu_Li_F_star << "\n";

        output.flush();

        start_Co = start_Co + step_Size;
        end_Co = start_Co + window_Size;
    }

    output.close();
}

// int fu_li::outgroup_Singleton(vector<string> &info, vector<string> &positions)
// {
//     functions function = functions();
//     string present = "false";
//     string MA;
//     string AA;
//     int return_Val = 0;
//     vector<string> line_partial;

//     for (string info_Data : info)
//     {
//         vector<string> AA_check;
//         function.split(AA_check, info_Data, '=');
//         if (AA_check[0] == "AA")
//         {
//             // cout << "AA: " << AA_check[1].at(0) << endl;
//             if ((toupper(AA_check[1].at(0)) == 'A') || (toupper(AA_check[1].at(0)) == 'T') || (toupper(AA_check[1].at(0)) == 'G') || (toupper(AA_check[1].at(0)) == 'C'))
//             {
//                 // cout << "AA Caught: " << AA_check[1].at(0) << endl;
//                 AA = AA_check[1].at(0);
//                 present = "true";
//             }
//             break;
//         }
//     }

//     if (present == "false")
//     {
//         AA = positions[3];
//     }

//     int REF_0 = 0;
//     int ALT_1 = 0;

//     vector<string> sample_1;
//     vector<string> sample_2;

//     function.split(sample_1, positions[9], '|');
//     function.split(sample_2, positions[10], '|');

//     for (string check : sample_1)
//     {
//         if (check == "0")
//         {
//             REF_0 = REF_0 + 1;
//         }
//         else
//         {
//             ALT_1 = ALT_1 + 1;
//         }
//     }

//     for (string check : sample_2)
//     {
//         if (check == "0")
//         {
//             REF_0 = REF_0 + 1;
//         }
//         else
//         {
//             ALT_1 = ALT_1 + 1;
//         }
//     }

//     if (REF_0 > ALT_1)
//     {
//         MA = positions[4];
//     }
//     else
//     {
//         MA = positions[3];
//     }

//     if (toupper(AA.at(0)) == toupper(MA.at(0)))
//     {
//         return_Val = 0;
//     }
//     else
//     {
//         return_Val = 1;
//     }

//     return return_Val;
// }

__global__ void fuli_Calculation(int N, float *a1_CUDA, float *a2_CUDA)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < N)
    {
        a1_CUDA[tid] = (float)1 / (tid + 1);
        a2_CUDA[tid] = (float)1 / ((tid + 1) * (tid + 1));
        tid += blockDim.x * gridDim.x;
    }
}

void fu_li::calc_Pre(int N_tot, float &an, float &vd, float &ud, float &vd_star, float &ud_star, float &uf, float &vf, float &uf_star, float &vf_star)
{
    /**
     * Calculates the prerequisite values required for Fu and Li
     **/

    // functions function = functions();
    // test
    // N_tot = 24;
    int N = N_tot - 1;
    float *a1_CUDA, *a2_CUDA;
    float *a1_partial, *a2_partial, bn;

    a1_partial = (float *)malloc(N * sizeof(float));
    a2_partial = (float *)malloc(N * sizeof(float));

    cudaMallocManaged(&a1_CUDA, N * sizeof(int));
    cudaMallocManaged(&a2_CUDA, N * sizeof(int));

    fuli_Calculation<<<tot_Blocks, tot_ThreadsperBlock>>>(N, a1_CUDA, a2_CUDA);
    cudaDeviceSynchronize();

    cudaMemcpy(a1_partial, a1_CUDA, N * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(a2_partial, a2_CUDA, N * sizeof(float), cudaMemcpyDeviceToHost);

    cudaFree(a1_CUDA);
    cudaFree(a2_CUDA);

    // an = function.add(N, a1_partial);
    // bn = function.add(N, a2_partial);

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

    float an_plus_1 = (float)an + (1 / N_float_tot);
    // cout << "an+1: " << an_plus_1 << endl;
    float cn = (float)2 * (((N_float_tot * an) - (2 * N_float)) / (N_float * (N_float_tot - 2)));
    float an_square = (float)pow(an, 2.0);

    // float test = (N_float_tot_float + 1) / N_float_float;
    // cout << "test: " << test << endl;

    /**
       * The D, D* and F statistics are calculated based on the original paper by Fu et al (1993).
              The F* statistic's vf* and uf* are calculated based on the corrected equations in Simonsen et al (1995).
    **/

    vd = (float)1 + ((an_square / (bn + an_square)) * (cn - (N_float_tot + 1) / N_float));
    ud = (float)an - 1 - vd;

    float dn = (float)cn + ((N_float_tot - 2) / pow(N_float, 2.0)) + ((2 / N_float) * (1.5 - (((2 * an_plus_1) - 3) / (N_float_tot - 2)) - (1 / N_float_tot)));
    // cout << "dn: " << dn << endl;
    vd_star = (float)((pow((N_float_tot / N_float), 2.0) * bn) + (pow(an, 2.0) * dn) - (2 * ((N_float_tot * an * (an + 1)) / (pow(N_float, 2.0))))) / (pow(an, 2.0) + bn);
    ud_star = (float)N_float_tot / N_float * (an - (N_float_tot / N_float)) - vd_star;

    vf = (float)(cn + ((2 * (pow(N_float_tot, 2) + N_float_tot + 3)) / (9 * N_float_tot * N_float)) - (2 / N_float)) / (pow(an, 2) + bn);
    uf = (float)((1 + ((N_float_tot + 1) / (3 * N_float)) - ((4 * ((N_float_tot + 1) / pow(N_float, 2))) * (an_plus_1 - ((2 * N_float_tot) / (N_float_tot + 1))))) / an) - vf;

    // Original papers incorrect equations due to paper's typos
    // vf_star = (float)(dn + ((2 * (pow(N_float_tot, 2) + N_float_tot + 3)) / (9 * N_float_tot * N_float)) - ((2 * (1 / N_float)) * ((4 * bn) - 6 + (8 / N_float_tot)))) / (pow(an, 2) + bn);
    // uf_star = (float)(((N_float_tot / N_float) + ((N_float_tot + 1) / (3 * N_float)) - (2 * (2 / (N_float_tot * N_float))) + ((2 * ((N_float_tot + 1) / pow(N_float, 2))) * (an_plus_1 - ((2 * N_float_tot) / (N_float_tot + 1))))) / an) - vf_star;

    // corrected equations Simonsen et al 1995
    vf_star = (float)((((2 * pow(N_float_tot, 3.0)) + (110 * pow(N_float_tot, 2.0)) - (225 * N_float_tot) + 153) / (9 * pow(N_float_tot, 2.0) * N_float)) + ((2 * N_float * an) / (pow(N_float_tot, 2.0))) - ((8 * bn) / N_float_tot)) / (pow(an, 2.0) + bn);
    uf_star = (float)((((4 * pow(N_float_tot, 2.0)) + (19 * N_float_tot) + 3 - (12 * (N_float_tot + 1) * an_plus_1)) / (3 * N_float_tot * (N_float))) / an) - vf_star;

    cout << "an: " << an << "\t"
         << "bn: " << bn << "\t"
         << "cn: " << cn << endl;

    cout << endl;

    cout << "ud\t: " << ud << "\t"
         << "vd\t: " << vd << endl;

    cout << "ud*\t: " << ud_star << "\t"
         << "vd*\t: " << vd_star << endl;

    cout << endl;

    cout << "an+1\t: " << an_plus_1 << "\t"
         << "dn\t: " << dn << endl;

    cout << endl;

    cout << "uf\t: " << uf << "\t"
         << "vf\t: " << vf << endl;

    cout << "uf*\t: " << uf_star << "\t"
         << "vf*\t: " << vf_star << endl;

    cout << endl;
}