#include "tajima.cuh"
#include "functions.cuh"
#include "prometheus.cuh"

tajima::tajima(string gene_List, string input_Folder, string ouput_Path, int cuda_ID, string intermediate_Path, int ploidy)
{
    /**
     * * Constructor Function
     * NORMAL - GENE MODE constructor
     **/

    cout << "Initiating CUDA powered Tajima's D calculator" << endl
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

tajima::tajima(string gene_List, string input_Folder, string ouput_Path, int cuda_ID, string intermediate_Path, int ploidy, string prometheus_Activate, string Multi_read, int number_of_genes, int CPU_cores, int SNPs_per_Run)
{
    /**
     * * Constructor Function
     * PROMETHEUS - GENE MODE constructor
     **/

    // PROMETHEUS Constructor gene file
    cout << "Initiating CUDA powered Tajima's D calculator on PROMETHEUS" << endl
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

    this->prometheus_Activate = "YES";
    this->CPU_cores = CPU_cores;
    this->SNPs_per_Run = SNPs_per_Run;
    /**
     * Multi_read is converted to uppercase to create uniformity prevent any user error.
     **/
    transform(Multi_read.begin(), Multi_read.end(), Multi_read.begin(), ::toupper);
    this->Multi_read = Multi_read;
    this->number_of_genes = number_of_genes;
}

tajima::tajima(string calc_Mode, int window_Size, int step_Size, string input_Folder, string ouput_Path, int cuda_ID, int ploidy, string prometheus_Activate, string Multi_read, int number_of_genes, int CPU_cores, int SNPs_per_Run)
{
    /**
     * * Constructor Function
     * PROMETHEUS - WINDOW MODE constructor
     **/

    // PROMETHEUS WINDOW MODE CONSTRUCTOR

    cout << "Initiating CUDA powered Tajima's D calculator on PROMETHEUS" << endl
         << endl;

    this->calc_Mode = "WINDOW";
    this->window_Size = window_Size;
    this->step_Size = step_Size;

    /**
     * gene_List and intermediate_Path variables are kept blank.
     * Because in WINDOW mode there is no requirement for a gene file and the resume function works off the indexed VCF files.
     **/
    set_Values("", input_Folder, ouput_Path, cuda_ID, "", ploidy);

    this->prometheus_Activate = "YES";
    this->CPU_cores = CPU_cores;
    this->SNPs_per_Run = SNPs_per_Run;
    /**
     * Multi_read is converted to uppercase to create uniformity prevent any user error.
     **/
    transform(Multi_read.begin(), Multi_read.end(), Multi_read.begin(), ::toupper);
    this->Multi_read = Multi_read;
    this->number_of_genes = number_of_genes;
}

tajima::tajima(string calc_Mode, int window_Size, int step_Size, string input_Folder, string ouput_Path, int cuda_ID, int ploidy)
{
    /**
     * * Constructor Function
     * NORMAL - WINDOW MODE constructor
     **/

    // NORMAL WINDOW CONSTRUCTOR

    cout << "Initiating CUDA powered Tajima's D calculator" << endl
         << endl;

    this->calc_Mode = "WINDOW";

    this->window_Size = window_Size;
    this->step_Size = step_Size;

    /**
     * gene_List and intermediate_Path variables are kept blank.
     * Because in WINDOW mode there is no requirement for a gene file and the resume function works off the indexed VCF files.
     **/
    set_Values("", input_Folder, ouput_Path, cuda_ID, "", ploidy);
}

void tajima::set_Values(string gene_List, string input_Folder, string ouput_Path, int cuda_ID, string intermediate_Path, int ploidy)
{
    /**
     * This function is used in conjunction with the constructor to set the common private variables.
     * Notifies the user if it is WINDOW mode or GENE (FILE) mode.
     * If WINDOW user is also notified if it is sliding window or normal step wise window mode.
     * Here the first call to the selected CUDA device occurs.
     **/

    if (this->calc_Mode == "WINDOW")
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
    else
    {
        cout << "Calculation mode: FILE" << endl;
        this->gene_List = gene_List;
        cout << "Gene list file path\t: " << gene_List << endl
             << endl;
    }

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

void tajima::ingress()
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
    cout << countries.size() << " population(s) were found: ";
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

        // for (auto file : folder_Index)
        // {
        //     cout << file.first << "\t" << file.second << endl;
        // }

        cout << "Completed indexing folder\t: " << country << endl;
        // string check_AF_country = country.substr(country.find_last_of("/") + 1, country.length()) + "_AF";
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
        cout << "Number of sequences in " << country.substr(country.find_last_of("/") + 1, country.length()) << " population [ " << samples << " x " << ploidy << " ] (N)\t: " << N << endl;

        /**
         * @param combinations defines number of total pairwise combinations being present.
         **/
        long int combinations = function.combos_N(N);
        cout << "Pairwise combinations\t: " << combinations << endl;
        cout << endl;

        /**
         * Pre-requisite values needed for determination of Tajima's D.
         * @param a1 is a weight estimator.
         * @param e1 is a weight estimator.
         * @param e2 is a weight estimator.
         * * reference: https://ocw.mit.edu/courses/hst-508-quantitative-genomics-fall-2005/0900020a633e85338e9510495c2e01a6_tajimad1.pdf
         **/
        float a1, e1, e2;
        calc_Pre(N, a1, e1, e2);

        /**
         * @param test is used by Prometheus, to tell it which test is being processed.
         * * T  = Tajima
         * FU   = Fu and Li
         * FA   = Fay and Wu
         * N    = All 3 Neutrality tests
         **/
        string test = "T";

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
                                 ".td";

            /**
             * Ensures if PROMETHEUS is being activated.
             **/
            if (prometheus_Activate == "YES")
            {
                /**
                 * If Prometheus is being ACTIVATED then it is initialised accordingly.
                 **/
                prometheus pro_Tajima_Window = prometheus(output_File, window_Size, step_Size, folder_Index, Multi_read, tot_Blocks, tot_ThreadsperBlock, combinations, a1, e1, e2, N, CPU_cores, SNPs_per_Run, number_of_genes);
                /**
                 * Ensures if it is NORMAL window or SLIDING window mode.
                 * If step_Size is = 0 then it is sliding window mode.
                 **/
                if (step_Size != 0)
                {
                    /**
                     * Initiates processing of Tajima on PROMETHEUS on step wise window mode.
                     **/
                    pro_Tajima_Window.process_Window(test);
                }
                else
                {
                    /**
                     * Initiates processing of Tajima on PROMETHEUS on sliding window mode.
                     **/
                    pro_Tajima_Window.process_C_sliding_Window(test);
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
                    window(output_File, a1, e1, e2, N, combinations, folder_Index);
                }
                else
                {
                    /**
                     * Initiates processing of Tajima on sliding window mode.
                     **/
                    window_Sliding(output_File, a1, e1, e2, N, combinations, folder_Index);
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
            string output_File = ouput_Path + "/" +
                                 country.substr(country.find_last_of("/") + 1, country.length()) + "_" +
                                 filesystem::path(gene_List).stem().string() +
                                 ".td";
            /**
             * Log file created in the intermediate folder for the population.
             * @param intermediate_File stores the log file's location.
             * ! This helps with the resume function. Automatically resumes from the last completely processed gene in the event of a program crash.
             **/
            string intermediate_File = intermediate_Path + "/" +
                                       country.substr(country.find_last_of("/") + 1, country.length()) + "_" +
                                       filesystem::path(gene_List).stem().string() +
                                       ".log_td";

            /**
             * Initiate the reading of the gene file.
             **/
            fstream gene_File;
            gene_File.open(gene_List, ios::in);
            cout << "Processing gene list:" << endl;

            cout << endl;

            cout << "Writing to file\t: " << output_File << endl;
            cout << endl;

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
                    function.createFile(output_File, "Gene_name\tCoordinates\tPi\tS\tTajimas_D");
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
                // cout << "CHECK: " << gene_Combo << endl;

                fstream output;
                fstream intermediate;
                output.open(output_File, ios::app);
                intermediate.open(intermediate_File, ios::app);

                /**
                 * Ensures if PROMETHEUS is being activated.
                 **/
                // ADD Prometheus HERE
                if (prometheus_Activate == "YES")
                {
                    cout << "Initializing Prometheus:" << endl
                         << endl;

                    // cout << "Processing on " << this->CPU_cores << " CPU cores" << endl;
                    // cout << "Processing " << this->number_of_genes << " genes at a time" << endl;
                    // cout << "Processing " << this->SNPs_per_Run << " SNPs at a time" << endl;
                    // if (this->Multi_read == "YES")
                    // {
                    //     cout << "Multi read: Available" << endl;
                    // }
                    // else
                    // {
                    //     cout << "Multi read: Unavailable" << endl;
                    // }
                    // cout << endl;

                    /**
                     * If Prometheus is being ACTIVATED then it is initialised accordingly.
                     **/
                    prometheus pro_Tajima = prometheus(folder_Index, Multi_read, this->tot_Blocks, this->tot_ThreadsperBlock, combinations, a1, e1, e2, N, CPU_cores, SNPs_per_Run, number_of_genes);

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
                            cout << "Prometheus batch initalized" << endl;
                            cout << "From: " << gene_Collect[0] << endl;
                            cout << "To  : " << gene_Collect[gene_Collect.size() - 1] << endl
                                 << endl;

                            /**
                             * LAUNCH Prometheus to process the collected query batch.
                             * @param write_Lines vector collects the lines that should be written to the output file.
                             */
                            // launch prometheus
                            vector<string> write_Lines = pro_Tajima.collection_Engine(gene_Collect, test);
                            // print
                            cout << "System is writing Tajima's D results" << endl;
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
                            pro_Tajima.erase();
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
                        cout << "Prometheus batch initalized" << endl;
                        cout << "From: " << gene_Collect[0] << endl;
                        cout << "To  : " << gene_Collect[gene_Collect.size() - 1] << endl
                             << endl;

                        /**
                         * LAUNCH Prometheus to process the collected query batch.
                         * @param write_Lines vector collects the lines that should be written to the output file.
                         */
                        vector<string> write_Lines = pro_Tajima.collection_Engine(gene_Collect, test);
                        // print

                        cout << "System is writing Tajima's D results" << endl;
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
                    pro_Tajima.erase();
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

                        // cout << gene_Combo << endl;
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
                        // int tot_pairwise_Differences_TEST = 0;
                        int segregating_Sites = 0;

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
                                        /**
                                         * If the SNP is between the range of the query region, it is collected.
                                         * Information from the SNP is extracted via the GPU.
                                         **/
                                        collect_Segregrating_sites.push_back(line);
                                        //     //cout << pos << endl;
                                        //     string check_0 = country.substr(country.find_last_of("/") + 1, country.length()) + "_AF=0";
                                        //     string GO = "GO";
                                        //     vector<string> info;
                                        //     split(info, positions[7], ";");
                                        //     for (string AF_check : info)
                                        //     {
                                        //         if (AF_check == check_0)
                                        //         {
                                        //             // cout << pos << endl;
                                        //             GO = "NO";
                                        //             break;
                                        //         }
                                        //     }
                                        //     if (GO == "GO")
                                        //     {
                                        //         //string check_AF_country = country.substr(country.find_last_of("/") + 1, country.length()) + "_AF";
                                        //         float MAF = 0.0000;
                                        //         for (string AF_check : info)
                                        //         {
                                        //             vector<string> split_info;
                                        //             split(split_info, AF_check, "=");
                                        //             if (split_info[0] == check_AF_country)
                                        //             {
                                        //                 MAF = stof(split_info[1]);
                                        //                 if (MAF > 0.5)
                                        //                 {
                                        //                     MAF = 1 - MAF;
                                        //                 }
                                        //                 //cout << split_info[0] << "\t: " << MAF << endl;
                                        //                 break;
                                        //             }
                                        //         }
                                        //         tot_pairwise_Differences = tot_pairwise_Differences + (MAF * (1 - MAF) * pow(N, 2));
                                        //         //cout << "pairwise differences: \t" << (MAF * (1 - MAF) * pow(N, 2)) << endl;
                                        //         //int pairwise_Differences = calc_Pairwise(line, N);
                                        //         //cout << pairwise_Differences << endl;
                                        //         //tot_pairwise_Differences = tot_pairwise_Differences + pairwise_Differences;
                                        //         segregating_Sites = segregating_Sites + 1;
                                        //         //break;
                                        //     }
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
                        function.process_Seg_sites_tajima(collect_Segregrating_sites, N, segregating_Sites, tot_pairwise_Differences, this->tot_Blocks, this->tot_ThreadsperBlock);

                        cout << endl;
                        // cout << "totDif " << tot_pairwise_Differences << endl;

                        /**
                         * @param pi is used to store the average pairwise polymorphisms
                         * @param D is used to store the Tajima's D value for the query region.
                         **/
                        float pi = 0;
                        // float pi_Test = 0;
                        float D = 0;
                        string Tajima_D;
                        cout << "Total segregating sites (S)\t: " << segregating_Sites << endl;

                        /**
                         * Ensure that there are segregating sites in the query region.
                         **/
                        if (segregating_Sites != 0)
                        {
                            // pi_Test = (float)tot_pairwise_Differences_TEST / combinations;
                            pi = (float)tot_pairwise_Differences / combinations;
                            // cout << "TEST Average pairwise polymorphisms (pi)\t: " << pi_Test << endl;
                            cout << "Average pairwise polymorphisms (pi)\t: " << pi << endl;
                            D = (float)(pi - (segregating_Sites / a1)) / sqrt(((e1 * segregating_Sites) + (e2 * segregating_Sites * (segregating_Sites - 1))));
                            cout << endl;
                            cout << "Tajima's D\t: " << D << endl;
                            Tajima_D = to_string(D);
                        }
                        else
                        {
                            cout << endl;
                            cout << "Tajima's D\t: "
                                 << "Not Available" << endl;
                            Tajima_D = "NA";
                        }

                        cout << endl;

                        //"Gene_name\tCoordinates\tPi\tS\tTajimas_D"
                        output << gene_Name << "\t"
                               << coordinates[0] << ":" << to_string(start_Co) << ":" << to_string(end_Co)
                               << "\t" << to_string(pi)
                               << "\t" << to_string(segregating_Sites)

                               << "\t" << Tajima_D << "\n";

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

void tajima::window_Sliding(string output_File, float a1, float e1, float e2, int N, long int combinations, vector<pair<string, string>> &folder_Index)
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
        function.createFile(output_File, "Coordinates\tPi\tS\tTajimas_D");
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
                         * @param tot_pairwise_Differences Tajima's D also requires the tot_pairwise_Differences in the query region to determine the average number of pairwise differences in the region.
                         * @param segregating_Sites Tajima's D requires the total number of segregating sites/ SNPS in the query region.
                         **/
                        float tot_pairwise_Differences = 0;
                        int segregating_Sites = 0;

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
                        function.process_Seg_sites_tajima(collect_Segregrating_sites, N, segregating_Sites, tot_pairwise_Differences, this->tot_Blocks, this->tot_ThreadsperBlock);

                        cout << endl;

                        /**
                         * @param pi is used to store the average pairwise polymorphisms
                         * @param D is used to store the Tajima's D value for the query region.
                         **/
                        float pi = 0;
                        float D = 0;

                        string Tajima_D;
                        cout << "Total segregating sites (S)\t: " << segregating_Sites << endl;

                        /**
                         * Ensure that there are segregating sites in the query region.
                         **/
                        if (segregating_Sites != 0)
                        {
                            pi = (float)tot_pairwise_Differences / combinations;
                            cout << "Average pairwise polymorphisms (pi)\t: " << pi << endl;
                            D = (float)(pi - (segregating_Sites / a1)) / sqrt(((e1 * segregating_Sites) + (e2 * segregating_Sites * (segregating_Sites - 1))));
                            cout << endl;
                            cout << "Tajima's D\t: " << D << endl;
                            Tajima_D = to_string(D);
                        }
                        else
                        {
                            cout << endl;
                            cout << "Tajima's D\t: "
                                 << "Not Available" << endl;
                            Tajima_D = "NA";
                        }

                        cout << endl;

                        output << to_string(start_Co) << ":" << to_string(end_Co)
                               << "\t" << to_string(pi)
                               << "\t" << to_string(segregating_Sites)

                               << "\t" << Tajima_D << "\n";

                        output.flush();
                    }

                    // REMOVE after test
                    // break;
                }
            }
            file_Main.close();
        }
        // REMOVE after test
        // break;
    }

    output.close();
}

void tajima::window(string output_File, float a1, float e1, float e2, int N, long int combinations, vector<pair<string, string>> &folder_Index)
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
        function.createFile(output_File, "Coordinates\tPi\tS\tTajimas_D");
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

        /**
         * Open the output file to be begin finding the resume point.
         **/
        // skipper
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

    // start_Co = 26820000;
    // end_Co = 26830000;

    while (start_Co <= end_Value)
    {
        cout << "Coordinates\t: Start: " << start_Co << " End: " << end_Co << endl;

        /**
         * @param tot_pairwise_Differences Tajima's D also requires the tot_pairwise_Differences in the query region to determine the average number of pairwise differences in the region.
         * @param segregating_Sites Tajima's D requires the total number of segregating sites/ SNPS in the query region.
         **/
        float tot_pairwise_Differences = 0;
        int segregating_Sites = 0;

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
        function.process_Seg_sites_tajima(collect_Segregrating_sites, N, segregating_Sites, tot_pairwise_Differences, this->tot_Blocks, this->tot_ThreadsperBlock);

        cout << endl;

        /**
         * @param pi is used to store the average pairwise polymorphisms
         * @param D is used to store the Tajima's D value for the query region.
         **/
        float pi = 0;
        float D = 0;

        string Tajima_D;
        cout << "Total segregating sites (S)\t: " << segregating_Sites << endl;

        if (segregating_Sites != 0)
        {
            pi = (float)tot_pairwise_Differences / combinations;
            cout << "Average pairwise polymorphisms (pi)\t: " << pi << endl;
            D = (float)(pi - (segregating_Sites / a1)) / sqrt(((e1 * segregating_Sites) + (e2 * segregating_Sites * (segregating_Sites - 1))));
            cout << endl;
            cout << "Tajima's D\t: " << D << endl;
            Tajima_D = to_string(D);
        }
        else
        {
            cout << endl;
            cout << "Tajima's D\t: "
                 << "Not Available" << endl;
            Tajima_D = "NA";
        }

        cout << endl;

        output << to_string(start_Co) << ":" << to_string(end_Co)
               << "\t" << to_string(pi)
               << "\t" << to_string(segregating_Sites)

               << "\t" << Tajima_D << "\n";

        output.flush();

        start_Co = start_Co + step_Size;
        end_Co = start_Co + window_Size;
    }

    output.close();
}

// void tajima::createFile(string path)
// {
//     fstream file;
//     file.open(path, ios::out);
//     file.close();
// }

// void tajima::createFile(string path, string text)
// {
//     fstream file;
//     file.open(path, ios::out);
//     file << text;
//     file << "\n";
//     file.close();
// }

__global__ void pairwise_Cuda(int N, int *SNP, int *differences)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < N)
    {
        int tot = 0;
        for (int i = tid + 1; i < N; i++)
        {
            // printf("snp1 is %d and snp %d is %d \n",SNP[tid],i,SNP[i]);
            if (SNP[tid] != SNP[i])
            {
                tot = tot + 1;
            }
        }
        // printf("tid is %d and diff %d \n",tid,tot);
        differences[tid] = tot;
        tid += blockDim.x * gridDim.x;
    }
}

// int tajima::calc_Pairwise(string &line, int N)
// {
//     int pairwise_Differences = 0;

//     // int *line_temp = (int *)malloc(N * sizeof(int));
//     int *line_temp = new int[N];
//     split_Convert(line_temp, line, "\t");
//     // for (int i = 0; i < N; i++)
//     // {
//     //     cout << i << "\t" << line_temp[i] << endl;
//     // }

//     int *cuda_line_Data;
//     cudaMallocManaged(&cuda_line_Data, N * sizeof(int));

//     int *differences, *cuda_Differences;
//     cudaMallocManaged(&cuda_Differences, N * sizeof(int));
//     differences = (int *)malloc(N * sizeof(int));

//     cudaMemcpy(cuda_line_Data, line_temp, (N * sizeof(int)), cudaMemcpyHostToDevice);

//     pairwise_Cuda<<<tot_Blocks, tot_ThreadsperBlock>>>(N, cuda_line_Data, cuda_Differences);
//     cudaDeviceSynchronize();

//     cudaMemcpy(differences, cuda_Differences, N * sizeof(int), cudaMemcpyDeviceToHost);

//     cudaFree(cuda_line_Data);
//     cudaFree(cuda_Differences);

//     for (int i = 0; i < N; i++)
//     {
//         pairwise_Differences = pairwise_Differences + differences[i];
//     }

//     // cout << "pairwise: " << pairwise_Differences << endl;

//     free(differences);
//     free(line_temp);

//     return pairwise_Differences;
// }

// void tajima::split_getPos(vector<string> &line_Data, string line, string delim)
// {
//     vector<string>().swap(line_Data);
//     char *convert;
//     string capture(line);
//     convert = &capture[0];
//     // cout<<convert;

//     char deliminator[delim.length() + 1];
//     strcpy(deliminator, delim.c_str());

//     char *split_data;
//     split_data = strtok(convert, deliminator);
//     int count = 0;

//     while (split_data != NULL)
//     {
//         // cout<<split_data<<endl;
//         string char2string;
//         char2string.append(split_data);
//         // cout << char2string << endl;
//         line_Data.push_back(char2string);
//         if (count == 7)
//         {
//             break;
//         }
//         split_data = strtok(NULL, deliminator);
//         count++;
//     }
// }

// long int tajima::combos_N(int count)
// {
//     long int combinations;

//     combinations = fact_half(count) / 2;

//     return combinations;
// }

// long int tajima::fact_half(int count)
// {
//     long int tot = 1;
//     for (int i = count; i > count - 2; i--)
//     {
//         // cout << tot;
//         tot = tot * i;
//     }
//     return tot;
// }

// vector<string> tajima::compound_interpolationSearch(vector<pair<string, string>> &folder_Index, int &start_Co, int &end_Co)
// {
//     vector<string> file_List;

//     vector<string> line_Data;
//     split(line_Data, folder_Index[0].first, "_");
//     int low_Value = stoi(line_Data[0]);
//     split(line_Data, folder_Index[folder_Index.size() - 1].first, "_");
//     int high_Value = stoi(line_Data[1]);
//     // cout << "first: " << low_Value << " last: " << high_Value << endl;

//     int start = 0;
//     int end = folder_Index.size() - 1;

//     while (start <= end && start_Co >= low_Value && start_Co <= high_Value)
//     {
//         vector<string> line_Data_get;

//         int pos = start + ((((double)(end - start) / (high_Value - low_Value)) * (start_Co - low_Value)));
//         // cout << pos << endl;

//         split(line_Data_get, folder_Index[pos].first, "_");
//         int low_Value_atpos = stoi(line_Data_get[0]);
//         int high_Value_atpos = stoi(line_Data_get[1]);

//         if ((start_Co >= low_Value_atpos && start_Co <= high_Value_atpos) || (start_Co <= low_Value_atpos && end_Co <= high_Value_atpos) || (start_Co <= low_Value_atpos && end_Co >= high_Value_atpos))
//         {
//             // cout << low_Value_atpos << "_" << high_Value_atpos << endl;
//             // backward_Search(pos, backward_get, folder_Index, start_Co, end_Co);

//             // thread backward(&tajima::backward_Search, this, pos, ref(backward_get), ref(folder_Index), ref(start_Co), ref(end_Co));
//             // thread forward(&tajima::forward_Search, this, pos, ref(forward_get), ref(folder_Index), ref(start_Co), ref(end_Co));

//             // backward.join();
//             // forward.join();

//             vector<int> backward_get;
//             vector<int> forward_get;

//             future<vector<int>> backward_thread = async(&tajima::backward_Search, this, pos, folder_Index, start_Co, end_Co);
//             future<vector<int>> forward_thread = async(&tajima::forward_Search, this, pos, folder_Index, start_Co, end_Co);

//             backward_get = backward_thread.get();
//             forward_get = forward_thread.get();

//             for (auto positions : backward_get)
//             {
//                 file_List.push_back(folder_Index[positions].second);
//             }

//             // cout << "caught :" << folder_Index[pos].second << endl;
//             file_List.push_back(folder_Index[pos].second);

//             for (auto positions : forward_get)
//             {
//                 file_List.push_back(folder_Index[positions].second);
//             }

//             break;
//         }
//         else if (start_Co > low_Value_atpos)
//         {
//             start = pos + 1;
//         }
//         else
//         {
//             end = pos - 1;
//         }

//         split(line_Data_get, folder_Index[start].first, "_");
//         low_Value = stoi(line_Data_get[0]);

//         split(line_Data_get, folder_Index[end].first, "_");
//         high_Value = stoi(line_Data_get[1]);
//     }

//     return file_List;
// }

// vector<int> tajima::backward_Search(int pos, vector<pair<string, string>> folder_Index, int start_Co, int end_Co)
// {
//     vector<int> backward_get;
//     pos = pos - 1;
//     vector<string> line_Data_get;
//     while (pos >= 0)
//     {
//         split(line_Data_get, folder_Index[pos].first, "_");
//         int low_Value_atpos = stoi(line_Data_get[0]);
//         int high_Value_atpos = stoi(line_Data_get[1]);

//         if (start_Co > high_Value_atpos)
//         {
//             break;
//         }

//         if ((start_Co >= low_Value_atpos && start_Co <= high_Value_atpos) || (start_Co <= low_Value_atpos && end_Co <= high_Value_atpos) || (start_Co <= low_Value_atpos && end_Co >= high_Value_atpos))
//         {
//             backward_get.push_back(pos);
//         }

//         pos = pos - 1;
//     }
//     return backward_get;
// }

// vector<int> tajima::forward_Search(int pos, vector<pair<string, string>> folder_Index, int start_Co, int end_Co)
// {
//     vector<int> forward_get;
//     pos = pos + 1;
//     vector<string> line_Data_get;
//     while (pos < folder_Index.size())
//     {
//         split(line_Data_get, folder_Index[pos].first, "_");
//         int low_Value_atpos = stoi(line_Data_get[0]);
//         int high_Value_atpos = stoi(line_Data_get[1]);

//         if (end_Co < low_Value_atpos)
//         {
//             break;
//         }

//         if ((start_Co >= low_Value_atpos && start_Co <= high_Value_atpos) || (start_Co <= low_Value_atpos && end_Co <= high_Value_atpos) || (start_Co <= low_Value_atpos && end_Co >= high_Value_atpos))
//         {
//             forward_get.push_back(pos);
//         }

//         pos = pos + 1;
//     }
//     return forward_get;
// }

__global__ void a_Calculation(int N, float *a1_CUDA, float *a2_CUDA)
{
    /**
     * Generate each value in the sequence for both a1 and a2
     **/
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < N)
    {
        a1_CUDA[tid] = (float)1 / (tid + 1);
        a2_CUDA[tid] = (float)1 / ((tid + 1) * (tid + 1));
        tid += blockDim.x * gridDim.x;
    }
}

__global__ void add_Cuda(const float *a, float *out, int arraySize)
{
    /**
     * CUDA based reduction add function for large arrays
     **/

    int idx = threadIdx.x;
    float sum = 0;
    for (int i = idx; i < arraySize; i += 1024)
        sum += a[i];
    __shared__ int r[1024];
    r[idx] = sum;
    __syncthreads();
    for (int size = 1024 / 2; size > 0; size /= 2)
    { // uniform
        if (idx < size)
            r[idx] += r[idx + size];
        __syncthreads();
    }
    if (idx == 0)
        *out = r[0];
}

void tajima::calc_Pre(int &N_tot, float &a1, float &e1, float &e2)
{
    /**
     * Calculates the prerequisite values required for Tajima's D
     **/

    int N = N_tot - 1;
    float *a1_CUDA, *a2_CUDA;
    float *a1_partial, *a2_partial, a2;

    a1_partial = (float *)malloc(N * sizeof(float));
    a2_partial = (float *)malloc(N * sizeof(float));

    cudaMallocManaged(&a1_CUDA, N * sizeof(int));
    cudaMallocManaged(&a2_CUDA, N * sizeof(int));

    a_Calculation<<<tot_Blocks, tot_ThreadsperBlock>>>(N, a1_CUDA, a2_CUDA);
    cudaDeviceSynchronize();

    cudaMemcpy(a1_partial, a1_CUDA, N * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(a2_partial, a2_CUDA, N * sizeof(float), cudaMemcpyDeviceToHost);

    cudaFree(a1_CUDA);
    cudaFree(a2_CUDA);

    /**
     * Summation of the harmonic means
     **/

    a1 = 0;
    a2 = 0;
    for (size_t i = 0; i < N; i++)
    {
        a1 += a1_partial[i];
        a2 += a2_partial[i];
    }

    cout << "a1\t: " << a1 << "\t"
         << "a2: " << a2 << endl;

    free(a1_partial);
    free(a2_partial);

    float b1 = (float)(N_tot + 1) / (3 * (N_tot - 1));
    float b2 = (float)(2 * (pow(N_tot, 2.0) + N_tot + 3)) / ((9 * N_tot) * (N_tot - 1));

    cout << "b1\t: " << b1 << "\t"
         << "b2: " << b2 << endl;

    float c1 = b1 - (1 / a1);
    float c2 = b2 - ((N_tot + 2) / (a1 * N_tot)) + (a2 / pow(a1, 2.0));

    cout << "c1\t: " << c1 << "\t"
         << "c2: " << c2 << endl;

    e1 = c1 / a1;
    e2 = c2 / (pow(a1, 2.0) + a2);

    cout << "e1\t: " << e1 << "\t"
         << "e2: " << e2 << endl;

    cout << endl;
}

// int tajima::getN_Split(string file)
// {
//     fstream file_nCount;
//     file_nCount.open(file);

//     string header;
//     getline(file_nCount, header);
//     file_nCount.close();

//     vector<string> header_Columns;
//     split(header_Columns, header, "\t");

//     int N = header_Columns.size() - 9;
//     // cout << N << endl;
//     return N;
// }

// vector<pair<string, string>> tajima::index_Folder(string &country)
// {
//     cout << "Initiating indexing folder\t: " << country << endl;
//     string country_Only = country.substr(country.find_last_of("/") + 1, country.length());

//     vector<string> index_pass_1;
//     vector<pair<string, string>> file_coordinate;

//     for (const auto &entry : filesystem::directory_iterator(country))
//     {
//         // cout << entry.path() << endl;
//         string coordinates = entry.path().string();
//         int trim_start = coordinates.find(country_Only + "_") + country_Only.length() + 1;
//         int trim_end = coordinates.find_last_of(".") - trim_start;
//         string trim = coordinates.substr(trim_start, trim_end);
//         index_pass_1.push_back(trim);
//         file_coordinate.push_back(make_pair(trim, coordinates));
//     }

//     vector<pair<int, int>> start_stop;
//     vector<int> starts;

//     for (string file : index_pass_1)
//     {
//         vector<string> file_start_end;
//         split(file_start_end, file, "_");
//         starts.push_back(stoi(file_start_end[0]));
//         start_stop.push_back(make_pair(stoi(file_start_end[0]), stoi(file_start_end[1])));
//     }

//     sort(starts.begin(), starts.end());

//     vector<pair<string, string>> sorted_Index;

//     for (int nums : starts)
//     {
//         // cout << nums << endl;
//         for (auto index_check : start_stop)
//         {
//             if (index_check.first == nums)
//             {
//                 string sort_Line = to_string(index_check.first) + "_" + to_string(index_check.second);
//                 for (auto coordinates : file_coordinate)
//                 {
//                     if (coordinates.first == sort_Line)
//                     {
//                         sorted_Index.push_back(make_pair(sort_Line, coordinates.second));
//                         break;
//                     }
//                 }
//                 break;
//             }
//         }
//     }

//     return sorted_Index;
// }

// vector<string> tajima::get_Countries()
// {
//     vector<string> folders;
//     for (auto &check : std::filesystem::recursive_directory_iterator(this->input_Folder))
//     {
//         if (check.is_directory())
//         {
//             folders.push_back(check.path().string());
//         }
//     }
//     return folders;
// }

// void tajima::split(vector<string> &line_Data, string line, string delim)
// {
//     vector<string>().swap(line_Data);
//     char *convert;
//     string capture(line);
//     convert = &capture[0];
//     // cout<<convert;

//     char deliminator[delim.length() + 1];
//     strcpy(deliminator, delim.c_str());

//     char *split_data;
//     split_data = strtok(convert, deliminator);

//     while (split_data != NULL)
//     {
//         // cout<<split_data<<endl;
//         string char2string;
//         char2string.append(split_data);
//         // cout << char2string << endl;
//         line_Data.push_back(char2string);
//         split_data = strtok(NULL, deliminator);
//     }

//     // delete convert;
//     // delete split_data;
// }

// void tajima::split_Convert(int *line_temp, string line, string delim)
// {
//     char *convert;
//     string capture(line);
//     convert = &capture[0];
//     // cout<<convert;

//     char deliminator[delim.length() + 1];
//     strcpy(deliminator, delim.c_str());

//     char *split_data;
//     split_data = strtok(convert, deliminator);
//     int count = 0;

//     while (split_data != NULL)
//     {
//         // cout<<split_data<<endl;
//         string char2string;
//         char2string.append(split_data);
//         // cout << char2string << "\t";
//         // cout << char2string << endl;
//         // line_Data.push_back(char2string);
//         if (count >= 9)
//         {
//             // cout << "count : " << char2string << endl;
//             if (char2string.substr(0, 3) == "0|0")
//             {
//                 // cout << "0|0" << count - 9 << endl;
//                 line_temp[count - 9] = 0;
//             }
//             else if (char2string.substr(0, 3) == "0|1")
//             {
//                 // cout << "0|1" << count - 9 << endl;
//                 line_temp[count - 9] = 1;
//             }
//             else if (char2string.substr(0, 3) == "1|0")
//             {
//                 line_temp[count - 9] = 2;
//             }
//             else if (char2string.substr(0, 3) == "1|1")
//             {
//                 line_temp[count - 9] = 3;
//             }
//             // cout << "count : " << count - 9 << "\t" << line_temp[count - 9] << endl;
//         }
//         count++;
//         split_data = strtok(NULL, deliminator);
//     }
//     // cout << endl;
//     // delete convert;
//     // delete split_data;
// }