#include "functions.cuh"
#include "vcf_splitter_2.cuh"

vcf_splitter_2::vcf_splitter_2(int cuda_ID, string input_vcf_Folder, string output_Folder, int cores, int SNPs_per_time_CPU, int SNPs_per_time_GPU, int allele_Count_REF, int allele_Count_ALT, int ploidy, int summary_Individuals)
{
    // ! THIS CONSTRUCTOR IS FOR THE BY CHROMOSOME SPLIT
    // trim to GT as well.
    // trim by REF and ALT allele count as well
    // split by chromosome

    cout << "Starting up VCF SPLITTER: Split by chromosome\n\nNOTE: Only the GT column will be extracted." << endl
         << endl;

    this->input_vcf_Folder = input_vcf_Folder;
    this->output_Folder = output_Folder;

    this->cores = cores;
    this->SNPs_per_time_CPU = SNPs_per_time_CPU;
    this->SNPs_per_time_GPU = SNPs_per_time_GPU;

    this->allele_Count_REF = allele_Count_REF;
    this->allele_Count_ALT = allele_Count_ALT;

    this->ploidy = ploidy;
    this->hap_Size = (2 * ploidy) - 1;

    this->summary_Individuals = summary_Individuals;

    cuda_Set_device(cuda_ID);
}

vcf_splitter_2::vcf_splitter_2(int cuda_ID, string input_vcf_Folder, string output_Folder, string population_File, int sampled_ID_col, int pop_ID_column, int cores, int SNPs_per_time_CPU, int SNPs_per_time_GPU, int ploidy, int max_SNPs_per_file, int logic_MAF, double MAF)
{
    cout << "Starting up VCF SPLITTER: Creating file hierarchy\n\nNOTE: VCF files should have only the GT columns and must be SNPs." << endl
         << endl;

    this->input_vcf_Folder = input_vcf_Folder;
    this->output_Folder = output_Folder;

    this->cores = cores;
    this->SNPs_per_time_CPU = SNPs_per_time_CPU;
    this->SNPs_per_time_GPU = SNPs_per_time_GPU;

    this->ploidy = ploidy;
    this->hap_Size = (2 * ploidy) - 1;

    this->population_File_path = population_File;
    this->column_Sample_ID = sampled_ID_col;
    this->column_Population_ID = pop_ID_column;
    this->SNP_count_per_File = max_SNPs_per_file;

    this->logic_MAF = logic_MAF;
    this->MAF = MAF;

    cout << "MAF: " << MAF << endl;
    cout << "Logical comparison: ";
    if (logic_MAF == 0)
    {
        cout << "equal\n"
             << endl;
    }
    else if (logic_MAF == 1)
    {
        cout << "greater than\n"
             << endl;
    }
    else if (logic_MAF == 2)
    {
        cout << "less than\n"
             << endl;
    }
    else if (logic_MAF == 10)
    {
        cout << "greater than or equal\n"
             << endl;
    }
    else if (logic_MAF == 20)
    {
        cout << "less than or equal\n"
             << endl;
    }

    cuda_Set_device(cuda_ID);
}

void vcf_splitter_2::cuda_Set_device(int cuda_ID)
{
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

void vcf_splitter_2::ingress_file_hierarchy()
{
    functions function = functions();

    fstream population_File;
    population_File.open(this->population_File_path, ios::in);

    vector<pair<string, string>> sample_population;
    set<string> population_Unique_IDs_SET;

    if (population_File.is_open())
    {
        cout << "Processing population file: " << this->population_File_path << endl;

        string line;
        getline(population_File, line);

        while (getline(population_File, line))
        {
            vector<string> data_Line;
            function.split(data_Line, line, '\t');

            sample_population.push_back(make_pair(data_Line[column_Sample_ID - 1], data_Line[column_Population_ID - 1]));
            population_Unique_IDs_SET.insert(data_Line[column_Population_ID - 1]);
        }

        population_File.close();
    }

    int num_Unique_populations = population_Unique_IDs_SET.size();
    cout << num_Unique_populations << " unique population(s) were found: ";

    for (string ID : population_Unique_IDs_SET)
    {
        vector<int> row;
        cout << ID;
        population_Unique_IDs.push_back(ID);
        population_sample_IDs.push_back(row);
        if (population_Unique_IDs.size() != num_Unique_populations)
        {
            cout << ", ";
        }
    }

    population_Unique_IDs_SET.clear();

    cout << endl
         << endl;

    cout << "Identifying VCF file(s): " << this->input_vcf_Folder << endl;

    vector<pair<string, string>> vcf_Files;

    for (const auto &entry : filesystem::directory_iterator(input_vcf_Folder))
    {
        string coordinates = entry.path().string();
        string extension = coordinates.substr(coordinates.find_last_of(".") + 1);
        transform(extension.begin(), extension.end(), extension.begin(), ::toupper);
        if (extension == "VCF")
        {
            vcf_Files.push_back(make_pair(coordinates, filesystem::path(coordinates).stem().string()));
        }
        // cout << coordinates << "\t"
        //      << extension << endl;
    }

    cout << vcf_Files.size() << " VCF file(s) have been found\n"
         << endl;

    for (int vcf_Index = 0; vcf_Index < vcf_Files.size(); vcf_Index++)
    {
        fstream file;
        file.open(vcf_Files[vcf_Index].first, ios::in);

        if (file.is_open())
        {
            cout << "Processing file: " << vcf_Files[vcf_Index].first << endl
                 << endl;
            string line;
            getline(file, line);
            while (line.substr(0, 2) == "##")
            {
                getline(file, line);
            }

            function.split(header_Data, line, '\t');

            for (int erase = 0; erase < 9; erase++)
            {
                header_Data.erase(header_Data.begin());
            }

            int N = header_Data.size();
            int augment = ((2 * ploidy) - 1) * N;

            cout << "Found " << N << " samples in VCF" << endl;
            cout << "Mapping samples to populations" << endl
                 << endl;

            sample_ID_population_ID = (int *)malloc(N * sizeof(int));

            vector<thread> threads_vec;

            int N_per_Thread = N / cores;
            int remainder = N % cores;

            for (int core_ID = 0; core_ID < cores; core_ID++)
            {
                int start_N = core_ID * N_per_Thread;
                int stop_N = start_N + N_per_Thread;

                threads_vec.push_back(thread{&vcf_splitter_2::individual_Map, this, start_N, stop_N, header_Data, sample_population, population_Unique_IDs});
            }

            if (remainder != 0)
            {
                int start_N = N - remainder;
                int stop_N = N;

                threads_vec.push_back(thread{&vcf_splitter_2::individual_Map, this, start_N, stop_N, header_Data, sample_population, population_Unique_IDs});
            }

            for (thread &t : threads_vec)
            {
                if (t.joinable())
                {
                    t.join();
                }
            }

            threads_vec.clear();
            sample_population.clear();

            // for (size_t i = 0; i < N; i++)
            // {
            //     cout << sample_ID_population_ID[i] << " ";
            // }
            // cout << endl;
            // for (size_t i = 0; i < N; i++)
            // {
            //     cout << sample_ID_population_ID[i] << ", ";
            // }

            cout << "Creating directories" << endl;
            string output_vcf_Folder = output_Folder + "/" + vcf_Files[vcf_Index].second;
            cout << "Creating root folder: " << output_vcf_Folder << endl;
            filesystem::create_directory(output_vcf_Folder);

            cout << "\nCreating population folders and calculating MAF: \n\n";

            cudaMallocManaged(&cuda_MAF_count_per_Population, population_Unique_IDs.size() * sizeof(int));
            MAF_count_per_Population = (int *)malloc(population_Unique_IDs.size() * sizeof(int));
            cudaMallocManaged(&cuda_sample_ID_population_ID, N * sizeof(int));

            for (int pop_ID = 0; pop_ID < population_Unique_IDs.size(); pop_ID++)
            {
                cout << "Population: " << population_Unique_IDs[pop_ID] << endl;
                string population_Folder = output_vcf_Folder + "/" + population_Unique_IDs[pop_ID];
                cout << "Folder created: " << population_Folder << endl;

                string header_String = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
                filesystem::create_directory(population_Folder);

                sort(population_sample_IDs[pop_ID].begin(), population_sample_IDs[pop_ID].end());

                for (int sample = 0; sample < population_sample_IDs[pop_ID].size(); sample++)
                {
                    header_String = header_String + header_Data[population_sample_IDs[pop_ID][sample]];

                    if (sample != population_sample_IDs[pop_ID].size() - 1)
                    {
                        header_String = header_String + "\t";
                    }
                }

                pop_Header.push_back(header_String);

                int MAF = population_sample_IDs[pop_ID].size() * this->ploidy * this->MAF;
                cout << "MAF count (" << population_sample_IDs[pop_ID].size()
                     << " x " << this->ploidy << " x " << this->MAF << "): "
                     << MAF << endl
                     << endl;
                MAF_count_per_Population[pop_ID] = MAF;
            }

            cudaMemcpy(cuda_MAF_count_per_Population, MAF_count_per_Population, population_Unique_IDs.size() * sizeof(int), cudaMemcpyHostToDevice);
            cudaMemcpy(cuda_sample_ID_population_ID, sample_ID_population_ID, N * sizeof(int), cudaMemcpyHostToDevice);

            // cudaError_t err = cudaGetLastError();
            // if (err != cudaSuccess)
            // {
            //     printf("CUDA Error 1: %s\n", cudaGetErrorString(err));

            //     // Possibly: exit(-1) if program cannot continue....
            // }
            // cudaDeviceSynchronize();

            while (getline(file, line))
            {
                all_Lines.push_back(line);

                if (all_Lines.size() == this->SNPs_per_time_CPU)
                {
                    process_SNPs_Hierarchy(N, population_Unique_IDs.size(), augment, output_vcf_Folder);
                }
            }

            file.close();

            // exit(3);

            if (all_Lines.size() != 0)
            {
                process_SNPs_Hierarchy(N, population_Unique_IDs.size(), augment, output_vcf_Folder);
            }

            for (int pop_ID = 0; pop_ID < population_Unique_IDs.size(); pop_ID++)
            {
                if (position_Complete_segs[pop_ID].size() > 0)
                {
                    string pop_Name = population_Unique_IDs[pop_ID];
                    string population_Folder = output_vcf_Folder + "/" + pop_Name;

                    sort(position_Complete_segs[pop_ID].begin(), position_Complete_segs[pop_ID].end());

                    vector<pair<int, string>> Segs_population = position_Complete_segs[pop_ID];

                    int pos_1 = Segs_population[0].first;
                    int pos_last = Segs_population[Segs_population.size() - 1].first;

                    string file_Name = population_Folder + "/" + file_CHR_value + "_" + pop_Name + "_" + to_string(pos_1) + "_" + to_string(pos_last) + ".vcf";
                    cout << "Writing file segment: " << file_Name << endl;
                    function.createFile(file_Name, pop_Header[pop_ID]);

                    fstream file_Segment;
                    file_Segment.open(file_Name, ios::app);

                    for (int line = 0; line < Segs_population.size(); line++)
                    {
                        file_Segment << Segs_population[line].second + "\n";
                    }

                    file_Segment.close();
                }
            }
        }

        free(sample_ID_population_ID);
        free(MAF_count_per_Population);

        cudaFree(cuda_sample_ID_population_ID);
        cudaFree(cuda_MAF_count_per_Population);

        header_Data.clear();
        pop_Header.clear();
        population_sample_IDs.clear();
        population_Unique_IDs.clear();
        position_Complete_segs.clear();
    }
    cout << endl;
}

void vcf_splitter_2::individual_Map(int start_N, int stop_N, vector<string> header_Data, vector<pair<string, string>> sample_population, vector<string> population_Unique_IDs)
{
    vector<pair<int, int>> sample_ID_population_ID_partial;

    for (int N_ID = start_N; N_ID < stop_N; N_ID++)
    {
        string sample_Name = header_Data[N_ID];
        for (int check_Sample = 0; check_Sample < sample_population.size(); check_Sample++)
        {
            if (sample_Name == sample_population[check_Sample].first)
            {
                string country = sample_population[check_Sample].second;
                for (int check_Country = 0; check_Country < population_Unique_IDs.size(); check_Country++)
                {
                    if (country == population_Unique_IDs[check_Country])
                    {
                        sample_ID_population_ID_partial.push_back(make_pair(N_ID, check_Country));
                        break;
                    }
                }
                sample_population.erase(sample_population.begin() + check_Sample);
                break;
            }
        }
    }

    unique_lock<shared_mutex> ul(g_mutex);
    for (size_t i = 0; i < sample_ID_population_ID_partial.size(); i++)
    {
        sample_ID_population_ID[sample_ID_population_ID_partial[i].first] = sample_ID_population_ID_partial[i].second;
        population_sample_IDs[sample_ID_population_ID_partial[i].second].push_back(sample_ID_population_ID_partial[i].first);
    }
}

__global__ void cuda_seg_Pop_process(char *sites, int *index, int num_Segregrating_sites, int ploidy, int tot_N_individuals, int *cuda_CHR_start_Index, int *cuda_CHR_end_Index, int *cuda_pos_start_Index, int *cuda_pos_end_Index, int *cuda_ID_start_Index, int *cuda_ID_end_Index, int *cuda_REF_start_Index, int *cuda_REF_end_Index, int *cuda_ALT_start_Index, int *cuda_ALT_end_Index, int *cuda_six_9_start_Index, int *cuda_six_9_end_Index, int *cuda_VALID_or_NOT, char *seg_Array, int num_pop, int **cuda_REF_populations, int **cuda_ALT_populations, int **cuda_VALID_or_NOT_populations, int *cuda_sample_ID_population_ID, int *cuda_MAF_count_per_Population, int logic_MAF)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < num_Segregrating_sites)
    {
        int column = 0;
        int site_Start = index[tid];
        int site_End = index[tid + 1];

        int i = site_Start;

        cuda_CHR_start_Index[tid] = site_Start;
        while (column < 1)
        {
            if (sites[i] == '\t')
            {
                column++;
            }
            i++;
        }
        cuda_CHR_end_Index[tid] = i - 1;

        // POS column
        cuda_pos_start_Index[tid] = i;
        while (column < 2)
        {
            if (sites[i] == '\t')
            {
                column++;
            }
            i++;
        }

        // will point to the tab but makes < easier later
        cuda_pos_end_Index[tid] = i - 1;

        cuda_ID_start_Index[tid] = i;
        while (column < 3)
        {
            if (sites[i] == '\t')
            {
                column++;
            }
            i++;
        }

        // will point to the tab but makes < easier later
        cuda_ID_end_Index[tid] = i - 1;

        cuda_REF_start_Index[tid] = i;
        while (column < 4)
        {
            if (sites[i] == '\t')
            {
                column++;
            }
            i++;
        }

        // will point to the tab but makes < easier later
        cuda_REF_end_Index[tid] = i - 1;

        int num_REF = cuda_REF_end_Index[tid] - cuda_REF_start_Index[tid];

        if (num_REF == 1)
        {
            cuda_ALT_start_Index[tid] = i;
            while (column < 5)
            {
                if (sites[i] == '\t')
                {
                    column++;
                }
                i++;
            }

            // will point to the tab but makes < easier later
            cuda_ALT_end_Index[tid] = i - 1;

            int num_ALT = cuda_ALT_end_Index[tid] - cuda_ALT_start_Index[tid];

            if (num_ALT == 1)
            {
                cuda_six_9_start_Index[tid] = i;
                while (column < 9)
                {
                    if (sites[i] == '\t')
                    {
                        column++;
                    }
                    i++;
                }

                // will point to the tab but makes < easier later
                cuda_six_9_end_Index[tid] = i - 1;

                // int all_pops = 0;
                cuda_VALID_or_NOT[tid] = 1;
                int sample_ID = 0;
                int pop_ID = cuda_sample_ID_population_ID[sample_ID];

                for (int pop = 0; pop < num_pop; pop++)
                {
                    cuda_REF_populations[tid][pop] = 0;
                    cuda_ALT_populations[tid][pop] = 0;
                    cuda_VALID_or_NOT_populations[tid][pop] = 1;
                }

                int seg_array_Start = (tid * ((2 * ploidy) - 1)) * tot_N_individuals;
                // // printf("%d\n", tid);

                while (i < site_End)
                {
                    if (sites[i] != '\t' && sites[i] != '\n')
                    {
                        seg_Array[seg_array_Start] = sites[i];
                        seg_array_Start++;
                        if (cuda_VALID_or_NOT_populations[tid][pop_ID] == 1)
                        {
                            if (sites[i] == '1')
                            {
                                cuda_ALT_populations[tid][pop_ID] = cuda_ALT_populations[tid][pop_ID] + 1;
                            }
                            else if (sites[i] == '0')
                            {
                                cuda_REF_populations[tid][pop_ID] = cuda_REF_populations[tid][pop_ID] + 1;
                            }
                            else if (sites[i] == '.')
                            {
                                cuda_VALID_or_NOT_populations[tid][pop_ID] = 0;
                            }
                        }
                    }
                    else
                    {
                        if (sites[i] == '\t')
                        {
                            sample_ID++;
                            pop_ID = cuda_sample_ID_population_ID[sample_ID];
                        }
                    }
                    i++;
                }

                int pop_Valid = 0;
                for (int pop = 0; pop < num_pop; pop++)
                {
                    if (cuda_VALID_or_NOT_populations[tid][pop] == 1)
                    {
                        int MA_Count;
                        if (cuda_REF_populations[tid][pop] > cuda_ALT_populations[tid][pop])
                        {
                            MA_Count = cuda_ALT_populations[tid][pop];
                        }
                        else
                        {
                            MA_Count = cuda_REF_populations[tid][pop];
                        }

                        // if (MA_Count != 0)
                        // {
                        if (logic_MAF == 0)
                        {
                            if (MA_Count == cuda_MAF_count_per_Population[pop])
                            {
                                // cuda_VALID_or_NOT_populations[tid][pop] = 1;
                            }
                            else
                            {
                                cuda_VALID_or_NOT_populations[tid][pop] = 0;
                                pop_Valid++;
                            }
                        }
                        else if (logic_MAF == 1)
                        {
                            if (MA_Count > cuda_MAF_count_per_Population[pop])
                            {
                                // cuda_VALID_or_NOT_populations[tid][pop] = 1;
                            }
                            else
                            {
                                cuda_VALID_or_NOT_populations[tid][pop] = 0;
                                pop_Valid++;
                            }
                        }
                        else if (logic_MAF == 2)
                        {
                            if (MA_Count < cuda_MAF_count_per_Population[pop])
                            {
                                //  cuda_VALID_or_NOT_populations[tid][pop] = 1;
                            }
                            else
                            {
                                cuda_VALID_or_NOT_populations[tid][pop] = 0;
                                pop_Valid++;
                            }
                        }
                        else if (logic_MAF == 10)
                        {
                            if (MA_Count >= cuda_MAF_count_per_Population[pop])
                            {
                                //  cuda_VALID_or_NOT_populations[tid][pop] = 1;
                            }
                            else
                            {
                                cuda_VALID_or_NOT_populations[tid][pop] = 0;
                                pop_Valid++;
                            }
                        }
                        else if (logic_MAF == 20)
                        {
                            if (MA_Count <= cuda_MAF_count_per_Population[pop])
                            {
                                // cuda_VALID_or_NOT_populations[tid][pop] = 1;
                            }
                            else
                            {
                                cuda_VALID_or_NOT_populations[tid][pop] = 0;
                                pop_Valid++;
                            }
                        }
                        // }
                        // else
                        // {
                        //     cuda_VALID_or_NOT_populations[tid][pop] = 0;
                        //     pop_Valid++;
                        // }
                    }
                    else
                    {
                        pop_Valid++;
                    }
                }
                if (pop_Valid == num_pop)
                {
                    cuda_VALID_or_NOT[tid] = 0;
                }
            }
            else
            {
                cuda_VALID_or_NOT[tid] = 0;
            }
        }
        else
        {
            cuda_VALID_or_NOT[tid] = 0;
        }

        if (cuda_VALID_or_NOT[tid] == 0)
        {
            int seg_array_Start = (tid * ((2 * ploidy) - 1)) * tot_N_individuals;
            int seg_array_Stop = seg_array_Start + (((2 * ploidy) - 1) * tot_N_individuals);

            for (size_t i = seg_array_Start; i < seg_array_Stop; i++)
            {
                seg_Array[i] = 'x';
            }
        }

        tid += blockDim.x * gridDim.x;
    }
}

void vcf_splitter_2::process_SNPs_Hierarchy(int N, int num_Populations, int augment, string output_vcf_Folder)
{
    functions function = functions();

    // cout << num_Populations << endl;

    /**
     * The GPU is permitted to handle only a certain max number of SNPs at a time.
     * Therefore the number of rounds of GPU processing and,
     * the range of SNPs that will be processed in each round will have to be determined.
     *
     * @param GPU_rounds_full rounds requiring the max set of SNPs to be processed.
     * @param GPU_rounds_partial rounds requiring the remaining set of SNPs to be processed.
     *
     * The start and stop range of each round is stored.
     **/

    int tot_Segs_Round = all_Lines.size();

    cout << "System is processing and filtering " << tot_Segs_Round << " segregating site(s)" << endl;

    int GPU_rounds_full = tot_Segs_Round / SNPs_per_time_GPU;
    int GPU_rounds_partial = tot_Segs_Round % SNPs_per_time_GPU;

    // vector<pair<int, int>> start_stop;
    for (int i = 0; i < GPU_rounds_full; i++)
    {
        int start = i * SNPs_per_time_GPU;
        int stop = start + SNPs_per_time_GPU;
        start_stop.push_back(make_pair(start, stop));
    }

    if (GPU_rounds_partial != 0)
    {
        int start = tot_Segs_Round - GPU_rounds_partial;
        int stop = tot_Segs_Round;
        start_stop.push_back(make_pair(start, stop));
    }

    vector<thread> threads_vec;

    /**
     * Concatenation of SNPs for GPU processing is also done in parallel,
     * Number of threads needed is based on the number of rounds needed.
     **/

    for (int rounds = 0; rounds < start_stop.size(); rounds++)
    {
        concat_Segs.push_back("");
        int *row;
        all_site_Index.push_back(row);
    }

    for (int rounds = 0; rounds < start_stop.size(); rounds++)
    {
        threads_vec.push_back(thread{&vcf_splitter_2::seg_Concat, this, rounds, start_stop[rounds].first, start_stop[rounds].second});
    }

    for (thread &t : threads_vec)
    {
        if (t.joinable())
        {
            t.join();
        }
    }

    threads_vec.clear();

    all_Lines.clear();

    for (int rounds = 0; rounds < start_stop.size(); rounds++)
    {
        string Seg_sites = concat_Segs[rounds];
        int *site_Index = all_site_Index[rounds];

        int start = start_stop[rounds].first;
        int stop = start_stop[rounds].second;

        int total_Segs = stop - start;

        char *full_Char;
        full_Char = (char *)malloc((Seg_sites.size() + 1) * sizeof(char));
        strcpy(full_Char, Seg_sites.c_str());

        char *cuda_full_Char;
        cudaMallocManaged(&cuda_full_Char, (Seg_sites.size() + 1) * sizeof(char));
        int *cuda_site_Index;
        cudaMallocManaged(&cuda_site_Index, (total_Segs + 1) * sizeof(int));
        cudaMemcpy(cuda_full_Char, full_Char, (Seg_sites.size() + 1) * sizeof(char), cudaMemcpyHostToDevice);
        cudaMemcpy(cuda_site_Index, site_Index, (total_Segs + 1) * sizeof(int), cudaMemcpyHostToDevice);

        int *chr_start_Index, *chr_end_Index, *cuda_chr_start_Index, *cuda_chr_end_Index;
        cudaMallocManaged(&cuda_chr_start_Index, total_Segs * sizeof(int));
        cudaMallocManaged(&cuda_chr_end_Index, total_Segs * sizeof(int));
        chr_start_Index = (int *)malloc(total_Segs * sizeof(int));
        chr_end_Index = (int *)malloc(total_Segs * sizeof(int));

        int *pos_start_Index, *pos_end_Index, *cuda_pos_start_Index, *cuda_pos_end_Index;
        cudaMallocManaged(&cuda_pos_start_Index, total_Segs * sizeof(int));
        cudaMallocManaged(&cuda_pos_end_Index, total_Segs * sizeof(int));
        pos_start_Index = (int *)malloc(total_Segs * sizeof(int));
        pos_end_Index = (int *)malloc(total_Segs * sizeof(int));

        int *ID_start_Index, *ID_end_Index, *cuda_ID_start_Index, *cuda_ID_end_Index;
        cudaMallocManaged(&cuda_ID_start_Index, total_Segs * sizeof(int));
        cudaMallocManaged(&cuda_ID_end_Index, total_Segs * sizeof(int));
        ID_start_Index = (int *)malloc(total_Segs * sizeof(int));
        ID_end_Index = (int *)malloc(total_Segs * sizeof(int));

        int *REF_start, *REF_stop, *ALT_start, *ALT_stop, *cuda_REF_start, *cuda_REF_stop, *cuda_ALT_start, *cuda_ALT_stop;

        cudaMallocManaged(&cuda_REF_start, total_Segs * sizeof(int));
        cudaMallocManaged(&cuda_REF_stop, total_Segs * sizeof(int));
        REF_start = (int *)malloc(total_Segs * sizeof(int));
        REF_stop = (int *)malloc(total_Segs * sizeof(int));

        cudaMallocManaged(&cuda_ALT_start, total_Segs * sizeof(int));
        cudaMallocManaged(&cuda_ALT_stop, total_Segs * sizeof(int));
        ALT_start = (int *)malloc(total_Segs * sizeof(int));
        ALT_stop = (int *)malloc(total_Segs * sizeof(int));

        int *six_9_start_Index, *six_9_stop_Index, *cuda_six_9_start_Index, *cuda_six_9_stop_Index;
        cudaMallocManaged(&cuda_six_9_start_Index, total_Segs * sizeof(int));
        cudaMallocManaged(&cuda_six_9_stop_Index, total_Segs * sizeof(int));
        six_9_start_Index = (int *)malloc(total_Segs * sizeof(int));
        six_9_stop_Index = (int *)malloc(total_Segs * sizeof(int));

        int *VALID_or_NOT, *cuda_VALID_or_NOT;
        VALID_or_NOT = (int *)malloc(total_Segs * sizeof(int));
        cudaMallocManaged(&cuda_VALID_or_NOT, total_Segs * sizeof(int));

        //  cudaError_t err4 = cudaGetLastError();
        // if (err4 != cudaSuccess)
        // {
        //     printf("CUDA Error 4: %s\n", cudaGetErrorString(err4));

        //     // Possibly: exit(-1) if program cannot continue....
        // }
        // cudaDeviceSynchronize();

        int **VALID_or_NOT_populations, **cuda_VALID_or_NOT_populations;
        int **cuda_REF_populations, **cuda_ALT_populations;

        VALID_or_NOT_populations = (int **)malloc(total_Segs * sizeof(int *));

        for (int i = 0; i < total_Segs; i++)
        {
            VALID_or_NOT_populations[i] = (int *)malloc((num_Populations + 1) * sizeof(int));
        }

        cudaMallocManaged(&cuda_VALID_or_NOT_populations, total_Segs * (num_Populations + 1) * sizeof(int));
        cudaMallocManaged(&cuda_REF_populations, total_Segs * (num_Populations + 1) * sizeof(int));
        cudaMallocManaged(&cuda_ALT_populations, total_Segs * (num_Populations + 1) * sizeof(int));

        // cudaError_t err5 = cudaGetLastError();
        // if (err5 != cudaSuccess)
        // {
        //     printf("CUDA Error 5: %s\n", cudaGetErrorString(err5));

        //     // Possibly: exit(-1) if program cannot continue....
        // }
        // cudaDeviceSynchronize();

        int **tmp = (int **)malloc(total_Segs * sizeof(tmp[0]));
        int **tmp_2 = (int **)malloc(total_Segs * sizeof(tmp_2[0]));
        int **tmp_3 = (int **)malloc(total_Segs * sizeof(tmp_3[0]));

        for (int i = 0; i < total_Segs; i++)
        {
            cudaMalloc((void **)&tmp[i], (num_Populations + 1) * sizeof(tmp[0][0]));
            cudaMalloc((void **)&tmp_2[i], (num_Populations + 1) * sizeof(tmp_2[0][0]));
            cudaMalloc((void **)&tmp_3[i], (num_Populations + 1) * sizeof(tmp_3[0][0]));
        }

        // cudaError_t err6 = cudaGetLastError();
        // if (err6 != cudaSuccess)
        // {
        //     printf("CUDA Error 6: %s\n", cudaGetErrorString(err6));

        //     // Possibly: exit(-1) if program cannot continue....
        // }
        // cudaDeviceSynchronize();

        cudaMemcpy(cuda_VALID_or_NOT_populations, tmp, total_Segs * sizeof(int *), cudaMemcpyHostToDevice);
        cudaMemcpy(cuda_REF_populations, tmp_2, total_Segs * sizeof(int *), cudaMemcpyHostToDevice);
        cudaMemcpy(cuda_ALT_populations, tmp_3, total_Segs * sizeof(int *), cudaMemcpyHostToDevice);

        // cudaError_t err2 = cudaGetLastError();
        // if (err2 != cudaSuccess)
        // {
        //     printf("CUDA Error 2: %s\n", cudaGetErrorString(err2));

        //     // Possibly: exit(-1) if program cannot continue....
        // }
        // cudaDeviceSynchronize();

        free(tmp);
        free(tmp_2);
        free(tmp_3);

        char *Hap_array, *cuda_Hap_array;
        cudaMallocManaged(&cuda_Hap_array, (((N * ((2 * ploidy) - 1)) * total_Segs) + 1) * sizeof(char));
        Hap_array = (char *)malloc((((N * ((2 * ploidy) - 1)) * total_Segs) + 1) * sizeof(char));

        // cout << N << endl;
        // cout << ploidy << endl;
        // cout << total_Segs << endl;

        cout << "\nRound " << rounds + 1 << ": System is extracting and filtering " << total_Segs << " segregating site(s)" << endl
             << endl;
        // MAF and no blank areas (samples with .) in seg for valid or not.
        // cuda_seg_Pop_process(char *sites, int *index, int num_Segregrating_sites, int ploidy, int tot_N_individuals, int *cuda_CHR_start_Index, int *cuda_CHR_end_Index, int *cuda_pos_start_Index, int *cuda_pos_end_Index, int *cuda_ID_start_Index, int *cuda_ID_end_Index, int *cuda_REF_start_Index, int *cuda_REF_end_Index, int *cuda_ALT_start_Index, int *cuda_ALT_end_Index, int *cuda_six_9_start_Index, int *cuda_six_9_end_Index, int *cuda_VALID_or_NOT, char *seg_Array, int num_pop, int **cuda_REF_populations, int **cuda_ALT_populations, int **cuda_VALID_or_NOT_populations, int *cuda_sample_ID_population_ID, int *cuda_MAF_count_per_Population, int logic_MAF)
        cuda_seg_Pop_process<<<tot_Blocks, tot_ThreadsperBlock>>>(cuda_full_Char, cuda_site_Index, total_Segs, ploidy, N, cuda_chr_start_Index, cuda_chr_end_Index, cuda_pos_start_Index, cuda_pos_end_Index, cuda_ID_start_Index, cuda_ID_end_Index, cuda_REF_start, cuda_REF_stop, cuda_ALT_start, cuda_ALT_stop, cuda_six_9_start_Index, cuda_six_9_stop_Index, cuda_VALID_or_NOT, cuda_Hap_array, num_Populations, cuda_REF_populations, cuda_ALT_populations, cuda_VALID_or_NOT_populations, cuda_sample_ID_population_ID, cuda_MAF_count_per_Population, logic_MAF);
        cudaError_t err3 = cudaGetLastError();
        if (err3 != cudaSuccess)
        {
            printf("CUDA Error: %s\n", cudaGetErrorString(err3));

            // Possibly: exit(-1) if program cannot continue....
        }
        cudaDeviceSynchronize();

        cout << "System is filtering valid segregating site(s)" << endl;

        cudaMemcpy(VALID_or_NOT, cuda_VALID_or_NOT, total_Segs * sizeof(int), cudaMemcpyDeviceToHost);

        cudaMemcpy(chr_start_Index, cuda_chr_start_Index, total_Segs * sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(chr_end_Index, cuda_chr_end_Index, total_Segs * sizeof(int), cudaMemcpyDeviceToHost);

        cudaMemcpy(pos_start_Index, cuda_pos_start_Index, total_Segs * sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(pos_end_Index, cuda_pos_end_Index, total_Segs * sizeof(int), cudaMemcpyDeviceToHost);

        cudaMemcpy(ID_start_Index, cuda_ID_start_Index, total_Segs * sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(ID_end_Index, cuda_ID_end_Index, total_Segs * sizeof(int), cudaMemcpyDeviceToHost);

        cudaMemcpy(REF_start, cuda_REF_start, total_Segs * sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(REF_stop, cuda_REF_stop, total_Segs * sizeof(int), cudaMemcpyDeviceToHost);

        cudaMemcpy(ALT_start, cuda_ALT_start, total_Segs * sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(ALT_stop, cuda_ALT_stop, total_Segs * sizeof(int), cudaMemcpyDeviceToHost);

        cudaMemcpy(six_9_start_Index, cuda_six_9_start_Index, total_Segs * sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(six_9_stop_Index, cuda_six_9_stop_Index, total_Segs * sizeof(int), cudaMemcpyDeviceToHost);

        cudaMemcpy(Hap_array, cuda_Hap_array, (((N * ((2 * ploidy) - 1)) * total_Segs) + 1) * sizeof(char), cudaMemcpyDeviceToHost);

        string seg_Full(Hap_array);

        for (size_t i = 0; i < total_Segs; i++)
        {
            seg_Collection.push_back("");
        }

        int segs_per_Thread = total_Segs / cores;
        int remainder = total_Segs % cores;

        // for (size_t i = 0; i < total_Segs; i++)
        // {
        //     cout << VALID_or_NOT[i] << " ";
        // }
        // cout << endl;

        // exit(0);

        for (int core_ID = 0; core_ID < cores; core_ID++)
        {
            int start_Seg = core_ID * segs_per_Thread;
            int stop_Seg = start_Seg + segs_per_Thread;

            threads_vec.push_back(thread{&vcf_splitter_2::seg_VALID_OR_NOT_list, this, start_Seg, stop_Seg, VALID_or_NOT, N, seg_Full, augment});
        }

        if (remainder != 0)
        {
            int start_Seg = total_Segs - remainder;
            int stop_Seg = total_Segs;

            threads_vec.push_back(thread{&vcf_splitter_2::seg_VALID_OR_NOT_list, this, start_Seg, stop_Seg, VALID_or_NOT, N, seg_Full, augment});
        }

        for (thread &t : threads_vec)
        {
            if (t.joinable())
            {
                t.join();
            }
        }

        threads_vec.clear();

        int valid_Count = VALID_List.size();

        cout << valid_Count << " valid segregating site(s) were found." << endl
             << endl;

        if (valid_Count > 0)
        {
            cout << "System is concatenating data from segregating site(s)" << endl;

            sort(VALID_List.begin(), VALID_List.end());

            for (size_t i = 0; i < valid_Count; i++)
            {
                write_CHR.push_back("");
                POS_valid.push_back(-1);
            }

            segs_per_Thread = valid_Count / cores;
            remainder = valid_Count % cores;

            for (int core_ID = 0; core_ID < cores; core_ID++)
            {
                int start_Seg = core_ID * segs_per_Thread;
                int stop_Seg = start_Seg + segs_per_Thread;
                // void concat_ALL_hierarchy(int start_Seg, int stop_Seg, vector<int> VALID_Segs, char *full_Char, int *chr_start_Index, int *chr_end_Index, int *pos_start_Index, int *pos_end_Index, int *ID_start_Index, int *ID_end_Index, int *REF_start, int *REF_stop, int *ALT_start, int *ALT_stop, int *six_9_start_Index, int *six_9_stop_Index);
                threads_vec.push_back(thread{&vcf_splitter_2::concat_ALL_hierarchy, this, start_Seg, stop_Seg, VALID_List, full_Char, chr_start_Index, chr_end_Index, pos_start_Index, pos_end_Index, ID_start_Index, ID_end_Index, REF_start, REF_stop, ALT_start, ALT_stop, six_9_start_Index, six_9_stop_Index});
            }

            if (remainder != 0)
            {
                int start_Seg = valid_Count - remainder;
                int stop_Seg = valid_Count;

                threads_vec.push_back(thread{&vcf_splitter_2::concat_ALL_hierarchy, this, start_Seg, stop_Seg, VALID_List, full_Char, chr_start_Index, chr_end_Index, pos_start_Index, pos_end_Index, ID_start_Index, ID_end_Index, REF_start, REF_stop, ALT_start, ALT_stop, six_9_start_Index, six_9_stop_Index});
            }

            for (thread &t : threads_vec)
            {
                if (t.joinable())
                {
                    t.join();
                }
            }

            threads_vec.clear();

            // for (string data : write_CHR)
            // {
            //     cout << data << endl;
            // }

            // write to population
            // copy population wise valid or NOT

            for (size_t i = 0; i < total_Segs; i++)
            {
                cudaMemcpy(VALID_or_NOT_populations[i], cuda_VALID_or_NOT_populations[i], (num_Populations + 1) * sizeof(cuda_VALID_or_NOT_populations[0][0]), cudaMemcpyDeviceToHost);
            }

            for (int pop_track = 0; pop_track < num_Populations; pop_track++)
            {
                vector<pair<int, string>> seg_Row;
                position_Complete_segs.push_back(seg_Row);
            }

            cout << "\nSystem is processing populations" << endl;

            for (int pop_track = 0; pop_track < num_Populations; pop_track++)
            {

                for (int core_ID = 0; core_ID < cores; core_ID++)
                {
                    int start_Seg = core_ID * segs_per_Thread;
                    int stop_Seg = start_Seg + segs_per_Thread;
                    // complete_Pop_segs(int pop_ID, int start_Seg, int stop_Seg, vector<int> VALID_Segs, int **VALID_or_NOT_populations, vector<int> sample_IDs_for_Population)
                    threads_vec.push_back(thread{&vcf_splitter_2::complete_Pop_segs, this, pop_track, start_Seg, stop_Seg, VALID_List, VALID_or_NOT_populations, population_sample_IDs[pop_track]});
                }

                if (remainder != 0)
                {
                    int start_Seg = valid_Count - remainder;
                    int stop_Seg = valid_Count;

                    threads_vec.push_back(thread{&vcf_splitter_2::complete_Pop_segs, this, pop_track, start_Seg, stop_Seg, VALID_List, VALID_or_NOT_populations, population_sample_IDs[pop_track]});
                }

                for (thread &t : threads_vec)
                {
                    if (t.joinable())
                    {
                        t.join();
                    }
                }

                threads_vec.clear();

                cout << "System has processed: " << population_Unique_IDs[pop_track] << "\t: " << position_Complete_segs[pop_track].size() << " valid sites." << endl;
            }

            write_CHR.clear();
            POS_valid.clear();
        }

        seg_Collection.clear();
        VALID_List.clear();

        free(site_Index);
        free(full_Char);

        cudaFree(cuda_full_Char);
        cudaFree(cuda_site_Index);

        //free(VALID_or_NOT_populations);

        for (int i = 0; i < total_Segs; i++)
        {
            free(VALID_or_NOT_populations[i]);
        }
        free(VALID_or_NOT_populations);
        
        cudaFree(cuda_VALID_or_NOT_populations);

        // for (int row = 0; row < total_Segs; row++)
        // {
        //     cudaFree(cuda_VALID_or_NOT_populations[row]);
        // }
        // // see
        // cudaFree(cuda_VALID_or_NOT_populations);

        free(chr_start_Index);
        free(chr_end_Index);
        cudaFree(cuda_chr_start_Index);
        cudaFree(cuda_chr_end_Index);

        free(pos_start_Index);
        free(pos_end_Index);
        cudaFree(cuda_pos_start_Index);
        cudaFree(cuda_pos_end_Index);

        free(ID_start_Index);
        free(ID_end_Index);
        cudaFree(cuda_ID_start_Index);
        cudaFree(cuda_ID_end_Index);

        free(REF_start);
        free(REF_stop);
        cudaFree(cuda_REF_start);
        cudaFree(cuda_REF_stop);

        free(ALT_start);
        free(ALT_stop);
        cudaFree(cuda_ALT_start);
        cudaFree(cuda_ALT_stop);

        free(six_9_start_Index);
        free(six_9_stop_Index);
        cudaFree(cuda_six_9_start_Index);
        cudaFree(cuda_six_9_stop_Index);

        // for (int row = 0; row < total_Segs; row++)
        // {
        //     cudaFree(cuda_REF_populations[row]);
        //     cudaFree(cuda_ALT_populations[row]);
        // }
        
        cudaFree(cuda_REF_populations);
        cudaFree(cuda_ALT_populations);

        free(Hap_array);
        cudaFree(cuda_Hap_array);

        // REMOVE after
        // break;
    }

    start_stop.clear();
    concat_Segs.clear();
    all_site_Index.clear();

    cout << endl;

    for (int pop_track = 0; pop_track < num_Populations; pop_track++)
    {
        if (position_Complete_segs[pop_track].size() >= this->SNP_count_per_File)
        {
            write_segs_to_Hierarchy(pop_track, output_vcf_Folder);
        }
    }

    cout << endl;
}

void vcf_splitter_2::write_segs_to_Hierarchy(int pop_ID, string output_vcf_Folder)
{
    functions function = functions();
    string pop_Name = population_Unique_IDs[pop_ID];
    string population_Folder = output_vcf_Folder + "/" + pop_Name;

    sort(position_Complete_segs[pop_ID].begin(), position_Complete_segs[pop_ID].end());

    vector<pair<int, string>> Segs_population = position_Complete_segs[pop_ID];

    int full_Rounds = Segs_population.size() / SNP_count_per_File;
    int partial_lines = Segs_population.size() % SNP_count_per_File;

    for (int round = 0; round < full_Rounds; round++)
    {
        int start_Value = round * SNP_count_per_File;
        int stop_Value = start_Value + SNP_count_per_File;

        int pos_1 = Segs_population[start_Value].first;
        int pos_last = Segs_population[stop_Value - 1].first;

        string file_Name = population_Folder + "/" + file_CHR_value + "_" + pop_Name + "_" + to_string(pos_1) + "_" + to_string(pos_last) + ".vcf";
        cout << "Writing file segment: " << file_Name << endl;
        function.createFile(file_Name, pop_Header[pop_ID]);

        fstream file_Segment;
        file_Segment.open(file_Name, ios::app);

        for (int line = start_Value; line < stop_Value; line++)
        {
            file_Segment << Segs_population[line].second + "\n";
        }

        file_Segment.close();
    }

    vector<pair<int, string>> remaining;

    if (partial_lines != 0)
    {
        for (size_t i = (Segs_population.size() - partial_lines); i < Segs_population.size(); i++)
        {
            remaining.push_back(Segs_population[i]);
        }
    }

    position_Complete_segs[pop_ID] = remaining;
}

void vcf_splitter_2::complete_Pop_segs(int pop_ID, int start_Seg, int stop_Seg, vector<int> VALID_Segs, int **VALID_or_NOT_populations, vector<int> sample_IDs_for_Population)
{
    vector<pair<int, string>> seg_Complete;
    int augment_per = (2 * ploidy) - 1;
    for (int seg = start_Seg; seg < stop_Seg; seg++)
    {
        int seg_Index = VALID_Segs[seg];
        if (VALID_or_NOT_populations[seg_Index][pop_ID] == 1)
        {
            string line = write_CHR[seg];
            string hap_Line = seg_Collection[seg_Index];
            for (int pop = 0; pop < sample_IDs_for_Population.size(); pop++)
            {
                int pop_Index = sample_IDs_for_Population[pop];
                int start_Pop = pop_Index * augment_per;
                for (int inc = start_Pop; inc < (start_Pop + augment_per); inc++)
                {
                    line = line + hap_Line.at(inc);
                }
                if (pop != sample_IDs_for_Population.size() - 1)
                {
                    line = line + "\t";
                }
            }
            seg_Complete.push_back(make_pair(POS_valid[seg], line));
        }
    }

    unique_lock<shared_mutex> ul(g_mutex);
    for (int i = 0; i < seg_Complete.size(); i++)
    {
        position_Complete_segs[pop_ID].push_back(seg_Complete[i]);
    }
}

void vcf_splitter_2::concat_ALL_hierarchy(int start_Seg, int stop_Seg, vector<int> VALID_Segs, char *full_Char, int *chr_start_Index, int *chr_end_Index, int *pos_start_Index, int *pos_end_Index, int *ID_start_Index, int *ID_end_Index, int *REF_start, int *REF_stop, int *ALT_start, int *ALT_stop, int *six_9_start_Index, int *six_9_stop_Index)
{
    vector<pair<int, string>> write_Line;
    string chr_value;
    vector<int> positions;
    for (int seg = start_Seg; seg < stop_Seg; seg++)
    {
        int index = VALID_Segs[seg];
        string chr_string = "";

        for (int i = chr_start_Index[index]; i < chr_end_Index[index]; i++)
        {
            chr_string = chr_string + full_Char[i];
        }

        string pos_string = "";
        for (int i = pos_start_Index[index]; i < pos_end_Index[index]; i++)
        {
            pos_string = pos_string + full_Char[i];
        }

        positions.push_back(stoi(pos_string));

        string ID_string = "";
        for (int i = ID_start_Index[index]; i < ID_end_Index[index]; i++)
        {
            ID_string = ID_string + full_Char[i];
        }

        string REF_string = "";
        for (int i = REF_start[index]; i < REF_stop[index]; i++)
        {
            REF_string = REF_string + full_Char[i];
        }

        string ALT_string = "";
        for (int i = ALT_start[index]; i < ALT_stop[index]; i++)
        {
            ALT_string = ALT_string + full_Char[i];
        }

        string six_9_string = "";
        for (int i = six_9_start_Index[index]; i < six_9_stop_Index[index]; i++)
        {
            six_9_string = six_9_string + full_Char[i];
        }

        if (seg == 0)
        {
            chr_value = chr_string;
        }

        write_Line.push_back(make_pair(seg, chr_string + "\t" + pos_string + "\t" + ID_string + "\t" + REF_string + "\t" + ALT_string + "\t" + six_9_string + "\t"));
    }

    unique_lock<shared_mutex> ul(g_mutex);
    for (size_t i = 0; i < write_Line.size(); i++)
    {
        write_CHR[write_Line[i].first] = write_Line[i].second;
        POS_valid[write_Line[i].first] = positions[i];

        if (write_Line[i].first == 0)
        {
            file_CHR_value = chr_value;
        }
    }
}

void vcf_splitter_2::ingress_chr_Split()
{
    // chr can be strings
    // provide a summary file

    functions function = functions();

    cout << "Identifying VCF file(s): " << this->input_vcf_Folder << endl;

    vector<pair<string, string>> vcf_Files;

    for (const auto &entry : filesystem::directory_iterator(input_vcf_Folder))
    {
        string coordinates = entry.path().string();
        string extension = coordinates.substr(coordinates.find_last_of(".") + 1);
        transform(extension.begin(), extension.end(), extension.begin(), ::toupper);
        if (extension == "VCF")
        {
            vcf_Files.push_back(make_pair(coordinates, filesystem::path(coordinates).stem().string()));
        }
        // cout << coordinates << "\t"
        //      << extension << endl;
    }

    cout << vcf_Files.size() << " VCF file(s) have been found\n"
         << endl;

    for (int vcf_Index = 0; vcf_Index < vcf_Files.size(); vcf_Index++)
    {
        // cout << vcf_Files[vcf_Index].first << "\t" << vcf_Files[vcf_Index].second << endl;

        fstream file;
        file.open(vcf_Files[vcf_Index].first, ios::in);

        if (file.is_open())
        {
            cout << "Processing file: " << vcf_Files[vcf_Index].first << endl
                 << endl;
            string line;
            getline(file, line);
            while (line.substr(0, 2) == "##")
            {
                getline(file, line);
            }

            // cout << line << endl;
            // process_Line for pools to get summary

            string head_Line = line;
            vector<string> header_Data;
            // vector<string> header_0_9_Columns;
            function.split(header_Data, head_Line, '\t');

            for (int erase = 0; erase < 9; erase++)
            {
                // header_0_9_Columns.push_back(header_Data[0]);
                header_Data.erase(header_Data.begin());
            }

            int N = header_Data.size();
            cout << "Found " << N << " sample(s) in VCF" << endl;

            int augment = ((2 * ploidy) - 1) * N;

            // for (string header : header_0_9_Columns)
            // {
            //     cout << header << " ";
            // }

            // create folder
            string output_vcf_Folder = output_Folder + "/" + vcf_Files[vcf_Index].second;
            cout << "Creating folder: " << output_vcf_Folder << endl;
            filesystem::create_directory(output_vcf_Folder);

            // vector<string> collect_Segregrating_sites;

            while (getline(file, line))
            {
                all_Lines.push_back(line);

                if (all_Lines.size() == this->SNPs_per_time_CPU)
                {
                    // process the sites

                    process_SNPs_CHR(N, augment, output_vcf_Folder, head_Line);

                    // break;

                    // clear everything before re run of loop here or in the process function
                    // collect_Segregrating_sites.clear();
                }
            }

            file.close();

            if (all_Lines.size() != 0)
            {
                process_SNPs_CHR(N, augment, output_vcf_Folder, head_Line);
            }

            if (summary_Individuals != 0)
            {
                string summary_Individuals_File = output_vcf_Folder + "/" + vcf_Files[vcf_Index].second + ".summary";
                cout << "\nWriting individuals summary file: " << summary_Individuals_File << endl
                     << endl;

                string header_Summary = "Individual_ID\t";

                for (int CHR_ID = 0; CHR_ID < CHR_individuals.size(); CHR_ID++)
                {
                    header_Summary = header_Summary + CHR_individuals[CHR_ID];
                    if (CHR_ID != (CHR_individuals.size() - 1))
                    {
                        header_Summary = header_Summary + "\t";
                    }
                }

                function.createFile(summary_Individuals_File, header_Summary);

                fstream summary_Out;
                summary_Out.open(summary_Individuals_File, ios::app);

                for (int individuals = 0; individuals < N; individuals++)
                {
                    summary_Out << header_Data[individuals] << "\t";
                    for (int chr = 0; chr < CHR_individuals.size(); chr++)
                    {
                        double ratio = (double)CHR_individuals_COUNT[chr][individuals] / (double)(CHR_full_Segs[chr] * this->ploidy);
                        summary_Out << to_string(ratio) << "\t";
                    }
                    summary_Out << "\n";
                }

                CHR_individuals.clear();
                CHR_individuals_COUNT.clear();
                CHR_full_Segs.clear();

                summary_Out.close();
            }

            // process if collect_Segregrating_sites has any remaining
        }
        // REMOVE AFTER TESTING
        // break;

        // unique_CHRs.clear();
    }
}

__global__ void cuda_seg_Info_extract(char *sites, int *index, int num_Segregrating_sites, int ploidy, int N_individuals, int *VALID_or_NOT, int REF_count, int ALT_count, int *cuda_CHR_start_Index, int *cuda_CHR_end_Index, int *cuda_pos_start_Index, int *cuda_pos_end_Index, int *cuda_ID_start_Index, int *cuda_ID_end_Index, int *cuda_REF_start_Index, int *cuda_REF_end_Index, int *cuda_ALT_start_Index, int *cuda_ALT_end_Index, int *cuda_six_8_start_Index, int *cuda_six_8_end_Index, int **sample_sequence_Tracker, char *seg_Array)
{

    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < num_Segregrating_sites)
    {
        int column = 0;
        int site_Start = index[tid];
        int site_End = index[tid + 1];

        int i = site_Start;

        cuda_CHR_start_Index[tid] = site_Start;
        while (column < 1)
        {
            if (sites[i] == '\t')
            {
                column++;
            }
            i++;
        }
        cuda_CHR_end_Index[tid] = i - 1;

        // POS column
        cuda_pos_start_Index[tid] = i;
        while (column < 2)
        {
            if (sites[i] == '\t')
            {
                column++;
            }
            i++;
        }

        // will point to the tab but makes < easier later
        cuda_pos_end_Index[tid] = i - 1;

        cuda_ID_start_Index[tid] = i;
        while (column < 3)
        {
            if (sites[i] == '\t')
            {
                column++;
            }
            i++;
        }

        // will point to the tab but makes < easier later
        cuda_ID_end_Index[tid] = i - 1;

        cuda_REF_start_Index[tid] = i;
        while (column < 4)
        {
            if (sites[i] == '\t')
            {
                column++;
            }
            i++;
        }

        // will point to the tab but makes < easier later
        cuda_REF_end_Index[tid] = i - 1;

        int num_REF = cuda_REF_end_Index[tid] - cuda_REF_start_Index[tid];

        if (REF_count == num_REF)
        {
            cuda_ALT_start_Index[tid] = i;
            while (column < 5)
            {
                if (sites[i] == '\t')
                {
                    column++;
                }
                i++;
            }

            // will point to the tab but makes < easier later
            cuda_ALT_end_Index[tid] = i - 1;

            int num_ALT = cuda_ALT_end_Index[tid] - cuda_ALT_start_Index[tid];

            if (ALT_count == num_ALT)
            {
                cuda_six_8_start_Index[tid] = i;

                while (column < 8)
                {
                    if (sites[i] == '\t')
                    {
                        column++;
                    }
                    i++;
                }

                cuda_six_8_end_Index[tid] = i - 1;

                int GT_found_column = -1;
                int GT_column = 0;

                while (column < 9)
                {
                    // GT pos
                    if (GT_found_column == -1)
                    {
                        if ((sites[i] == 'G' && sites[i + 1] == 'T') || (sites[i] == 'g' && sites[i + 1] == 't'))
                        {
                            GT_found_column = GT_column;
                        }

                        if (sites[i] == ':')
                        {
                            GT_column++;
                        }
                    }

                    if (sites[i] == '\t')
                    {
                        column++;
                    }
                    i++;
                }

                if (GT_found_column != 1)
                {
                    VALID_or_NOT[tid] = 1;

                    int within_column_Sample = 0;
                    int Sample_no = 0;
                    int sample_Sequences = 0;

                    int seg_array_Start = (tid * ((2 * ploidy) - 1)) * N_individuals;

                    while (i < site_End)
                    {
                        // GET GT data only and individuals data

                        if (sites[i] == ':')
                        {
                            within_column_Sample++;
                        }
                        else if (within_column_Sample == GT_found_column)
                        {
                            if (sites[i] != '\t')
                            {
                                seg_Array[seg_array_Start] = sites[i];
                                if (sites[i] != '|' && sites[i] != '/' && sites[i] != '.')
                                {
                                    sample_Sequences++;
                                }
                                seg_array_Start++;
                            }
                        }
                        if (sites[i] == '\t' || (i == (site_End - 1)))
                        {
                            sample_sequence_Tracker[tid][Sample_no] = sample_Sequences;

                            Sample_no++;

                            sample_Sequences = 0;
                            within_column_Sample = 0;
                        }
                        i++;
                    }
                }
                else
                {
                    VALID_or_NOT[tid] = 0;
                }
            }
            else
            {
                VALID_or_NOT[tid] = 0;
            }
        }
        else
        {
            VALID_or_NOT[tid] = 0;
        }

        if (VALID_or_NOT[tid] == 0)
        {
            int seg_array_Start = (tid * ((2 * ploidy) - 1)) * N_individuals;
            int seg_array_Stop = seg_array_Start + (((2 * ploidy) - 1) * N_individuals);

            for (size_t i = seg_array_Start; i < seg_array_Stop; i++)
            {
                seg_Array[i] = 'x';
            }
        }

        tid += blockDim.x * gridDim.x;
    }
}

void vcf_splitter_2::process_SNPs_CHR(int N, int augment, string &output_vcf_Folder, string &head_Line)
{
    functions function = functions();
    // rounds process

    /**
     * The GPU is permitted to handle only a certain max number of SNPs at a time.
     * Therefore the number of rounds of GPU processing and,
     * the range of SNPs that will be processed in each round will have to be determined.
     *
     * @param GPU_rounds_full rounds requiring the max set of SNPs to be processed.
     * @param GPU_rounds_partial rounds requiring the remaining set of SNPs to be processed.
     *
     * The start and stop range of each round is stored.
     **/

    int tot_Segs_Round = all_Lines.size();

    cout << "\nSystem is processing and filtering " << tot_Segs_Round << " segregating site(s)" << endl;

    int GPU_rounds_full = tot_Segs_Round / SNPs_per_time_GPU;
    int GPU_rounds_partial = tot_Segs_Round % SNPs_per_time_GPU;

    // vector<pair<int, int>> start_stop;
    for (int i = 0; i < GPU_rounds_full; i++)
    {
        int start = i * SNPs_per_time_GPU;
        int stop = start + SNPs_per_time_GPU;
        start_stop.push_back(make_pair(start, stop));
    }

    if (GPU_rounds_partial != 0)
    {
        int start = tot_Segs_Round - GPU_rounds_partial;
        int stop = tot_Segs_Round;
        start_stop.push_back(make_pair(start, stop));
    }

    vector<thread> threads_vec;

    /**
     * Concatenation of SNPs for GPU processing is also done in parallel,
     * Number of threads needed is based on the number of rounds needed.
     **/

    for (int rounds = 0; rounds < start_stop.size(); rounds++)
    {
        concat_Segs.push_back("");
        int *row;
        all_site_Index.push_back(row);
    }

    for (int rounds = 0; rounds < start_stop.size(); rounds++)
    {
        threads_vec.push_back(thread{&vcf_splitter_2::seg_Concat, this, rounds, start_stop[rounds].first, start_stop[rounds].second});
    }

    for (thread &t : threads_vec)
    {
        if (t.joinable())
        {
            t.join();
        }
    }

    threads_vec.clear();

    all_Lines.clear();

    for (int rounds = 0; rounds < start_stop.size(); rounds++)
    {
        string Seg_sites = concat_Segs[rounds];
        int *site_Index = all_site_Index[rounds];

        int start = start_stop[rounds].first;
        int stop = start_stop[rounds].second;

        int total_Segs = stop - start;

        char *full_Char;
        full_Char = (char *)malloc((Seg_sites.size() + 1) * sizeof(char));
        strcpy(full_Char, Seg_sites.c_str());

        char *cuda_full_Char;
        cudaMallocManaged(&cuda_full_Char, (Seg_sites.size() + 1) * sizeof(char));
        int *cuda_site_Index;
        cudaMallocManaged(&cuda_site_Index, (total_Segs + 1) * sizeof(int));
        cudaMemcpy(cuda_full_Char, full_Char, (Seg_sites.size() + 1) * sizeof(char), cudaMemcpyHostToDevice);
        cudaMemcpy(cuda_site_Index, site_Index, (total_Segs + 1) * sizeof(int), cudaMemcpyHostToDevice);

        int *cuda_VALID_or_NOT, *VALID_or_NOT;
        cudaMallocManaged(&cuda_VALID_or_NOT, total_Segs * sizeof(int));
        VALID_or_NOT = (int *)malloc(total_Segs * sizeof(int));

        int *chr_start_Index, *chr_end_Index, *cuda_chr_start_Index, *cuda_chr_end_Index;
        cudaMallocManaged(&cuda_chr_start_Index, total_Segs * sizeof(int));
        cudaMallocManaged(&cuda_chr_end_Index, total_Segs * sizeof(int));
        chr_start_Index = (int *)malloc(total_Segs * sizeof(int));
        chr_end_Index = (int *)malloc(total_Segs * sizeof(int));

        int *pos_start_Index, *pos_end_Index, *cuda_pos_start_Index, *cuda_pos_end_Index;
        cudaMallocManaged(&cuda_pos_start_Index, total_Segs * sizeof(int));
        cudaMallocManaged(&cuda_pos_end_Index, total_Segs * sizeof(int));
        pos_start_Index = (int *)malloc(total_Segs * sizeof(int));
        pos_end_Index = (int *)malloc(total_Segs * sizeof(int));

        int *ID_start_Index, *ID_end_Index, *cuda_ID_start_Index, *cuda_ID_end_Index;
        cudaMallocManaged(&cuda_ID_start_Index, total_Segs * sizeof(int));
        cudaMallocManaged(&cuda_ID_end_Index, total_Segs * sizeof(int));
        ID_start_Index = (int *)malloc(total_Segs * sizeof(int));
        ID_end_Index = (int *)malloc(total_Segs * sizeof(int));

        int *REF_start, *REF_stop, *ALT_start, *ALT_stop, *cuda_REF_start, *cuda_REF_stop, *cuda_ALT_start, *cuda_ALT_stop;

        cudaMallocManaged(&cuda_REF_start, total_Segs * sizeof(int));
        cudaMallocManaged(&cuda_REF_stop, total_Segs * sizeof(int));
        REF_start = (int *)malloc(total_Segs * sizeof(int));
        REF_stop = (int *)malloc(total_Segs * sizeof(int));

        cudaMallocManaged(&cuda_ALT_start, total_Segs * sizeof(int));
        cudaMallocManaged(&cuda_ALT_stop, total_Segs * sizeof(int));
        ALT_start = (int *)malloc(total_Segs * sizeof(int));
        ALT_stop = (int *)malloc(total_Segs * sizeof(int));

        int *six_8_start_Index, *six_8_stop_Index, *cuda_six_8_start_Index, *cuda_six_8_stop_Index;
        cudaMallocManaged(&cuda_six_8_start_Index, total_Segs * sizeof(int));
        cudaMallocManaged(&cuda_six_8_stop_Index, total_Segs * sizeof(int));
        six_8_start_Index = (int *)malloc(total_Segs * sizeof(int));
        six_8_stop_Index = (int *)malloc(total_Segs * sizeof(int));

        char *Hap_array, *cuda_Hap_array;
        cudaMallocManaged(&cuda_Hap_array, (((N * ((2 * ploidy) - 1)) * total_Segs) + 1) * sizeof(char));
        Hap_array = (char *)malloc((((N * ((2 * ploidy) - 1)) * total_Segs) + 1) * sizeof(char));

        int **cuda_sample_sequence_Tracker;
        cudaMallocManaged(&cuda_sample_sequence_Tracker, (N + 1) * total_Segs * sizeof(int));
        int **tmp = (int **)malloc(total_Segs * sizeof(tmp[0]));
        for (int i = 0; i < total_Segs; i++)
        {
            cudaMalloc((void **)&tmp[i], (N + 1) * sizeof(tmp[0][0]));
        }
        cudaMemcpy(cuda_sample_sequence_Tracker, tmp, total_Segs * sizeof(int *), cudaMemcpyHostToDevice);
        free(tmp);

        cout << "\nRound " << rounds + 1 << ": System is extracting and filtering " << total_Segs << " segregating site(s)" << endl
             << endl;

        // cuda_seg_Info_extract(char *sites, int *index, int num_Segregrating_sites, int ploidy, int N_individuals, int *VALID_or_NOT, int REF_count, int ALT_count, int *cuda_CHR_start_Index, int *cuda_CHR_end_Index, int *cuda_pos_start_Index, int *cuda_pos_end_Index, int *cuda_ID_start_Index, int *cuda_ID_end_Index, int *cuda_REF_start_Index, int *cuda_REF_end_Index, int *cuda_ALT_start_Index, int *cuda_ALT_end_Index, int *cuda_six_8_start_Index, int *cuda_six_8_end_Index, int **sample_sequence_Tracker, char *seg_Array)
        cuda_seg_Info_extract<<<tot_Blocks, tot_ThreadsperBlock>>>(cuda_full_Char, cuda_site_Index, total_Segs, ploidy, N, cuda_VALID_or_NOT, allele_Count_REF, allele_Count_ALT, cuda_chr_start_Index, cuda_chr_end_Index, cuda_pos_start_Index, cuda_pos_end_Index, cuda_ID_start_Index, cuda_ID_end_Index, cuda_REF_start, cuda_REF_stop, cuda_ALT_start, cuda_ALT_stop, cuda_six_8_start_Index, cuda_six_8_stop_Index, cuda_sample_sequence_Tracker, cuda_Hap_array);

        cudaError_t err = cudaGetLastError();

        if (err != cudaSuccess)
        {
            printf("CUDA Error: %s\n", cudaGetErrorString(err));

            // Possibly: exit(-1) if program cannot continue....
        }
        cudaDeviceSynchronize();

        cout << "System is filtering valid segregating site(s)" << endl;

        cudaMemcpy(VALID_or_NOT, cuda_VALID_or_NOT, total_Segs * sizeof(int), cudaMemcpyDeviceToHost);

        cudaMemcpy(chr_start_Index, cuda_chr_start_Index, total_Segs * sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(chr_end_Index, cuda_chr_end_Index, total_Segs * sizeof(int), cudaMemcpyDeviceToHost);

        cudaMemcpy(pos_start_Index, cuda_pos_start_Index, total_Segs * sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(pos_end_Index, cuda_pos_end_Index, total_Segs * sizeof(int), cudaMemcpyDeviceToHost);

        cudaMemcpy(ID_start_Index, cuda_ID_start_Index, total_Segs * sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(ID_end_Index, cuda_ID_end_Index, total_Segs * sizeof(int), cudaMemcpyDeviceToHost);

        cudaMemcpy(REF_start, cuda_REF_start, total_Segs * sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(REF_stop, cuda_REF_stop, total_Segs * sizeof(int), cudaMemcpyDeviceToHost);

        cudaMemcpy(ALT_start, cuda_ALT_start, total_Segs * sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(ALT_stop, cuda_ALT_stop, total_Segs * sizeof(int), cudaMemcpyDeviceToHost);

        cudaMemcpy(six_8_start_Index, cuda_six_8_start_Index, total_Segs * sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(six_8_stop_Index, cuda_six_8_stop_Index, total_Segs * sizeof(int), cudaMemcpyDeviceToHost);

        cudaMemcpy(Hap_array, cuda_Hap_array, (((N * ((2 * ploidy) - 1)) * total_Segs) + 1) * sizeof(char), cudaMemcpyDeviceToHost);

        int **sample_sequence_Tracker;
        if (this->summary_Individuals != 0)
        {
            sample_sequence_Tracker = (int **)malloc(total_Segs * sizeof(int *));

            for (size_t i = 0; i < total_Segs; i++)
            {
                sample_sequence_Tracker[i] = (int *)malloc((N + 1) * sizeof(int));
            }

            for (size_t i = 0; i < total_Segs; i++)
            {
                // cout << i << endl;
                cudaMemcpy(sample_sequence_Tracker[i], cuda_sample_sequence_Tracker[i], (N + 1) * sizeof(cuda_sample_sequence_Tracker[0][0]), cudaMemcpyDeviceToHost);
            }

            // for (size_t r = 0; r < total_Segs; r++)
            // {
            //     for (size_t c = 0; c < N; c++)
            //     {
            //         cout << sample_sequence_Tracker[r][c] << "\t";
            //     }
            //     break;
            // }
        }

        // exit(1);

        string seg_Full(Hap_array);
        // this->trim_Seg_full = seg_Full;

        for (size_t i = 0; i < total_Segs; i++)
        {
            seg_Collection.push_back("");
        }

        int segs_per_Thread = total_Segs / cores;
        int remainder = total_Segs % cores;

        for (int core_ID = 0; core_ID < cores; core_ID++)
        {
            int start_Seg = core_ID * segs_per_Thread;
            int stop_Seg = start_Seg + segs_per_Thread;

            threads_vec.push_back(thread{&vcf_splitter_2::seg_VALID_OR_NOT_list, this, start_Seg, stop_Seg, VALID_or_NOT, N, seg_Full, augment});
        }

        if (remainder != 0)
        {
            int start_Seg = total_Segs - remainder;
            int stop_Seg = total_Segs;

            threads_vec.push_back(thread{&vcf_splitter_2::seg_VALID_OR_NOT_list, this, start_Seg, stop_Seg, VALID_or_NOT, N, seg_Full, augment});
        }

        for (thread &t : threads_vec)
        {
            if (t.joinable())
            {
                t.join();
            }
        }

        threads_vec.clear();

        // sort(VALID_List.begin(), VALID_List.end());

        // for (size_t i = 0; i < VALID_List.size(); i++)
        // {
        //     cout << VALID_List[i] << endl;
        //     cout << seg_Collection[VALID_List[i]] << endl
        //          << endl;
        // }

        // exit(1);

        int valid_Count = VALID_List.size();

        cout << valid_Count << " valid segregating site(s) were found." << endl
             << endl;

        cout << "System is extracting chromosomes of segregating site(s)" << endl;

        segs_per_Thread = valid_Count / cores;
        remainder = valid_Count % cores;

        for (int core_ID = 0; core_ID < cores; core_ID++)
        {
            int start_Seg = core_ID * segs_per_Thread;
            int stop_Seg = start_Seg + segs_per_Thread;

            threads_vec.push_back(thread{&vcf_splitter_2::seg_CHR_get, this, start_Seg, stop_Seg, full_Char, chr_start_Index, chr_end_Index});
        }

        if (remainder != 0)
        {
            int start_Seg = valid_Count - remainder;
            int stop_Seg = valid_Count;

            threads_vec.push_back(thread{&vcf_splitter_2::seg_CHR_get, this, start_Seg, stop_Seg, full_Char, chr_start_Index, chr_end_Index});
        }

        for (thread &t : threads_vec)
        {
            if (t.joinable())
            {
                t.join();
            }
        }

        threads_vec.clear();

        // cout << "Round " << rounds + 1 << ": System is extracting chromosomes of segregating site(s)" << endl;

        for (size_t i = 0; i < CHR_unbound.size(); i++)
        {
            int found = 0;
            for (int CHR = 0; CHR < unique_CHRs.size(); CHR++)
            {
                if (CHR_unbound[i].first == unique_CHRs[CHR])
                {
                    CHR_collected_SEGs[CHR].push_back(CHR_unbound[i].second);
                    found = 1;
                }
            }
            if (found == 0)
            {
                unique_CHRs.push_back(CHR_unbound[i].first);
                vector<int> seg_ID;
                seg_ID.push_back(CHR_unbound[i].second);
                CHR_collected_SEGs.push_back(seg_ID);
            }
        }

        CHR_unbound.clear();

        cout << unique_CHRs.size() << " unique chromosome(s) were found." << endl
             << endl;

        for (size_t i = 0; i < unique_CHRs.size(); i++)
        {
            vector<int> collected_Segs = CHR_collected_SEGs[i];
            sort(collected_Segs.begin(), collected_Segs.end());

            int CHR_seg_size = collected_Segs.size();
            string CHR_name = unique_CHRs[i];
            cout << "Processing chromosome: " << CHR_name << "\nNumber of segregating site(s): " << CHR_seg_size << endl;

            for (int initialize = 0; initialize < CHR_seg_size; initialize++)
            {
                write_CHR.push_back("");
            }

            segs_per_Thread = CHR_seg_size / cores;
            remainder = CHR_seg_size % cores;

            for (int core_ID = 0; core_ID < cores; core_ID++)
            {
                int start_Seg = core_ID * segs_per_Thread;
                int stop_Seg = start_Seg + segs_per_Thread;

                threads_vec.push_back(thread{&vcf_splitter_2::concat_ALL, this, CHR_name, start_Seg, stop_Seg, collected_Segs, full_Char, pos_start_Index, pos_end_Index, ID_start_Index, ID_end_Index, REF_start, REF_stop, ALT_start, ALT_stop, six_8_start_Index, six_8_stop_Index});
            }

            if (remainder != 0)
            {
                int start_Seg = CHR_seg_size - remainder;
                int stop_Seg = CHR_seg_size;

                threads_vec.push_back(thread{&vcf_splitter_2::concat_ALL, this, CHR_name, start_Seg, stop_Seg, collected_Segs, full_Char, pos_start_Index, pos_end_Index, ID_start_Index, ID_end_Index, REF_start, REF_stop, ALT_start, ALT_stop, six_8_start_Index, six_8_stop_Index});
            }

            for (thread &t : threads_vec)
            {
                if (t.joinable())
                {
                    t.join();
                }
            }

            threads_vec.clear();

            // for (string lines_write : write_CHR)
            // {
            //     cout << lines_write << endl;
            // }

            string CHR_file_name = output_vcf_Folder + "/" + CHR_name + ".vcf";
            if (filesystem::exists(CHR_file_name) == 0)
            {
                function.createFile(CHR_file_name, head_Line);
            }

            fstream write_CHR_File;
            write_CHR_File.open(CHR_file_name, ios::app);

            cout << "Writing data to file: " << CHR_file_name << endl;

            for (int chr_seg = 0; chr_seg < CHR_seg_size; chr_seg++)
            {
                write_CHR_File << write_CHR[chr_seg];

                int index = collected_Segs[chr_seg];
                // string per_seg_Haps = seg_Collection[index];

                char *seg_Char;
                seg_Char = (char *)malloc((seg_Collection[index].size() + 1) * sizeof(char));
                strcpy(seg_Char, seg_Collection[index].c_str());

                int count_Size = 0;

                for (int stride = 0; stride < (seg_Collection[index].size()); stride++)
                {
                    if (count_Size == hap_Size)
                    {
                        write_CHR_File << "\t";
                        count_Size = 0;
                    }
                    write_CHR_File << seg_Char[stride];
                    count_Size++;
                }

                write_CHR_File << "\n";
                // cout << "check" << endl;

                free(seg_Char);
                // break;
                write_CHR_File.flush();
            }

            write_CHR.clear();
            write_CHR_File.close();

            if (this->summary_Individuals != 0)
            {
                cout << "Calculating individual summary" << endl;

                int found = -1;

                for (size_t chr_ind = 0; chr_ind < CHR_individuals.size(); chr_ind++)
                {
                    if (CHR_name == CHR_individuals[chr_ind])
                    {
                        found = chr_ind;
                    }
                }

                if (found == -1)
                {
                    CHR_individuals.push_back(CHR_name);
                    CHR_full_Segs.push_back(0);
                    vector<int> individuals_COUNT;

                    for (int fill = 0; fill < N; fill++)
                    {
                        individuals_COUNT.push_back(0);
                    }

                    CHR_individuals_COUNT.push_back(individuals_COUNT);

                    found = CHR_individuals.size() - 1;
                }

                CHR_full_Segs[found] = CHR_full_Segs[found] + CHR_seg_size;

                segs_per_Thread = N / cores;
                remainder = N % cores;

                for (int core_ID = 0; core_ID < cores; core_ID++)
                {
                    int start_Individual = core_ID * segs_per_Thread;
                    int stop_Individual = start_Individual + segs_per_Thread;

                    threads_vec.push_back(thread{&vcf_splitter_2::summary_Individuals_process, this, found, start_Individual, stop_Individual, sample_sequence_Tracker, collected_Segs});
                }

                if (remainder != 0)
                {
                    int start_Individual = N - remainder;
                    int stop_Individual = N;

                    threads_vec.push_back(thread{&vcf_splitter_2::summary_Individuals_process, this, found, start_Individual, stop_Individual, sample_sequence_Tracker, collected_Segs});
                }

                for (thread &t : threads_vec)
                {
                    if (t.joinable())
                    {
                        t.join();
                    }
                }

                threads_vec.clear();
            }
        }

        // exit(1);
        VALID_List.clear();
        CHR_collected_SEGs.clear();
        unique_CHRs.clear();
        seg_Collection.clear();

        // PROCESS CHR at a time.

        // for (int VL : VALID_List)
        // {
        //     cout << VL << "\t";
        // }
        // cout << endl;

        // cout << VALID_List.size() << endl;

        // first extract chr
        // then assign chr to each seg
        // if chrs have been tracked before then use them to tag and tag anything new later

        // Since the rounds are processed in order of read the lines are in read order. No need to sort the lines.

        free(site_Index);
        free(full_Char);

        cudaFree(cuda_full_Char);
        cudaFree(cuda_site_Index);

        free(VALID_or_NOT);
        cudaFree(cuda_VALID_or_NOT);

        free(chr_start_Index);
        free(chr_end_Index);
        cudaFree(cuda_chr_start_Index);
        cudaFree(cuda_chr_end_Index);

        free(pos_start_Index);
        free(pos_end_Index);
        cudaFree(cuda_pos_start_Index);
        cudaFree(cuda_pos_end_Index);

        free(ID_start_Index);
        free(ID_end_Index);
        cudaFree(cuda_ID_start_Index);
        cudaFree(cuda_ID_end_Index);

        free(REF_start);
        free(REF_stop);
        cudaFree(cuda_REF_start);
        cudaFree(cuda_REF_stop);

        free(ALT_start);
        free(ALT_stop);
        cudaFree(cuda_ALT_start);
        cudaFree(cuda_ALT_stop);

        free(six_8_start_Index);
        free(six_8_stop_Index);
        cudaFree(cuda_six_8_start_Index);
        cudaFree(cuda_six_8_stop_Index);

        free(Hap_array);
        cudaFree(cuda_Hap_array);

        cudaFree(cuda_sample_sequence_Tracker);
        if (summary_Individuals != 0)
        {
            free(sample_sequence_Tracker);
        }

        // // ! REMOVE after testing
        // cout << "DOne" << endl;
        // break;
    }

    start_stop.clear();
    concat_Segs.clear();
    all_site_Index.clear();
}

void vcf_splitter_2::summary_Individuals_process(int CHR_ID, int start_Individual, int stop_Individual, int **sample_sequence_Tracker, vector<int> collected_Segs)
{
    vector<pair<int, int>> individual_data;
    for (int col = start_Individual; col < stop_Individual; col++)
    {
        int count = 0;
        for (size_t i = 0; i < collected_Segs.size(); i++)
        {
            int row = collected_Segs[i];
            count = count + sample_sequence_Tracker[row][col];
        }
        individual_data.push_back(make_pair(col, count));
    }

    unique_lock<shared_mutex> ul(g_mutex);
    for (size_t i = 0; i < individual_data.size(); i++)
    {
        CHR_individuals_COUNT[CHR_ID][individual_data[i].first] = CHR_individuals_COUNT[CHR_ID][individual_data[i].first] + individual_data[i].second;
    }
}

void vcf_splitter_2::concat_ALL(string CHR_name, int start_Seg, int stop_Seg, vector<int> collected_Segs, char *full_Char, int *pos_start_Index, int *pos_end_Index, int *ID_start_Index, int *ID_end_Index, int *REF_start, int *REF_stop, int *ALT_start, int *ALT_stop, int *six_8_start_Index, int *six_8_stop_Index)
{
    vector<pair<int, string>> write_Line;
    for (int seg = start_Seg; seg < stop_Seg; seg++)
    {
        // cout << seg << endl;
        int index = collected_Segs[seg];

        string pos_string = "";
        for (int i = pos_start_Index[index]; i < pos_end_Index[index]; i++)
        {
            pos_string = pos_string + full_Char[i];
        }

        string ID_string = "";
        for (int i = ID_start_Index[index]; i < ID_end_Index[index]; i++)
        {
            ID_string = ID_string + full_Char[i];
        }

        string REF_string = "";
        for (int i = REF_start[index]; i < REF_stop[index]; i++)
        {
            REF_string = REF_string + full_Char[i];
        }

        string ALT_string = "";
        for (int i = ALT_start[index]; i < ALT_stop[index]; i++)
        {
            ALT_string = ALT_string + full_Char[i];
        }

        string six_8_string = "";
        for (int i = six_8_start_Index[index]; i < six_8_stop_Index[index]; i++)
        {
            six_8_string = six_8_string + full_Char[i];
        }

        write_Line.push_back(make_pair(seg, CHR_name + "\t" + pos_string + "\t" + ID_string + "\t" + REF_string + "\t" + ALT_string + "\t" + six_8_string + "\tGT\t"));
    }

    unique_lock<shared_mutex> ul(g_mutex);
    for (size_t i = 0; i < write_Line.size(); i++)
    {
        write_CHR[write_Line[i].first] = write_Line[i].second;
    }
}

void vcf_splitter_2::seg_CHR_get(int start_Seg, int stop_Seg, char *full_Char, int *CHR_start_Index, int *CHR_end_Index)
{
    vector<pair<string, int>> CHRs_collected;

    for (int seg = start_Seg; seg < stop_Seg; seg++)
    {
        int index = VALID_List[seg];
        string CHR_string = "";

        for (int i = CHR_start_Index[index]; i < CHR_end_Index[index]; i++)
        {
            /**
             * Extracting the position through concatenation of the data in the Position column of the SNP data.
             **/
            CHR_string = CHR_string + full_Char[i];
        }

        // int CHR = stoi(CHR_string);
        CHRs_collected.push_back(make_pair(CHR_string, index));
    }

    unique_lock<shared_mutex> ul(g_mutex);
    for (size_t i = 0; i < CHRs_collected.size(); i++)
    {
        CHR_unbound.push_back(make_pair(CHRs_collected[i].first, CHRs_collected[i].second));
        // unique_CHRs.push_back(CHRs_collected[i].first);
    }
}

void vcf_splitter_2::seg_VALID_OR_NOT_list(int start_Seg, int stop_Seg, int *VALID_or_NOT, int N_individuals, string seg_Full, int augment)
{
    vector<int> VALID_only;
    vector<string> segs_Haps;
    for (size_t i = start_Seg; i < stop_Seg; i++)
    {
        if (VALID_or_NOT[i] == 1)
        {
            VALID_only.push_back(i);
            int seg_array_Start = (i * ((2 * ploidy) - 1)) * N_individuals;
            segs_Haps.push_back(seg_Full.substr(seg_array_Start, augment));
        }
    }
    unique_lock<shared_mutex> ul(g_mutex);
    // for (int i : VALID_only)
    // {
    //     VALID_List.push_back(i);
    // }
    for (size_t i = 0; i < VALID_only.size(); i++)
    {
        VALID_List.push_back(VALID_only[i]);
        seg_Collection[VALID_only[i]] = segs_Haps[i];
    }
}

void vcf_splitter_2::seg_Concat(int round_ID, int start_Seg, int stop_Seg)
{
    /**
     * ! This is a multithreaded function.
     * Will spawn threads based on the umber of GPU rounds needed.
     * Will concat the segments for GPU processing per GPU rounds.
     **/

    int *site_Index;
    string Seg_sites = "";

    int total_Segs = stop_Seg - start_Seg;

    site_Index = (int *)malloc((total_Segs + 1) * sizeof(int));
    site_Index[0] = 0;

    // vector<string> total_Segregrating_sites = all_Lines;
    for (size_t i = start_Seg; i < stop_Seg; i++)
    {
        Seg_sites.append(all_Lines[i]);
        site_Index[i + 1 - start_Seg] = site_Index[i - start_Seg] + all_Lines[i].size();
        // if (round_ID == 1)
        // {
        //     cout << site_Index[i] << endl;
        // }
    }

    unique_lock<shared_mutex> ul(g_mutex);
    all_site_Index[round_ID] = site_Index;
    concat_Segs[round_ID] = Seg_sites;
}