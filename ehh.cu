#include "functions.cuh"
#include "ehh.cuh"

ehh::ehh(string range_Mode, string file_Mode_path, string fixed_Mode_value, string input_Folder, string output_Path, int cuda_ID, string intermediate_Path, int ploidy)
{
    cout << "Initiating CUDA powered Extended Haplotype Homozygosity (EHH) calculator" << endl
         << endl;

    this->range_Mode = range_Mode;
    this->file_Mode_path = file_Mode_path;

    if (range_Mode == "FIXED")
    {
        this->fixed_mode_add_OR_minus = fixed_Mode_value.at(0);
        // cout << this->fixed_mode_add_OR_minus << endl;
        if (this->fixed_mode_add_OR_minus != '+' && this->fixed_mode_add_OR_minus != '-')
        {
            cout << "ERROR IN FIXED RANGE FORMAT. TERMINATING PROGRAM." << endl;
            exit(3);
        }
        string value = fixed_Mode_value.substr(1, fixed_Mode_value.length());
        this->fixed_Mode_value = stoi(value);
    }

    this->input_Folder = input_Folder;
    this->ouput_Path = output_Path;
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

void ehh::ingress()
{
    functions function = functions();

    cout << "Range mode: " << this->range_Mode << endl;
    if (this->range_Mode == "FIXED")
    {
        cout << "Core haplotype region will be augmented by " << this->fixed_mode_add_OR_minus << this->fixed_Mode_value << " to obtain the Extended haplotype region." << endl;
        cout << "NOTE: Augmentation will be performed to the core region's start marker." << endl;
    }
    cout << endl;

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

        fstream gene_File;
        gene_File.open(file_Mode_path, ios::in);

        cout << "Processing gene list:" << endl;
        string output_File = ouput_Path + "/" +
                             country.substr(country.find_last_of("/") + 1, country.length()) + "_" +
                             filesystem::path(file_Mode_path).stem().string() +
                             ".ehh";
        string intermediate_File = intermediate_Path + "/" +
                                   country.substr(country.find_last_of("/") + 1, country.length()) + "_" +
                                   filesystem::path(file_Mode_path).stem().string() +
                                   ".log_ehh";

        cout << endl;
        cout << "Writing to file\t: " << output_File << endl;
        cout << endl;

        if (gene_File.is_open())
        {
            string gene_Combo;

            if (filesystem::exists(output_File) == 0)
            {
                function.createFile(output_File, "Gene_name\tCore_coordinates\tExtended_coordinates\tCore_Haplotype_Number\tCt\tTotal_Et\tEHH");
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

            while (getline(gene_File, gene_Combo))
            {
                vector<string> split_Data;
                function.split(split_Data, gene_Combo, "\t");
                string gene_Name = split_Data[0];
                cout << "Gene name\t: " << gene_Name << endl;

                vector<string> Core_coordinates;
                function.split(Core_coordinates, split_Data[1], ":");
                int core_start_Co = stoi(Core_coordinates[1]);
                int core_end_Co = stoi(Core_coordinates[2]);

                int ext_start_Co = 0;
                int ext_end_Co = 0;
                // ADD IF TO CHECK IF FIXED OR NOT
                if (range_Mode == "FILE")
                {
                    if (split_Data[2].at(0) == '+')
                    {
                        string value = split_Data[2].substr(1, split_Data[2].length());
                        ext_start_Co = core_start_Co;
                        ext_end_Co = ext_start_Co + stoi(value);
                    }
                    else if (split_Data[2].at(0) == '-')
                    {
                        string value = split_Data[2].substr(1, split_Data[2].length());
                        ext_start_Co = core_start_Co - stoi(value);
                        ext_end_Co = core_end_Co;
                    }
                    else
                    {
                        vector<string> Ext_coordinates;
                        function.split(Ext_coordinates, split_Data[2], ":");
                        ext_start_Co = stoi(Ext_coordinates[0]);
                        ext_end_Co = stoi(Ext_coordinates[1]);
                    }
                }
                else
                {
                    if (fixed_mode_add_OR_minus == '+')
                    {
                        ext_start_Co = core_start_Co;
                        ext_end_Co = ext_start_Co + fixed_Mode_value;
                    }
                    else
                    {
                        ext_start_Co = core_start_Co - fixed_Mode_value;
                        ext_end_Co = core_end_Co;
                    }
                }

                cout << "Coordinates:" << endl;
                cout << "Chromosome: " << Core_coordinates[0] << endl;
                cout << "Core region: "
                     << "Start: " << core_start_Co << " End: " << core_end_Co << endl;
                cout << "Extended region: "
                     << "Start: " << ext_start_Co << " End: " << ext_end_Co << endl;

                cout << endl;

                int VALID_or_NOT = 0;

                if (core_start_Co < ext_start_Co)
                {
                    VALID_or_NOT = VALID_or_NOT + 1;
                }
                if (core_end_Co > ext_end_Co)
                {
                    VALID_or_NOT = VALID_or_NOT + 2;
                }

                if (VALID_or_NOT == 0)
                {

                    vector<string> file_List;
                    cout << "System is retrieving file(s)" << endl;
                    if (folder_Index.size() > 1)
                    {
                        file_List = function.compound_interpolationSearch(folder_Index, ext_start_Co, ext_end_Co);
                    }
                    else
                    {
                        file_List.push_back(folder_Index[0].second);
                    }
                    cout << "System has retrieved all file(s)" << endl;
                    cout << endl;

                    cout << "System is collecting segregrating site(s)" << endl;

                    vector<string> collect_Segregrating_sites;
                    // core = 0 ext = 1
                    vector<int> core_OR_ext;
                    int core_Count = 0;
                    // int ext_Count = 0;

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
                                function.split_getPos_ONLY(positions, line, "\t");
                                int pos = stoi(positions[1]);

                                if (pos >= ext_start_Co && pos <= ext_end_Co)
                                {
                                    collect_Segregrating_sites.push_back(line);

                                    if (pos >= core_start_Co && pos <= core_end_Co)
                                    {
                                        core_OR_ext.push_back(0);
                                        core_Count = core_Count + 1;
                                    }
                                    else
                                    {
                                        core_OR_ext.push_back(1);
                                        // ext_Count = ext_Count + 1;
                                    }
                                }
                                else if (pos > ext_end_Co)
                                {
                                    break;
                                }
                            }
                            file.close();
                        }
                    }
                    // cout << collect_Segregrating_sites.size() << "\t" << (core_Count + ext_Count) << endl;
                    // vector<string> &total_Segregrating_sites, vector<int> core_OR_ext, int core_Count, float N
                    // int ext_Count = collect_Segregrating_sites.size() - core_Count;

                    vector<int> core_Haplotype_Collection;
                    vector<int> extended_Haplotype_Sums;
                    process_EHH(collect_Segregrating_sites, core_OR_ext, core_Count, N_float, core_Haplotype_Collection, extended_Haplotype_Sums);
                    cout << endl;

                    // Gene_name\tCore_coordinates\tExtended_coordinates\tCore_Haplotype_Number\tCt\tTotal_Et\tEHH

                    for (int core = 0; core < core_Haplotype_Collection.size(); core++)
                    {
                        cout << "Processing core haplotype " << core + 1 << ": " << endl;
                        int core_Hap_Sum = (core_Haplotype_Collection[core] * (core_Haplotype_Collection[core] - 1)) / 2;
                        cout << "Core combinations: " << core_Hap_Sum << endl;
                        int ext_Sum = extended_Haplotype_Sums[core];
                        cout << "Sum of extended haplotype combinations: " << ext_Sum << endl;

                        float EHH = (float)ext_Sum / (float)core_Hap_Sum;
                        string EHH_string = to_string(EHH);

                        if (isnan(EHH))
                        {
                            EHH_string = "DIV_0";
                        }

                        cout << "EHH: " << EHH_string << endl;

                        output << gene_Name << "\t"
                               << Core_coordinates[0] << ":" << to_string(core_start_Co) << ":" << to_string(core_end_Co)
                               << "\t" << Core_coordinates[0] << ":" << to_string(ext_start_Co) << ":" << to_string(ext_end_Co)
                               << "\t" << to_string(core + 1)
                               << "\t" << to_string(core_Hap_Sum)
                               << "\t" << to_string(ext_Sum)
                               << "\t" << EHH_string << "\n";

                        cout << endl;
                    }
                    output.flush();
                }
                else
                {
                    cout << "ERROR IN EXTENDED REGION'S COORDINATES: " << gene_Name << endl;
                    if (VALID_or_NOT == 1)
                    {
                        cout << "ERROR IN START COORDINATES. EXTENDED REGION'S START APPEARS TO FALL INSIDE THE CORE REGION. PLEASE CHECK" << endl;
                    }
                    else if (VALID_or_NOT == 2)
                    {
                        cout << "ERROR IN END COORDINATES. EXTENDED REGION'S TAIL APPEARS TO FALL INSIDE THE CORE REGION. PLEASE CHECK" << endl;
                    }
                    else if (VALID_or_NOT == 3)
                    {
                        cout << "ERROR IN BOTH START AND END COORDINATES. ENTIRE EXTENDED REGION APPEARS TO LIE INSIDE THE CORE REGION. PLEASE CHECK" << endl;
                    }
                    cout << "GENE " << gene_Name << " WILL BE SKIPPED." << endl;

                    cout << endl;
                }

                intermediate << gene_Combo << "\n";
                intermediate.flush();
            }
            output.close();
            intermediate.close();
            gene_File.close();
        }
        // REMOVE BREAK AFTER TESTING
        // break;
    }
}

__global__ void cuda_SNP_grid(int total_Segs, char *sites, int *index, char **grid)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < total_Segs)
    {
        int column = 0;
        int site_Start = index[tid];
        int site_End = index[tid + 1];

        int i = site_Start;

        while (column < 9)
        {
            if (sites[i] == '\t')
            {
                column++;
            }
            i++;
        }

        int grid_Row = tid;
        int grid_Column = 0;

        while (i < site_End)
        {
            if (sites[i] == '0' || sites[i] == '1')
            {
                grid[grid_Row][grid_Column] = sites[i];
                grid_Column++;
            }
            i++;
            // if (grid_Row == 0 && grid_Column == 1321)
            // {
            //     printf("%d\t%c\n", grid_Column, sites[i]);
            // }
        }

        tid += blockDim.x * gridDim.x;
    }
}

__global__ void cuda_Core_Haplotype_concat(int N, int total_Segs, char **grid, int *core_OR_ext, int core_Count, char *core_Hap_array, char *ext_Hap_array)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < N)
    {
        int start_Ext = tid * total_Segs;

        int start_Core = tid * core_Count;
        int stride_Core = 0;

        for (size_t stride = 0; stride < total_Segs; stride++)
        {
            // grid row = stride
            char value = grid[stride][tid];
            ext_Hap_array[start_Ext + stride] = value;

            int current_CorE = core_OR_ext[stride];

            if (current_CorE == 0)
            {
                core_Hap_array[start_Core + stride_Core] = value;
                stride_Core++;
            }
        }

        tid += blockDim.x * gridDim.x;
    }
}

__global__ void cuda_search_Extended_Haplotypes(int core_Count, char **grid, int **core_Map_grid, int *core_Sizes, int total_Segs, int *cores_Hap_Sums)
{
    // FUNCTION IS NOT USED
    // WAS REPLACED BY A CPU ALTERNATIVE WHICH PROVED TO BE FASTER DUE TO THE COMPARISON OF LONG STRINGS AND
    // SMALL CORE REQUIEMENT.

    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < core_Count)
    {
        int count = core_Sizes[tid];
        int *found = new int[count];
        int found_Count = 0;

        int sum = 0;

        for (int query = 0; query < count; query++)
        {
            int query_check_Found = 0;
            for (int i = 0; i < found_Count; i++)
            {
                if (query == found[i])
                {
                    query_check_Found = 1;
                    break;
                }
            }

            if (query_check_Found == 0)
            {
                int query_Count = 0;

                found[found_Count] = query;
                found_Count++;
                query_Count++;

                // core map  grid = columns = position of N and rows = core_Count
                int query_Col = core_Map_grid[tid][query];

                for (int subject = query + 1; subject < count; subject++)
                {

                    int subject_check_Found = 0;

                    for (int i = 0; i < found_Count; i++)
                    {
                        if (subject == found[i])
                        {
                            subject_check_Found = 1;
                            break;
                        }
                    }

                    if (subject_check_Found == 0)
                    {

                        int subject_Col = core_Map_grid[tid][subject];
                        int same_OR_not = 0;

                        for (int row = 0; row < total_Segs; row++)
                        {
                            int query_Cell = grid[row][query_Col];
                            int subject_Cell = grid[row][subject_Col];

                            if (query_Cell != subject_Cell)
                            {
                                same_OR_not = 1;
                                break;
                            }
                        }
                        if (same_OR_not == 0)
                        {
                            found[found_Count] = subject;
                            found_Count++;
                            query_Count++;
                        }
                    }
                    if (found_Count == count)
                    {
                        break;
                    }
                }
                // find sum
                // n(n-1)/2
                int sum_Partial = (query_Count * (query_Count - 1)) / 2;
                sum = sum + sum_Partial;
            }
        }
        cores_Hap_Sums[tid] = sum;
        tid += blockDim.x * gridDim.x;
    }
}

void ehh::process_EHH(vector<string> &total_Segregrating_sites, vector<int> &core_OR_ext, int core_Count, float N, vector<int> &core_Haplotype_Collection, vector<int> &extended_Haplotype_Sums)
{
    cout << "\nSystem is conducting Haplotype(s) construction" << endl;
    int num_segregrating_Sites = total_Segregrating_sites.size();

    string Seg_sites = "";
    int site_Index[num_segregrating_Sites + 1];
    site_Index[0] = 0;

    int *core_OR_ext_Array, *cuda_core_OR_ext_Array;
    core_OR_ext_Array = (int *)malloc(num_segregrating_Sites * sizeof(int));

    for (size_t i = 0; i < num_segregrating_Sites; i++)
    {
        Seg_sites.append(total_Segregrating_sites[i]);
        site_Index[i + 1] = site_Index[i] + total_Segregrating_sites[i].size();

        core_OR_ext_Array[i] = core_OR_ext[i];
    }
    char *full_Char;
    full_Char = (char *)malloc((Seg_sites.size() + 1) * sizeof(char));
    strcpy(full_Char, Seg_sites.c_str());

    total_Segregrating_sites.clear();
    core_OR_ext.clear();

    // columns = N
    // rows = seg_sites

    // int **snp_N_grid;
    char **cuda_snp_N_grid;

    // snp_N_grid = (int **)malloc(num_segregrating_Sites * sizeof(int *));
    // for (int i = 0; i < num_segregrating_Sites; i++)
    // {
    //     snp_N_grid[i] = (int *)malloc(N * sizeof(int));
    // }

    cudaMallocManaged(&cuda_snp_N_grid, N * num_segregrating_Sites * sizeof(char));
    char **tmp = (char **)malloc(num_segregrating_Sites * sizeof(tmp[0]));
    for (int i = 0; i < num_segregrating_Sites; i++)
    {
        cudaMalloc((void **)&tmp[i], N * sizeof(tmp[0][0]));
    }
    cudaMemcpy(cuda_snp_N_grid, tmp, num_segregrating_Sites * sizeof(char *), cudaMemcpyHostToDevice);
    free(tmp);

    char *cuda_full_Char;
    cudaMallocManaged(&cuda_full_Char, (Seg_sites.size() + 1) * sizeof(char));
    int *cuda_site_Index;
    cudaMallocManaged(&cuda_site_Index, (num_segregrating_Sites + 1) * sizeof(int));

    cudaMemcpy(cuda_full_Char, full_Char, (Seg_sites.size() + 1) * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(cuda_site_Index, site_Index, (num_segregrating_Sites + 1) * sizeof(int), cudaMemcpyHostToDevice);

    // cuda_SNP_grid(int total_Segs, char *sites, int *index, int **grid)
    cout << "STEP 1 OF 3: Sorting segegrating sites into an array" << endl;
    cuda_SNP_grid<<<tot_Blocks, tot_ThreadsperBlock>>>(num_segregrating_Sites, cuda_full_Char, cuda_site_Index, cuda_snp_N_grid);
    cudaDeviceSynchronize();

    cudaFree(cuda_full_Char);
    cudaFree(cuda_site_Index);

    free(full_Char);

    // for (size_t i = 0; i < num_segregrating_Sites; i++)
    // {
    //     cudaMemcpy(snp_N_grid[i], cuda_snp_N_grid[i], N * sizeof(cuda_snp_N_grid[0][0]), cudaMemcpyDeviceToHost);
    // }

    // // print grid and test
    // cout << "CHECK: " << snp_N_grid[99][1321] << "\n";
    // for (size_t c = 0; c < N; c++)
    // {
    //     for (size_t r = 0; r < num_segregrating_Sites; r++)
    //     {
    //         if (snp_N_grid[r][c] != 0 && snp_N_grid[r][c] != 1)
    //         {
    //             cout << snp_N_grid[r][c] << endl;
    //         }
    //     }
    // }
    // cout << endl;

    char *core_Hap_array, *cuda_core_Hap_array;
    char *ext_Hap_array, *cuda_ext_Hap_array;

    core_Hap_array = (char *)malloc(((N * core_Count) + 1) * sizeof(char));
    cudaMallocManaged(&cuda_core_Hap_array, ((N * core_Count) + 1) * sizeof(char));
    ext_Hap_array = (char *)malloc(((N * num_segregrating_Sites) + 1) * sizeof(char));
    cudaMallocManaged(&cuda_ext_Hap_array, ((N * num_segregrating_Sites) + 1) * sizeof(char));

    cudaMallocManaged(&cuda_core_OR_ext_Array, num_segregrating_Sites * sizeof(int));
    cudaMemcpy(cuda_core_OR_ext_Array, core_OR_ext_Array, num_segregrating_Sites * sizeof(int), cudaMemcpyHostToDevice);

    // cuda_Haplotype_core_ext(int N, int total_Segs, int **grid, int *core_OR_ext, int core_Count, int ext_Count, int *core_Hap_array, int *ext_Hap_array)
    cout << "STEP 2 OF 3: Forging core and extended haplotypes" << endl;
    cuda_Core_Haplotype_concat<<<tot_Blocks, tot_ThreadsperBlock>>>(N, num_segregrating_Sites, cuda_snp_N_grid, cuda_core_OR_ext_Array, core_Count, cuda_core_Hap_array, cuda_ext_Hap_array);
    cudaDeviceSynchronize();

    // get two Hap arrays out
    cudaMemcpy(core_Hap_array, cuda_core_Hap_array, ((N * core_Count) + 1) * sizeof(char), cudaMemcpyDeviceToHost);
    cudaMemcpy(ext_Hap_array, cuda_ext_Hap_array, ((N * num_segregrating_Sites) + 1) * sizeof(char), cudaMemcpyDeviceToHost);

    // PRINT AND TEST THE HAPLOTYPES FORMED
    // cout << "TEST" << endl;
    // cout << ext_Hap_array[0] << endl;

    string ext_Haplotypes = ext_Hap_array;
    string core_Haplotypes = core_Hap_array;

    // CORE HAPLOTYPES in VECTOR

    vector<string> ext_Haplotypes_All;
    vector<string> core_Haplotypes_All;

    // int check = 0;
    for (size_t i = 0; i < (num_segregrating_Sites * N); i = i + num_segregrating_Sites)
    {
        // cout << ext_Haplotypes.substr(i, num_segregrating_Sites) << endl;
        ext_Haplotypes_All.push_back(ext_Haplotypes.substr(i, num_segregrating_Sites));
    }

    // cout << check << endl;

    // int check = 0;
    for (size_t i = 0; i < (core_Count * N); i = i + core_Count)
    {
        // check++;
        core_Haplotypes_All.push_back(core_Haplotypes.substr(i, core_Count));
        // cout << core_Haplotypes.substr(i, core_Count) << endl;
        //   break;
    }

    // cout << core_Haplotypes_All.size() << endl;
    // cout << ext_Haplotypes_All.size() << endl;

    // FIND UNIQUE CORE HAPLOTYPES

    vector<int> Unique_haplotypes_sum_Partials;
    vector<int> Unique_Ext_haplotypes_SUM;

    // int max_Count = 0;
    // int hap_Count = 0;

    // track all found haplotype locations
    vector<int> found;

    cout << "STEP 3 OF 3: Detecting unique core haplotypes with their unique extended haplotypes" << endl;
    for (size_t query = 0; query < core_Haplotypes_All.size(); query++)
    {
        // int present_OR_not = 0;

        if (binary_search(found.begin(), found.end(), query) == false)
        {
            int query_core_Count = 0;
            query_core_Count++;
            // vector<int> locations;
            // locations.push_back(query);
            found.push_back(query);
            // hap_Count++;

            vector<pair<string, int>> ext_Haplotypes_count;
            // vector<int> ext_Haplotypes_count;

            // cout << ext_Haplotypes_All[query] << endl;
            ext_Haplotypes_count.push_back(make_pair(ext_Haplotypes_All[query], 1));
            // ext_Haplotypes.push_back(ext_Haplotypes_All[query]);
            // ext_Haplotypes_count.push_back(1);

            // track found positions
            string query_Hap = core_Haplotypes_All[query];

            for (size_t subject = query + 1; subject < core_Haplotypes_All.size(); subject++)
            {
                if (binary_search(found.begin(), found.end(), subject) == false)
                {
                    string subject_Hap = core_Haplotypes_All[subject];
                    // cout << query_Hap << "\t" << subject_Hap << endl;
                    if (query_Hap.compare(subject_Hap) == 0)
                    {
                        // cout << query_Hap << "\t" << subject_Hap << endl;
                        // cout << ext_Haplotypes_All[query] << "\t" << ext_Haplotypes_All[subject] << endl;
                        // locations.push_back(subject);
                        query_core_Count++;
                        found.push_back(subject);
                        // hap_Count++;

                        // string ext_hap_Check = ext_Haplotypes_All[subject];
                        int search_Pos = -1;
                        for (int search_Index = 0; search_Index < ext_Haplotypes_count.size(); search_Index++)
                        {

                            if (ext_Haplotypes_count[search_Index].first.compare(ext_Haplotypes_All[subject]) == 0)
                            {
                                // cout << ext_Haplotypes_count[search_Index].first << "\t" << ext_Haplotypes_All[subject] << endl;
                                search_Pos = search_Index;
                                break;
                            }
                        }

                        if (search_Pos != -1)
                        {
                            ext_Haplotypes_count[search_Pos].second = ext_Haplotypes_count[search_Pos].second + 1;
                            // cout << ext_Haplotypes_count[search_Pos].second << endl;
                        }
                        else
                        {
                            ext_Haplotypes_count.push_back(make_pair(ext_Haplotypes_All[subject], 1));
                        }
                    }
                }
            }

            Unique_haplotypes_sum_Partials.push_back(query_core_Count);

            int Ext_sum = 0;
            for (size_t unique_Haps = 0; unique_Haps < ext_Haplotypes_count.size(); unique_Haps++)
            {
                int sum_Partial = (ext_Haplotypes_count[unique_Haps].second * (ext_Haplotypes_count[unique_Haps].second - 1)) / 2;
                Ext_sum = Ext_sum + sum_Partial;
            }

            // cout << Ext_sum << endl;

            Unique_Ext_haplotypes_SUM.push_back(Ext_sum);

            // if (locations.size() > max_Count)
            // {
            //     max_Count = locations.size();
            // }
        }

        if (found.size() == core_Haplotypes_All.size())
        {
            break;
        }
        else
        {
            sort(found.begin(), found.end());
        }
    }

    cout << "             " << Unique_haplotypes_sum_Partials.size() << " unique core haplotypes were found" << endl;
    core_Haplotype_Collection = Unique_haplotypes_sum_Partials;
    extended_Haplotype_Sums = Unique_Ext_haplotypes_SUM;

    // cout << Unique_Ext_haplotypes_SUM.size() << endl;

    // cout << ext_Haplotypes_All[N - 1] << endl;

    // PROCESS EACH CORE HAPLOTYPE SEPERATELY IN THE GPU

    // CREATE 2D ARRAY
    // core map  grid = columns = position of N and rows = core_Count
    // int **core_Map_grid, **cuda_core_Map_grid;
    // int *core_Sizes, *cuda_core_Sizes;

    // core_Sizes = (int *)malloc(unique_Haplotypes_num * sizeof(int));

    // core_Map_grid = (int **)malloc(unique_Haplotypes_num * sizeof(int *));

    // for (int i = 0; i < unique_Haplotypes_num; i++)
    // {
    //     core_Map_grid[i] = (int *)malloc(max_Count * sizeof(int));
    // }
    // // cout << "1"<< endl;
    // for (size_t core_ID = 0; core_ID < unique_Haplotypes_num; core_ID++)
    // {
    //     vector<int> locations = Unique_haplotypes_locations[core_ID];
    //     core_Sizes[core_ID] = locations.size();

    //     for (size_t position = 0; position < locations.size(); position++)
    //     {
    //         core_Map_grid[core_ID][position] = locations[position];
    //     }
    // }

    // cudaMallocManaged(&cuda_core_Sizes, unique_Haplotypes_num * sizeof(int));
    // cudaMemcpy(cuda_core_Sizes, core_Sizes, unique_Haplotypes_num * sizeof(int), cudaMemcpyHostToDevice);

    // cudaMallocManaged(&cuda_core_Map_grid, unique_Haplotypes_num * max_Count * sizeof(int));
    // int **tmp_2 = (int **)malloc(unique_Haplotypes_num * sizeof(tmp_2[0]));

    // for (size_t i = 0; i < unique_Haplotypes_num; i++)
    // {
    //     cudaMalloc((void **)&tmp_2[i], max_Count * sizeof(tmp_2[0][0]));
    // }
    // cudaMemcpy(cuda_core_Map_grid, tmp_2, unique_Haplotypes_num * sizeof(int *), cudaMemcpyHostToDevice);

    // for (size_t i = 0; i < unique_Haplotypes_num; i++)
    // {
    //     cudaMemcpy(tmp_2[i], core_Map_grid[i], max_Count * sizeof(cuda_core_Map_grid[0][0]), cudaMemcpyHostToDevice);
    // }
    // free(tmp_2);

    // int *core_Hap_sums, *cuda_core_Hap_sums;
    // cudaMallocManaged(&cuda_core_Hap_sums, unique_Haplotypes_num * sizeof(int));
    // core_Hap_sums = (int *)malloc(unique_Haplotypes_num * sizeof(int));

    // // search_Extended_Haplotypes(int core_Count, char **grid, int **core_Map_grid, int *core_Sizes, int total_Segs, int *cores_Hap_Sums)
    // search_Extended_Haplotypes<<<tot_Blocks, tot_ThreadsperBlock>>>(unique_Haplotypes_num, cuda_snp_N_grid, cuda_core_Map_grid, cuda_core_Sizes, num_segregrating_Sites, cuda_core_Hap_sums);
    // cudaDeviceSynchronize();

    // cudaMemcpy(core_Hap_sums, cuda_core_Hap_sums, unique_Haplotypes_num * sizeof(int), cudaMemcpyDeviceToHost);

    free(core_Hap_array);
    free(ext_Hap_array);

    cudaFree(cuda_snp_N_grid);
    cudaFree(cuda_core_Hap_array);
    cudaFree(cuda_ext_Hap_array);
}