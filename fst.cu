#include "functions.cuh"
#include "fst.cuh"
#include "fst_test_pop.cuh"

fst::fst(string gene_List, string input_Folder, string output_Path, int cuda_ID, string intermediate_Path, int ploidy, string pop_Index_path, string pop_List)
{
    cout << "Initiating CUDA powered Fst (Fixation Index) calculator" << endl
         << endl;

    this->pop_Index_path = pop_Index_path;
    this->pop_List = pop_List;

    set_Values(gene_List, input_Folder, output_Path, cuda_ID, intermediate_Path, ploidy);

    // this->gene_List = gene_List;
    // cout << "Gene list file path\t: " << gene_List << endl
    //      << endl;
    // this->input_Folder = input_Folder;
    // this->output_Path = output_Path;
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

fst::fst(string calc_Mode, int window_Size, int step_Size, string input_Folder, string output_Path, int cuda_ID, int ploidy, string pop_Index_path, string pop_List)
{
    // WINDOW MODE
    cout << "Initiating CUDA powered Fst (Fixation Index) calculator" << endl
         << endl;

    this->pop_Index_path = pop_Index_path;
    this->pop_List = pop_List;

    this->calc_Mode = "WINDOW";
    this->window_Size = window_Size;
    this->step_Size = step_Size;

    set_Values("", input_Folder, output_Path, cuda_ID, "", ploidy);
}

void fst::set_Values(string gene_List, string input_Folder, string output_Path, int cuda_ID, string intermediate_Path, int ploidy)
{
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
            cout << "ERROR STEP SIZE CANNOT BE \"0\" FOR Fst" << endl;
            exit(1);
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

void fst::ingress()
{
    vector<string> test_Pops;
    vector<vector<string>> super_Pop_per_ID_full;
    vector<vector<string>> sample_IDs;
    population_Processing(test_Pops, super_Pop_per_ID_full, sample_IDs);

    functions function = functions();

    vector<string> countries = function.get_Countries(this->input_Folder);
    cout << countries.size() << " populations were found in storage: ";
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

    // read gene line and collect lines popid wise
    // create class array for each super_pop and it creates class for each pop_ID

    // CLASS for each POP ID
    fst_test_pop *pop_ID = new fst_test_pop[test_Pops.size()];

    for (size_t i = 0; i < test_Pops.size(); i++)
    {
        pop_ID[i] = fst_test_pop(test_Pops[i], super_Pop_per_ID_full[i], sample_IDs[i], countries, this->ploidy);
        cout << endl;
        // remove BREAK: ONLY FOR TESTING
        // break;
    }

    // commit to memory the sample_ID column locations in the vcfs
    sample_location_Index(test_Pops.size(), pop_ID);
    cout << endl;

    fstream gene_File;
    gene_File.open(gene_List, ios::in);
    cout << "Processing gene list:" << endl;

    string pop_List_mod = this->pop_List;
    replace(pop_List_mod.begin(), pop_List_mod.end(), ',', '_');

    // CALC MODE SEPERATION
    if (calc_Mode != "FILE")
    {
        string output_File = output_Path + "/" + pop_List_mod + "_" + to_string(window_Size) + "_" + to_string(step_Size) + ".fst";
        // cout << output_File << endl;
        cout << endl;
        cout << "Writing to file\t: " << output_File << endl;
        cout << endl;

        int max, min;
        for (size_t i = 0; i < super_Pop_Unique_vec.size(); i++)
        {
            string folder_Path = this->input_Folder + "/" + super_Pop_Unique_vec[i];
            vector<pair<string, string>> folder_Index = function.index_Folder(folder_Path);
            if (i == 0)
            {
                max = stoi(folder_Index[folder_Index.size() - 1].first.substr(folder_Index[folder_Index.size() - 1].first.find('_') + 1));
                min = stoi(folder_Index[0].first.substr(0, folder_Index[0].first.find('_')));
            }
            else
            {
                int current_Max, current_Min;
                current_Max = stoi(folder_Index[folder_Index.size() - 1].first.substr(folder_Index[folder_Index.size() - 1].first.find('_') + 1));
                current_Min = stoi(folder_Index[0].first.substr(0, folder_Index[0].first.find('_')));

                if (current_Max > max)
                {
                    max = current_Max;
                }
                if (current_Min < min)
                {
                    min = current_Min;
                }
            }
        }

        int start_Value = min;
        int end_Value = max;

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
            function.createFile(output_File, "Coordinates\tTotal_Fst\tTotal_seg_sites\tAvg_numerator\tAvg_denominator\tAverage_Fst");
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

            for (size_t i = 0; i < test_Pops.size(); i++)
            {
                pop_ID[i].folder_Search(start_Co, end_Co);
                cout << endl;
            }

            for (size_t i = 0; i < test_Pops.size(); i++)
            {
                cout << "System is collecting segregating site(s) for " << test_Pops[i] << endl;
                if (i == 0)
                {
                    pop_ID[i].seg_Retrival(start_Co, end_Co);
                    cout << endl;
                }
                else
                {
                    // prevent searching for the same seg list

                    vector<string> query_Super_pops = super_Pop_per_ID_full[i];
                    vector<pair<string, int>> super_Pop_pop_ID;

                    for (int check = 0; check < i; check++)
                    {

                        vector<string> CHECK_Super_pops = super_Pop_per_ID_full[check];
                        vector<string> found_Pops;

                        for (string query : query_Super_pops)
                        {
                            for (string check_pop : CHECK_Super_pops)
                            {
                                if (check_pop == query)
                                {
                                    super_Pop_pop_ID.push_back(make_pair(query, check));
                                    found_Pops.push_back(query);
                                }
                            }
                        }

                        for (string remover : found_Pops)
                        {
                            query_Super_pops.erase(remove(query_Super_pops.begin(), query_Super_pops.end(), remover), query_Super_pops.end());
                        }
                    }

                    vector<string> super_Pops_FOUND;
                    vector<vector<pair<int, string>>> collect_Segregrating_sites_FOUND;

                    for (int find = 0; find < super_Pop_pop_ID.size(); find++)
                    {
                        vector<pair<int, string>> get_Segs_Found;
                        get_Segs_Found = pop_ID[super_Pop_pop_ID[find].second].return_Seg_site(super_Pop_pop_ID[find].first);
                        super_Pops_FOUND.push_back(super_Pop_pop_ID[find].first);
                        collect_Segregrating_sites_FOUND.push_back(get_Segs_Found);
                    }

                    // collect segs with found
                    pop_ID[i].seg_Retrival_with_Found(start_Co, end_Co, super_Pops_FOUND, collect_Segregrating_sites_FOUND);
                    cout << endl;
                }
            }

            // check matching pos in each and combine if present
            cout << "System is concatenating each population's segregating sites:" << endl;
            for (size_t i = 0; i < test_Pops.size(); i++)
            {
                cout << "Processing " << test_Pops[i];
                pop_ID[i].combine_Segs();
            }
            cout << "System has completed concatenating each population's segregating sites" << endl;
            cout << endl;

            // check for matching pos in all
            // make a list and use in the gpu to validate if we want to process it or not based on pos
            float Fst_Total, numerator_Avg, denominator_Avg, ratio_of_Avg;
            int Seg_count;

            process_FST(pop_ID, test_Pops.size(), Seg_count, Fst_Total, numerator_Avg, denominator_Avg, ratio_of_Avg);
            cout << endl
                 << "Fst has been calculated" << endl;

            //string FST_out = to_string(Fst_Avg);
            string ratio_of_Avg_string = to_string(ratio_of_Avg);
            string numerator_String = to_string(numerator_Avg);
            string denominator_String = to_string(denominator_Avg);

            if (Fst_Total==0)
            {
               // FST_out = "NaN";
                ratio_of_Avg_string = "NaN";
                numerator_String = "NaN";
                denominator_String = "NaN";
            }
            // if (isnan(Fst_Avg))
            // {
            //     FST_out = "NaN";
            // }
            //\tAvg_numerator\tAvg_denominator\tRatio_of_Averages
            output << to_string(start_Co) << ":" << to_string(end_Co)
                   << "\t" << Fst_Total
                   << "\t" << Seg_count
                   //<< "\t" << FST_out
                   << "\t" << numerator_String
                   << "\t" << denominator_String
                   << "\t" << ratio_of_Avg_string << "\n";

            cout << endl;

            output.flush();

            start_Co = start_Co + step_Size;
            end_Co = start_Co + window_Size;
        }
        output.close();
    }
    else
    {
        string output_File = output_Path + "/" + filesystem::path(gene_List).stem().string() + "_" + pop_List_mod + ".fst";
        // cout << output_File << endl;
        cout << endl;
        cout << "Writing to file\t: " << output_File << endl;
        cout << endl;
        string intermediate_File = intermediate_Path + "/" + filesystem::path(gene_List).stem().string() + "_" + pop_List_mod + ".log_fst";
        // cout << intermediate_File << endl;

        if (gene_File.is_open())
        {
            string gene_Combo;

            if (filesystem::exists(output_File) == 0)
            {
                // CHANGE TO FST
                function.createFile(output_File, "Gene_name\tCoordinates\tTotal_Fst\tTotal_seg_sites\tAvg_numerator\tAvg_denominator\tAverage_Fst");
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
                function.split(split_Data, gene_Combo, '\t');
                string gene_Name = split_Data[0];
                cout << "Gene name\t: " << gene_Name << endl;
                vector<string> coordinates;
                function.split(coordinates, split_Data[1], ':');
                int start_Co = stoi(coordinates[1]);
                int end_Co = stoi(coordinates[2]);
                cout << "Coordinates\t: Chromosome: " << coordinates[0] << " Start: " << start_Co << " End: " << end_Co << endl
                     << endl;

                // Process gene combo at a time.

                for (size_t i = 0; i < test_Pops.size(); i++)
                {
                    pop_ID[i].folder_Search(start_Co, end_Co);
                    cout << endl;
                }

                for (size_t i = 0; i < test_Pops.size(); i++)
                {
                    cout << "System is collecting segregating site(s) for " << test_Pops[i] << endl;
                    if (i == 0)
                    {
                        pop_ID[i].seg_Retrival(start_Co, end_Co);
                        cout << endl;
                    }
                    else
                    {
                        // prevent searching for the same seg list

                        vector<string> query_Super_pops = super_Pop_per_ID_full[i];
                        vector<pair<string, int>> super_Pop_pop_ID;

                        for (int check = 0; check < i; check++)
                        {

                            vector<string> CHECK_Super_pops = super_Pop_per_ID_full[check];
                            vector<string> found_Pops;

                            for (string query : query_Super_pops)
                            {
                                for (string check_pop : CHECK_Super_pops)
                                {
                                    if (check_pop == query)
                                    {
                                        super_Pop_pop_ID.push_back(make_pair(query, check));
                                        found_Pops.push_back(query);
                                    }
                                }
                            }

                            // remove the caught ones to prevent redundancy

                            for (string remover : found_Pops)
                            {
                                query_Super_pops.erase(remove(query_Super_pops.begin(), query_Super_pops.end(), remover), query_Super_pops.end());
                            }
                        }

                        vector<string> super_Pops_FOUND;
                        vector<vector<pair<int, string>>> collect_Segregrating_sites_FOUND;

                        for (int find = 0; find < super_Pop_pop_ID.size(); find++)
                        {
                            vector<pair<int, string>> get_Segs_Found;
                            get_Segs_Found = pop_ID[super_Pop_pop_ID[find].second].return_Seg_site(super_Pop_pop_ID[find].first);
                            super_Pops_FOUND.push_back(super_Pop_pop_ID[find].first);
                            collect_Segregrating_sites_FOUND.push_back(get_Segs_Found);
                        }

                        // collect segs with found
                        pop_ID[i].seg_Retrival_with_Found(start_Co, end_Co, super_Pops_FOUND, collect_Segregrating_sites_FOUND);
                        cout << endl;
                    }
                }

                // check matching pos in each and combine if present
                cout << "System is concatenating each population's segregating sites:" << endl;
                for (size_t i = 0; i < test_Pops.size(); i++)
                {
                    cout << "Processing " << test_Pops[i];
                    pop_ID[i].combine_Segs();
                }
                cout << "System has completed concatenating each population's segregating sites" << endl;
                cout << endl;

                // check for matching pos in all
                // make a list and use in the gpu to validate if we want to process it or not based on pos
                float Fst_Total, numerator_Avg, denominator_Avg, ratio_of_Avg;
                int Seg_count;
                // process_FST(fst_test_pop pop_IDs[], int num_Pop_Ids, int Segs_count_All, float Fst_All, float Avg_Fst)
                process_FST(pop_ID, test_Pops.size(), Seg_count, Fst_Total, numerator_Avg, denominator_Avg, ratio_of_Avg);
                cout << endl
                     << "Fst has been calculated" << endl;

                // write to file
                // Gene_name\tCoordinates\tTotal_Fst\tTotal_seg_sites\tAvg_Fst
                //string FST_out = to_string(Fst_Avg);
                string ratio_of_Avg_string = to_string(ratio_of_Avg);
                string numerator_String = to_string(numerator_Avg);
                string denominator_String = to_string(denominator_Avg);

                if (Fst_Total==0)
                {
                   // FST_out = "NaN";
                    ratio_of_Avg_string = "NaN";
                    numerator_String = "NaN";
                    denominator_String = "NaN";
                }

                // if (isnan(ratio_of_Avg))
                // {
                //     ratio_of_Avg_string = "NaN";
                // }
                output << gene_Name
                       << "\t" << coordinates[0] << ":" << to_string(start_Co) << ":" << to_string(end_Co)
                       << "\t" << Fst_Total
                       << "\t" << Seg_count
                       //<< "\t" << FST_out
                       << "\t" << numerator_String
                       << "\t" << denominator_String
                       << "\t" << ratio_of_Avg_string << "\n";

                cout << endl;

                intermediate << gene_Combo << "\n";
                output.flush();
                intermediate.flush();
                // Remove BREAK BEFORE RUN: ONLY FOR TESTING
                // break;
            }
            output.close();
            intermediate.close();
            gene_File.close();
        }
    }

    free(sample_Location_array);
    free(locations_Size);
    free(pop_seqeunce_Size_Array);

    cudaFree(cuda_sample_Location_array);
    cudaFree(cuda_locations_Size);
    cudaFree(cuda_pop_seqeunce_Size_Array);
}

void fst::sample_location_Index(int num_Pop_Ids, fst_test_pop pop_IDs[])
{
    cout << "Committing sample locations to memory" << endl;
    max_Location_Size = 0;
    vector<vector<int>> all_Sample_Locations;
    // int *locations_Size;

    // int *pop_seqeunce_Size_Array;
    //  pop_Sample_Size_Array = (int *)malloc(num_Pop_Ids * sizeof(int));
    //  cudaMallocManaged(&cuda_pop_Sample_Size_Array, num_Pop_Ids * sizeof(int));
    pop_seqeunce_Size_Array = (int *)malloc(num_Pop_Ids * sizeof(int));

    locations_Size = (int *)malloc(num_Pop_Ids * sizeof(int));

    for (int i = 0; i < num_Pop_Ids; i++)
    {
        vector<int> query_Location = pop_IDs[i].get_sample_Location();
        all_Sample_Locations.push_back(query_Location);

        // pop_Sample_Size_Array[i] = pop_IDs[i].get_Sample_Size();
        pop_seqeunce_Size_Array[i] = pop_IDs[i].get_Sequence_Size();

        locations_Size[i] = query_Location.size();

        if (locations_Size[i] > max_Location_Size)
        {
            max_Location_Size = locations_Size[i];
        }

        pop_IDs[i].clear_sample_Location();
    }

    // rows
    // now global var
    // int **sample_Location_array;
    sample_Location_array = (int **)malloc(max_Location_Size * sizeof(int *));
    // columns
    for (int i = 0; i < max_Location_Size; i++)
    {
        // PRESENT_or_NOT[i] = (int *)malloc(num_Pop_Ids * sizeof(int));
        sample_Location_array[i] = (int *)malloc((num_Pop_Ids + 1) * sizeof(int));
    }

    for (size_t c = 0; c < num_Pop_Ids; c++)
    {
        vector<int> query_Location = all_Sample_Locations[c];
        for (size_t r = 0; r < query_Location.size(); r++)
        {
            sample_Location_array[r][c] = query_Location[r];
        }
    }

    cudaMallocManaged(&cuda_pop_seqeunce_Size_Array, num_Pop_Ids * sizeof(int));
    cudaMallocManaged(&cuda_locations_Size, num_Pop_Ids * sizeof(int));
    cudaMemcpy(cuda_locations_Size, locations_Size, num_Pop_Ids * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(cuda_pop_seqeunce_Size_Array, pop_seqeunce_Size_Array, num_Pop_Ids * sizeof(int), cudaMemcpyHostToDevice);

    cudaMallocManaged(&cuda_sample_Location_array, max_Location_Size * (num_Pop_Ids + 1) * sizeof(int));
    int **tmp_3 = (int **)malloc(max_Location_Size * sizeof(tmp_3[0]));
    for (size_t i = 0; i < max_Location_Size; i++)
    {
        cudaMalloc((void **)&tmp_3[i], (num_Pop_Ids + 1) * sizeof(tmp_3[0][0]));
    }
    cudaMemcpy(cuda_sample_Location_array, tmp_3, max_Location_Size * sizeof(int *), cudaMemcpyHostToDevice);

    for (size_t i = 0; i < max_Location_Size; i++)
    {
        cudaMemcpy(tmp_3[i], sample_Location_array[i], (num_Pop_Ids + 1) * sizeof(cuda_sample_Location_array[0][0]), cudaMemcpyHostToDevice);
    }
    free(tmp_3);
}

__global__ void cuda_position_Filter(int threads_Needed, int *sites_Size_Array, int pop_ID_count, int **seg_Position_Array, int *Present_or_Not, int **found_Relationships)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    // printf("works\n");
    while (tid < threads_Needed)
    {
        int found_All = 1;

        found_Relationships[tid][0] = tid;

        int query_Pos = seg_Position_Array[tid][0];
        // printf("%d\n", query_Pos);

        for (int c = 1; c < pop_ID_count; c++)
        {
            // binary search
            int top = 0;
            int bottom = sites_Size_Array[c] - 1;
            int middle = top + ((bottom - top) / 2);
            int found = 0;

            while (top <= bottom)
            {
                if (query_Pos == seg_Position_Array[middle][c])
                {
                    found_All = found_All + 1;
                    found_Relationships[tid][c] = middle;
                    found = 1;
                    break;
                }
                else if (seg_Position_Array[middle][c] < query_Pos)
                {
                    top = middle + 1;
                }
                else
                {
                    bottom = middle - 1;
                }
                middle = top + ((bottom - top) / 2);
            }
            // if at least one is not found break
            if (found == 0)
            {
                break;
            }
        }
        // printf("%d\n", found_All);
        // printf("Hello\n");
        if (found_All == pop_ID_count)
        {
            // Fill present or not
            Present_or_Not[tid] = 1;
        }
        else
        {
            Present_or_Not[tid] = 0;
        }
        // printf("works\n");
        tid += blockDim.x * gridDim.x;
    }
}

__global__ void cuda_process_FST(int segs_Seperate, int num_Pop_Ids, int *pop_Sample_Size_Array, int *pop_seqeunce_Size_Array, int *VALID_or_NOT_ALL, int *VALID_or_NOT_FST, int *REF_Count_all, int *ALT_Count_all, float *Fst, float *CUDA_numerator, float *CUDA_denominators)
{
    // printf("works\n");
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    // printf("works\n");
    while (tid < segs_Seperate)
    {
        int start = tid * num_Pop_Ids;
        // stride
        // int end = start + num_Pop_Ids;

        int VALID = 1;

        for (int i = 0; i < num_Pop_Ids; i++)
        {
            VALID = VALID * VALID_or_NOT_ALL[start + i];
            // printf("%d\n", VALID_or_NOT_ALL[start + i]);
            if (VALID == 0)
            {
                break;
            }
        }

        if (VALID == 1)
        {
            VALID_or_NOT_FST[tid] = 1;

            float p_bar = 0;
            float q_bar = 0;
            int numerator_p_bar = 0;
            int denominator_p_bar = 0;

            float numerator_Hs = 0;
            float denominator_Hs = 0;

            for (int i = 0; i < num_Pop_Ids; i++)
            {
                float p_REF_freq = 0;
                float q_ALT_freq = 0;

                p_REF_freq = (float)REF_Count_all[start + i] / (float)pop_seqeunce_Size_Array[i];
                q_ALT_freq = (float)ALT_Count_all[start + i] / (float)pop_seqeunce_Size_Array[i];

                numerator_p_bar = numerator_p_bar + REF_Count_all[start + i];
                denominator_p_bar = denominator_p_bar + pop_seqeunce_Size_Array[i];

                float p_Square = p_REF_freq * p_REF_freq;
                float q_Square = q_ALT_freq * q_ALT_freq;

                float H_exp_Pop_Array = 1 - (p_Square + q_Square);

                numerator_Hs = numerator_Hs + (H_exp_Pop_Array * pop_Sample_Size_Array[i]);
                denominator_Hs = denominator_Hs + pop_Sample_Size_Array[i];
            }

            p_bar = (float)numerator_p_bar / (float)denominator_p_bar;
            q_bar = 1 - p_bar;

            float Hs_calc = numerator_Hs / denominator_Hs;

            float p_bar_Square = p_bar * p_bar;
            float q_bar_Square = q_bar * q_bar;
            float Ht_calc = 1 - (p_bar_Square + q_bar_Square);

            float Fst_Calc = (Ht_calc - Hs_calc) / Ht_calc;

            // Ht[tid] = Ht_calc;
            // Hs[tid] = Hs_calc;
            Fst[tid] = Fst_Calc;
            CUDA_numerator[tid] = Ht_calc - Hs_calc;
            CUDA_denominators[tid] = Ht_calc;
        }
        else
        {
            VALID_or_NOT_FST[tid] = 0;
        }
        // printf("works\n");
        tid += blockDim.x * gridDim.x;
    }
}

__global__ void cuda_process_Segs(int total_Segs, char *sites, int *index, int *pop_Sample_Size_Array, int *seg_Site_pop_ID, int **sample_Location_array, int *VALID_or_NOT, int *REF_Count_all, int *ALT_Count_all)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < total_Segs)
    {
        int pop_ID = seg_Site_pop_ID[tid];

        int column = 0;
        int site_Start = index[tid];
        int site_End = index[tid + 1];

        // printf("%d\t%d\n", site_Start, site_End);

        int i = site_Start;

        while (column < 9)
        {
            if (sites[i] == '\t')
            {
                column++;
            }
            i++;
        }

        int ALT_count = 0;
        int REF_count = 0;

        // REF = 0 ALT = 1;
        // int REF_REF_count = 0;
        // int hetero_count = 0;
        // int ALT_ALT_count = 0;

        int gen_Count = 0;
        // int site_1, site_2;

        int sample_Stride = 0;
        int pop_Sample_Size = pop_Sample_Size_Array[pop_ID];

        while (i < site_End)
        {

            if (sites[i] == '\t')
            {
                column++;
            }
            else
            {
                if (sample_Stride < pop_Sample_Size)
                {
                    if (column == sample_Location_array[sample_Stride][pop_ID])
                    {
                        if (sites[i] == '1')
                        {

                            ALT_count = ALT_count + 1;
                            gen_Count = gen_Count + 1;

                            // if (gen_Count == 1)
                            // {
                            //     site_1 = 1;
                            // }
                            // else
                            // {
                            //     site_2 = 1;
                            // }
                        }
                        else if (sites[i] == '0')
                        {
                            REF_count = REF_count + 1;
                            gen_Count = gen_Count + 1;

                            // if (gen_Count == 1)
                            // {
                            //     site_1 = 0;
                            // }
                            // else
                            // {
                            //     site_2 = 0;
                            // }
                        }

                        // TEST the CATCH
                        if (gen_Count == 2)
                        {

                            // if (site_1 == 0 && site_2 == 0)
                            // {
                            //     REF_REF_count = REF_REF_count + 1;
                            // }
                            // else if (site_1 == 1 && site_2 == 1)
                            // {
                            //     ALT_ALT_count = ALT_ALT_count + 1;
                            // }
                            // else
                            // {
                            //     hetero_count = hetero_count + 1;
                            // }

                            // check next location
                            sample_Stride = sample_Stride + 1;
                            gen_Count = 0;
                        }
                    }
                }
                else
                {
                    break;
                }
            }

            i++;
        }
        // printf("%d\t%d\n", ALT_count, REF_count);
        if (ALT_count == 0 || REF_count == 0)
        {
            VALID_or_NOT[tid] = 0;
        }
        else
        {
            VALID_or_NOT[tid] = 1;
            REF_Count_all[tid] = REF_count;
            ALT_Count_all[tid] = ALT_count;

            // genotype_REF_REF[tid] = REF_REF_count;
            // genotype_REF_ALT[tid] = hetero_count;
            // genotype_ALT_ALT[tid] = ALT_ALT_count
        }
        tid += blockDim.x * gridDim.x;
    }
}

void fst::process_FST(fst_test_pop pop_IDs[], int num_Pop_Ids, int &Segs_count_All, float &Fst_All, float &numerator_Avg, float &denominator_Avg, float &ratio_of_Avg)
{
    cout << "System is calculating Fst:" << endl
         << endl;
    cout << "STEP 1 OF 3: System is performing global merge of " << num_Pop_Ids << " population's segregrating sites" << endl;

    int max_Segs = 0;
    vector<vector<pair<int, string>>> pop_IDs_Segs;
    int *pop_Seg_size_Array, *cuda_pop_Seg_size_Array;
    pop_Seg_size_Array = (int *)malloc(num_Pop_Ids * sizeof(int));

    for (size_t i = 0; i < num_Pop_Ids; i++)
    {
        vector<pair<int, string>> final_Seg_Collection = pop_IDs[i].get_final_Seg_Collection();
        pop_IDs[i].clear_final_Seg_Collection();
        // cout << final_Seg_Collection.size() << endl;
        pop_Seg_size_Array[i] = final_Seg_Collection.size();
        pop_IDs_Segs.push_back(final_Seg_Collection);

        if (pop_Seg_size_Array[i] > max_Segs)
        {
            max_Segs = pop_Seg_size_Array[i];
        }
    }

    int first_Seg_sites = pop_Seg_size_Array[0];

    // rows
    int **seg_Positions = (int **)malloc(max_Segs * sizeof(int *));
    int *PRESENT_or_NOT = (int *)malloc(first_Seg_sites * sizeof(int));
    int **first_match_Relationships = (int **)malloc(first_Seg_sites * sizeof(int *));

    // columns
    for (int i = 0; i < max_Segs; i++)
    {
        // PRESENT_or_NOT[i] = (int *)malloc(num_Pop_Ids * sizeof(int));
        seg_Positions[i] = (int *)malloc((num_Pop_Ids + 1) * sizeof(int));
    }

    for (int i = 0; i < first_Seg_sites; i++)
    {
        first_match_Relationships[i] = (int *)malloc((num_Pop_Ids + 1) * sizeof(int));
    }

    int **cuda_seg_Positions, *cuda_PRESENT_or_NOT, **cuda_first_match_Relationships;

    // columns = num of Pop IDs
    // rows = max seg site pos
    // cudaMallocManaged(&cuda_site_Index, (num_segregrating_Sites + 1) * sizeof(int));

    cudaMallocManaged(&cuda_pop_Seg_size_Array, num_Pop_Ids * sizeof(int));

    // seg_Positions = (int *)malloc(num_Pop_Ids * max_Segs * sizeof(int));
    cudaMallocManaged(&cuda_seg_Positions, (num_Pop_Ids + 1) * max_Segs * sizeof(int));
    int **tmp_2 = (int **)malloc(max_Segs * sizeof(tmp_2[0]));
    for (size_t i = 0; i < max_Segs; i++)
    {
        cudaMalloc((void **)&tmp_2[i], (num_Pop_Ids + 1) * sizeof(tmp_2[0][0]));
    }
    cudaMemcpy(cuda_seg_Positions, tmp_2, max_Segs * sizeof(int *), cudaMemcpyHostToDevice);

    // PRESENT_or_NOT = (int *)malloc(num_Pop_Ids * max_Segs * sizeof(int));
    cudaMallocManaged(&cuda_PRESENT_or_NOT, first_Seg_sites * sizeof(int));
    // int **tmp_3 = (int **)malloc(max_Segs * sizeof(tmp_3[0]));
    // for (size_t i = 0; i < max_Segs; i++)
    // {
    //     cudaMalloc((void **)&tmp_3[i], num_Pop_Ids * sizeof(tmp_2[0][0]));
    // }
    // cudaMemcpy(cuda_PRESENT_or_NOT, tmp_3, max_Segs * sizeof(int *), cudaMemcpyHostToDevice);

    // first_match_Relationships = (int *)malloc(num_Pop_Ids * first_Seg_sites * sizeof(int));
    cudaMallocManaged(&cuda_first_match_Relationships, (num_Pop_Ids + 1) * first_Seg_sites * sizeof(int));
    // cudaMalloc((void **)&cuda_first_match_Relationships, first_Seg_sites * sizeof(int *));
    // cout << "run" << endl;
    int **tmp = (int **)malloc(first_Seg_sites * sizeof(tmp[0]));
    // cout << "run "<<first_Seg_sites << endl;
    for (int i = 0; i < first_Seg_sites; i++)
    {
        // cout << i << endl;
        cudaMalloc((void **)&tmp[i], (num_Pop_Ids + 1) * sizeof(tmp[0][0]));
    }
    // cout << "run" << endl;
    cudaMemcpy(cuda_first_match_Relationships, tmp, first_Seg_sites * sizeof(int *), cudaMemcpyHostToDevice);

    // for (int i = 0; i < first_Seg_sites; i++)
    // {
    //     cudaMalloc((void **)&cuda_first_match_Relationships[i], num_Pop_Ids * sizeof(int));
    // }

    for (int c = 0; c < num_Pop_Ids; c++)
    {
        // cout << pop_IDs_Segs.size() << endl;
        vector<pair<int, string>> final_Seg_Collection = pop_IDs_Segs[c];
        for (int r = 0; r < pop_Seg_size_Array[c]; r++)
        {
            seg_Positions[r][c] = final_Seg_Collection[r].first;
            // PRESENT_or_NOT[r][c] = 0;
        }

        // for (int r = 0; r < first_Seg_sites; r++)
        // {
        //     first_match_Relationships[r][c] = -1;
        // }
    }

    // for (size_t i = 0; i < first_Seg_sites; i++)
    // {
    //     PRESENT_or_NOT[i] = 0;
    // }

    // cudaMemcpy(cuda_full_Char, full_Char, (Seg_sites.size() + 1) * sizeof(char), cudaMemcpyHostToDevice);
    // gpu
    cudaMemcpy(cuda_pop_Seg_size_Array, pop_Seg_size_Array, num_Pop_Ids * sizeof(int), cudaMemcpyHostToDevice);

    // cudaMemcpy(cuda_seg_Positions, seg_Positions, num_Pop_Ids * max_Segs * sizeof(int), cudaMemcpyHostToDevice);
    for (size_t i = 0; i < max_Segs; i++)
    {
        cudaMemcpy(tmp_2[i], seg_Positions[i], (num_Pop_Ids + 1) * sizeof(cuda_seg_Positions[0][0]), cudaMemcpyHostToDevice);
    }

    free(tmp_2);

    // cudaMemcpy(cuda_PRESENT_or_NOT, PRESENT_or_NOT, num_Pop_Ids * max_Segs * sizeof(int), cudaMemcpyHostToDevice);
    // for (size_t i = 0; i < max_Segs; i++)
    // {
    //     cudaMemcpy(tmp_3[i], PRESENT_or_NOT[i], num_Pop_Ids * sizeof(cuda_PRESENT_or_NOT[0][0]), cudaMemcpyHostToDevice);
    // }

    // free(tmp_3);

    // for (size_t i = 0; i < first_Seg_sites; i++)
    // {
    //     cudaMemcpy(tmp[i], first_match_Relationships[i], num_Pop_Ids * sizeof(cuda_first_match_Relationships[0][0]), cudaMemcpyHostToDevice);
    // }
    free(tmp);

    // cuda_position_Filter(int *sites_Size_Array, int pop_ID_count, int *seg_Position_Array, int *Present_or_Not, int *found_Relationships);
    cuda_position_Filter<<<tot_Blocks, tot_ThreadsperBlock>>>(pop_Seg_size_Array[0], cuda_pop_Seg_size_Array, num_Pop_Ids, cuda_seg_Positions, cuda_PRESENT_or_NOT, cuda_first_match_Relationships);
    cudaDeviceSynchronize();

    cudaMemcpy(PRESENT_or_NOT, cuda_PRESENT_or_NOT, first_Seg_sites * sizeof(int), cudaMemcpyDeviceToHost);
    // cudaMemcpy(first_match_Relationships, cuda_first_match_Relationships, num_Pop_Ids * first_Seg_sites * sizeof(int), cudaMemcpyDeviceToHost);

    // for (size_t i = 0; i < max_Segs; i++)
    // {
    //     cudaMemcpy(PRESENT_or_NOT[i], cuda_PRESENT_or_NOT[i], num_Pop_Ids * sizeof(cuda_PRESENT_or_NOT[0][0]), cudaMemcpyDeviceToHost);
    // }

    for (size_t i = 0; i < first_Seg_sites; i++)
    {
        // cout << i << endl;
        cudaMemcpy(first_match_Relationships[i], cuda_first_match_Relationships[i], (num_Pop_Ids + 1) * sizeof(cuda_first_match_Relationships[0][0]), cudaMemcpyDeviceToHost);
    }

    cout << "           : System has completed performing global merge of " << num_Pop_Ids << " population's segregrating sites" << endl;

    // convert string to char for seg processing

    int num_segregrating_Sites = 0;
    for (size_t i = 0; i < first_Seg_sites; i++)
    {
        num_segregrating_Sites = num_segregrating_Sites + PRESENT_or_NOT[i];
    }

    // num of pop_Ids times the Seg sites of first pop_Id
    int tot_num_segregrating_Sites = num_segregrating_Sites * num_Pop_Ids;

    if (tot_num_segregrating_Sites > 0)
    {
        cout << "STEP 2 OF 3: System is processing and filtering all segregrating site(s)" << endl;
        int *seg_Site_pop_ID, *cuda_seg_Site_pop_ID;
        seg_Site_pop_ID = (int *)malloc(tot_num_segregrating_Sites * sizeof(int));
        cudaMallocManaged(&cuda_seg_Site_pop_ID, tot_num_segregrating_Sites * sizeof(int));

        string Seg_sites = "";
        int site_Index[tot_num_segregrating_Sites + 1];
        site_Index[0] = 0;

        // concat seg sites with pops one after the other
        int Seg_Count = 0;
        for (size_t r = 0; r < first_Seg_sites; r++)
        {
            if (PRESENT_or_NOT[r] == 1)
            {
                for (int pop = 0; pop < num_Pop_Ids; pop++)
                {
                    int array_pos = first_match_Relationships[r][pop];
                    string line = pop_IDs_Segs[pop][array_pos].second;
                    // cout << line << endl;
                    Seg_sites.append(line);
                    site_Index[Seg_Count + 1] = site_Index[Seg_Count] + line.size();
                    seg_Site_pop_ID[Seg_Count] = pop;
                    Seg_Count++;
                }
            }
        }

        free(PRESENT_or_NOT);
        cudaFree(cuda_PRESENT_or_NOT);
        cudaFree(cuda_pop_Seg_size_Array);
        free(pop_Seg_size_Array);
        cudaFree(cuda_seg_Positions);
        free(first_match_Relationships);
        cudaFree(cuda_first_match_Relationships);

        char *full_Char;
        full_Char = (char *)malloc((Seg_sites.size() + 1) * sizeof(char));
        strcpy(full_Char, Seg_sites.c_str());

        pop_IDs_Segs.clear();

        char *cuda_full_Char;
        cudaMallocManaged(&cuda_full_Char, (Seg_sites.size() + 1) * sizeof(char));
        int *cuda_site_Index;
        cudaMallocManaged(&cuda_site_Index, (tot_num_segregrating_Sites + 1) * sizeof(int));

        cudaMemcpy(cuda_full_Char, full_Char, (Seg_sites.size() + 1) * sizeof(char), cudaMemcpyHostToDevice);
        cudaMemcpy(cuda_site_Index, site_Index, (tot_num_segregrating_Sites + 1) * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(cuda_seg_Site_pop_ID, seg_Site_pop_ID, tot_num_segregrating_Sites * sizeof(int), cudaMemcpyHostToDevice);

        // int *cuda_pop_seqeunce_Size_Array, *cuda_locations_Size;
        // int **cuda_sample_Location_array;
        // cudaMallocManaged(&cuda_pop_seqeunce_Size_Array, num_Pop_Ids * sizeof(int));
        // cudaMallocManaged(&cuda_locations_Size, num_Pop_Ids * sizeof(int));

        // free(site_Index);
        // move location from here to be universally available

        // int max_Location_Size = 0;
        // vector<vector<int>> all_Sample_Locations;
        // int *locations_Size, *cuda_locations_Size;

        // int *pop_seqeunce_Size_Array, *cuda_pop_seqeunce_Size_Array;
        // // pop_Sample_Size_Array = (int *)malloc(num_Pop_Ids * sizeof(int));
        // // cudaMallocManaged(&cuda_pop_Sample_Size_Array, num_Pop_Ids * sizeof(int));
        // pop_seqeunce_Size_Array = (int *)malloc(num_Pop_Ids * sizeof(int));
        // cudaMallocManaged(&cuda_pop_seqeunce_Size_Array, num_Pop_Ids * sizeof(int));

        // locations_Size = (int *)malloc(num_Pop_Ids * sizeof(int));
        // cudaMallocManaged(&cuda_locations_Size, num_Pop_Ids * sizeof(int));

        // for (int i = 0; i < num_Pop_Ids; i++)
        // {
        //     vector<int> query_Location = pop_IDs[i].get_sample_Location();
        //     all_Sample_Locations.push_back(query_Location);

        //     // pop_Sample_Size_Array[i] = pop_IDs[i].get_Sample_Size();
        //     pop_seqeunce_Size_Array[i] = pop_IDs[i].get_Sequence_Size();

        //     locations_Size[i] = query_Location.size();

        //     if (locations_Size[i] > max_Location_Size)
        //     {
        //         max_Location_Size = locations_Size[i];
        //     }

        //     // pop_IDs[i].clear_sample_Location();
        // }

        // cudaMemcpy(cuda_locations_Size, locations_Size, num_Pop_Ids * sizeof(int), cudaMemcpyHostToDevice);
        // cudaMemcpy(cuda_pop_seqeunce_Size_Array, pop_seqeunce_Size_Array, num_Pop_Ids * sizeof(int), cudaMemcpyHostToDevice);

        // // rows
        // int **sample_Location_array, **cuda_sample_Location_array;
        // sample_Location_array = (int **)malloc(max_Location_Size * sizeof(int *));
        // // columns
        // for (int i = 0; i < max_Location_Size; i++)
        // {
        //     // PRESENT_or_NOT[i] = (int *)malloc(num_Pop_Ids * sizeof(int));
        //     sample_Location_array[i] = (int *)malloc(num_Pop_Ids * sizeof(int));
        // }

        // for (size_t c = 0; c < num_Pop_Ids; c++)
        // {
        //     vector<int> query_Location = all_Sample_Locations[c];
        //     for (size_t r = 0; r < query_Location.size(); r++)
        //     {
        //         sample_Location_array[r][c] = query_Location[r];
        //     }
        // }

        // cudaMemcpy(cuda_locations_Size, locations_Size, num_Pop_Ids * sizeof(int), cudaMemcpyHostToDevice);
        // cudaMemcpy(cuda_pop_seqeunce_Size_Array, pop_seqeunce_Size_Array, num_Pop_Ids * sizeof(int), cudaMemcpyHostToDevice);

        // // copy location array to CUDA kernal
        // cudaMallocManaged(&cuda_sample_Location_array, max_Location_Size * num_Pop_Ids * sizeof(int));
        // int **tmp_3 = (int **)malloc(max_Segs * sizeof(tmp_3[0]));
        // for (size_t i = 0; i < max_Location_Size; i++)
        // {
        //     cudaMalloc((void **)&tmp_3[i], num_Pop_Ids * sizeof(tmp_3[0][0]));
        // }
        // cudaMemcpy(cuda_sample_Location_array, tmp_3, max_Location_Size * sizeof(int *), cudaMemcpyHostToDevice);

        // for (size_t i = 0; i < max_Location_Size; i++)
        // {
        //     cudaMemcpy(tmp_3[i], sample_Location_array[i], num_Pop_Ids * sizeof(cuda_sample_Location_array[0][0]), cudaMemcpyHostToDevice);
        // }
        // free(tmp_3);

        // move location upto here to be universally available

        int *cuda_VALID_or_NOT;
        cudaMallocManaged(&cuda_VALID_or_NOT, tot_num_segregrating_Sites * sizeof(int));
        // VALID_or_NOT = (int *)malloc(tot_num_segregrating_Sites * sizeof(int));

        int *cuda_VALID_or_NOT_FST, *VALID_or_NOT_FST;
        cudaMallocManaged(&cuda_VALID_or_NOT_FST, num_segregrating_Sites * sizeof(int));
        VALID_or_NOT_FST = (int *)malloc(num_segregrating_Sites * sizeof(int));

        int *cuda_REF_Count, *cuda_ALT_Count;
        float *Fst_per_Seg, *cuda_Fst_per_Seg;

        float *numerator, *denominators;
        float *CUDA_numerator, *CUDA_denominators;

        Fst_per_Seg = (float *)malloc(num_segregrating_Sites * sizeof(float));
        numerator = (float *)malloc(num_segregrating_Sites * sizeof(float));
        denominators = (float *)malloc(num_segregrating_Sites * sizeof(float));
        // Ht_per_Seg = (float *)malloc(num_segregrating_Sites * sizeof(float));
        // Hs_per_Seg = (float *)malloc(num_segregrating_Sites * sizeof(float));

        cudaMallocManaged(&cuda_REF_Count, tot_num_segregrating_Sites * sizeof(int));
        cudaMallocManaged(&cuda_ALT_Count, tot_num_segregrating_Sites * sizeof(int));

        cudaMallocManaged(&cuda_Fst_per_Seg, num_segregrating_Sites * sizeof(float));
        cudaMallocManaged(&CUDA_numerator, num_segregrating_Sites * sizeof(float));
        cudaMallocManaged(&CUDA_denominators, num_segregrating_Sites * sizeof(float));
        // cudaMallocManaged(&cuda_Ht_per_Seg, num_segregrating_Sites * sizeof(float));
        // cudaMallocManaged(&cuda_Hs_per_Seg, num_segregrating_Sites * sizeof(float));

        // cuda_process_segs
        // cuda_process_Segs(int total_Segs, char *sites, int *index, int *pop_Sample_Size_Array, int *seg_Site_pop_ID, int **sample_Location_array, int *VALID_or_NOT, int *REF_Count_all, int *ALT_Count_all)
        cuda_process_Segs<<<tot_Blocks, tot_ThreadsperBlock>>>(tot_num_segregrating_Sites, cuda_full_Char, cuda_site_Index, cuda_locations_Size, cuda_seg_Site_pop_ID, cuda_sample_Location_array, cuda_VALID_or_NOT, cuda_REF_Count, cuda_ALT_Count);
        cudaDeviceSynchronize();
        cout << "           : System has completed processing and filtering all segregrating site(s)" << endl;
        // cout << "complete" << endl;

        // cudaMemcpy(VALID_or_NOT, cuda_VALID_or_NOT, tot_num_segregrating_Sites * sizeof(int), cudaMemcpyDeviceToHost);
        // for (size_t i = 0; i < Seg_Count; i++)
        // {
        //     cout << VALID_or_NOT[i] << endl;
        // }
        // cout << "complete" << endl;

        // cudaError_t error = cudaGetLastError();
        // if (error != cudaSuccess)
        // {
        //     // print the CUDA error message and exit
        //     printf("CUDA error: %s\n", cudaGetErrorString(error));
        //     exit(-1);
        // }

        // cout << num_segregrating_Sites << endl;
        //__global__ void cuda_process_FST(int segs_Seperate, int num_Pop_Ids, int *pop_Sample_Size_Array, int *pop_seqeunce_Size_Array, int *VALID_or_NOT_ALL, int *VALID_or_NOT_FST, int *REF_Count_all, int *ALT_Count_all, float *Fst)
        cout << "STEP 3 OF 3: System is calculating Fst per segregrating site(s)" << endl;
        cuda_process_FST<<<tot_Blocks, tot_ThreadsperBlock>>>(num_segregrating_Sites, num_Pop_Ids, cuda_locations_Size, cuda_pop_seqeunce_Size_Array, cuda_VALID_or_NOT, cuda_VALID_or_NOT_FST, cuda_REF_Count, cuda_ALT_Count, cuda_Fst_per_Seg,
                                                              CUDA_numerator, CUDA_denominators);
        cudaDeviceSynchronize();

        cudaMemcpy(VALID_or_NOT_FST, cuda_VALID_or_NOT_FST, num_segregrating_Sites * sizeof(int), cudaMemcpyDeviceToHost);
        // cudaMemcpy(Ht_per_Seg, cuda_Ht_per_Seg, num_segregrating_Sites * sizeof(float), cudaMemcpyDeviceToHost);
        // cudaMemcpy(Hs_per_Seg, cuda_Hs_per_Seg, num_segregrating_Sites * sizeof(float), cudaMemcpyDeviceToHost);
        cudaMemcpy(Fst_per_Seg, cuda_Fst_per_Seg, num_segregrating_Sites * sizeof(float), cudaMemcpyDeviceToHost);
        cudaMemcpy(numerator, CUDA_numerator, num_segregrating_Sites * sizeof(float), cudaMemcpyDeviceToHost);
        cudaMemcpy(denominators, CUDA_denominators, num_segregrating_Sites * sizeof(float), cudaMemcpyDeviceToHost);

        cout << "           : System has calculated Fst per segregrating site(s)" << endl;

        Segs_count_All = 0;
        Fst_All = 0;
        float numerator_Total = 0;
        float denominator_Total = 0;

        for (size_t i = 0; i < num_segregrating_Sites; i++)
        {
            // cout << VALID_or_NOT_FST[i] << endl;
            if (VALID_or_NOT_FST[i] == 1)
            {
                Segs_count_All = Segs_count_All + 1;
                Fst_All = Fst_All + Fst_per_Seg[i];
                numerator_Total = numerator_Total + numerator[i];
                denominator_Total = denominator_Total + denominators[i];
            }
        }

        //Avg_Fst = Fst_All / (float)Segs_count_All;
        numerator_Avg = numerator_Total / (float)Segs_count_All;
        denominator_Avg = denominator_Total / (float)Segs_count_All;
        ratio_of_Avg = numerator_Avg / denominator_Avg;
        cout << endl
             << "Total Overall Fst: " << Fst_All
             << "\nAverages Fst (Ratio of Averages): " << ratio_of_Avg << endl;

        // Seg_count is equal to seg number minus 1!

        // for (size_t r = 0; r < first_Seg_sites; r++)
        // {
        //     for (size_t c = 0; c < num_Pop_Ids; c++)
        //     {
        //         cout << first_match_Relationships[r][c] << "\t";
        //     }
        //     cout << endl;
        // }

        cudaFree(cuda_VALID_or_NOT_FST);
        cudaFree(cuda_VALID_or_NOT);
        cudaFree(cuda_Fst_per_Seg);
        cudaFree(cuda_REF_Count);
        cudaFree(cuda_ALT_Count);
        cudaFree(cuda_full_Char);
        cudaFree(cuda_site_Index);
        cudaFree(cuda_seg_Site_pop_ID);
        cudaFree(CUDA_numerator);
        cudaFree(CUDA_denominators);

        free(seg_Positions);
        free(VALID_or_NOT_FST);
        free(Fst_per_Seg);
        free(full_Char);
        free(seg_Site_pop_ID);
        free(numerator);
        free(denominators);
    }
    else
    {
        cout << "           : No segregrating sites in target region" << endl;

        free(PRESENT_or_NOT);
        cudaFree(cuda_PRESENT_or_NOT);
        cudaFree(cuda_pop_Seg_size_Array);
        free(pop_Seg_size_Array);
        cudaFree(cuda_seg_Positions);
        free(first_match_Relationships);
        cudaFree(cuda_first_match_Relationships);

        Segs_count_All = 0;
        Fst_All = 0;
        //Avg_Fst = Fst_All / (float)Segs_count_All;
    }
}

void fst::population_Processing(vector<string> &test_Pops, vector<vector<string>> &super_Pop_per_ID_full, vector<vector<string>> &sample_IDs)
{
    // make pair of sampleID and super pop
    functions function = functions();

    // vector<string> test_Pops;
    function.split(test_Pops, this->pop_List, ',');
    // vector<string> super_Pops_per_ID[test_Pops.size()];
    vector<string> super_Pops_per_ID[test_Pops.size()];

    cout << "Indexing populations:" << endl
         << endl;

    vector<string> sample_ID[test_Pops.size()];
    // vector<string> super_Pop_Unique_vec;
    //  vector<string> population_ID[test_Pops.size()];
    //  vector<string> super_Pops[test_Pops.size()];

    // vector<pair<string, string>> test_Pop_super_Pop;

    // set<string> super_Pop_Unique;
    // vector<string> super_Pop_Unique_vec;

    fstream pop_Index;
    pop_Index.open(this->pop_Index_path, ios::in);
    if (pop_Index.is_open())
    {
        cout << "Reading and mapping population index file: " << this->pop_Index_path << endl
             << endl;
        string line;

        // skip first line
        getline(pop_Index, line);

        while (getline(pop_Index, line))
        {
            vector<string> index_Split;
            function.split(index_Split, line, '\t');

            for (int i = 0; i < test_Pops.size(); i++)
            {
                string test_Pop = test_Pops[i];

                if (test_Pop == index_Split[1])
                {

                    // population_ID[i].push_back(index_Split[1]);
                    // auto it = std::remove_if(index_Split[2].begin(), index_Split[2].end(), [](char const &c)
                    //                          { return !std::isalnum(c); });

                    // index_Split[2].erase(it, index_Split[2].end());
                    if (index_Split[2].at(index_Split[2].length() - 1) == '\r')
                    {
                        // cout << "caught" << endl;
                        index_Split[2] = index_Split[2].substr(0, index_Split[2].length() - 1);
                    }
                    // sample_ID[i].push_back(index_Split[0] + "\t" + index_Split[2]);
                    sample_ID[i].push_back(index_Split[0]);
                    // super_Pops.push_back(index_Split[2]);

                    // super_Pop_Unique.insert(index_Split[2]);
                    if (std::find(super_Pop_Unique_vec.begin(), super_Pop_Unique_vec.end(), index_Split[2]) != super_Pop_Unique_vec.end())
                    {
                    }
                    else
                    {
                        super_Pop_Unique_vec.push_back(index_Split[2]);
                    }

                    if (count(super_Pops_per_ID[i].begin(), super_Pops_per_ID[i].end(), index_Split[2]))
                    {
                    }
                    else
                    {
                        super_Pops_per_ID[i].push_back(index_Split[2]);
                    }

                    // if (test_Pop_super_Pop.size() != test_Pops.size())
                    // {
                    //     int found = 0;

                    //     for (size_t i = 0; i < test_Pop_super_Pop.size(); i++)
                    //     {
                    //         if (index_Split[1] == test_Pop_super_Pop[i].first)
                    //         {
                    //             found = 1;
                    //             break;
                    //         }
                    //     }

                    //     if (found == 0)
                    //     {
                    //         test_Pop_super_Pop.push_back(make_pair(index_Split[1], index_Split[2]));
                    //     }
                    // }

                    break;
                }
            }
        }

        pop_Index.close();
    }

    cout << test_Pops.size() << " populations under study: ";
    for (size_t i = 0; i < test_Pops.size(); i++)
    {
        cout << test_Pops[i];
        if (i != test_Pops.size() - 1)
        {
            cout << ", ";
        }
    }

    cout << endl;

    cout << super_Pop_Unique_vec.size() << " super populations: ";
    // Num_super_pops = super_Pop_Unique_vec.size();

    for (int i = 0; i < super_Pop_Unique_vec.size(); i++)
    {
        cout << super_Pop_Unique_vec[i];
        if (i != super_Pop_Unique_vec.size() - 1)
        {
            cout << ", ";
        }
    }

    cout << endl
         << endl;

    cout << "Query population's super-population relationship:" << endl;

    for (size_t i = 0; i < test_Pops.size(); i++)
    {
        cout << test_Pops[i] << ": ";
        int count = 0;
        vector<string> super_pops_Collect;
        for (auto it = super_Pops_per_ID[i].begin(); it != super_Pops_per_ID[i].end(); it++)
        {
            cout << *it;
            super_pops_Collect.push_back(*it);
            if (count != super_Pops_per_ID[i].size() - 1)
            {
                cout << ", ";
            }
            count++;
        }
        super_Pop_per_ID_full.push_back(super_pops_Collect);
        cout << endl;

        vector<string> sample_Collect;
        for (auto it = sample_ID[i].begin(); it != sample_ID[i].end(); it++)
        {
            // cout << *it << endl;
            sample_Collect.push_back(*it);
        }
        sample_IDs.push_back(sample_Collect);
        // cout << test_Pop_super_Pop[i].first << "\t: " << test_Pop_super_Pop[i].second << endl;
    }

    cout << endl;

    // cout << super_Pop_per_ID_full[1][0] << endl;

    // ACCOUNT FOR POPidS SPANNING MULITPLE SUPER POPS - DONE
}
