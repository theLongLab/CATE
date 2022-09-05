#include "functions.cuh"
#include "hap_extract.cuh"

hap_extract::hap_extract(string gene_List, string input_Folder, string output_Path, int cuda_ID, string intermediate_Path, int ploidy, string reference_File, string pop_Out)
{
    cout << "Initiating CUDA powered Haplotype extractor" << endl
         << endl;

    set_Values(gene_List, input_Folder, output_Path, cuda_ID, intermediate_Path, ploidy);

    this->reference_File = reference_File;

    transform(pop_Out.begin(), pop_Out.end(), pop_Out.begin(), ::toupper);
    if (pop_Out != "NO")
    {
        this->pop_Out = "YES";
    }
}

void hap_extract::set_Values(string gene_List, string input_Folder, string output_Path, int cuda_ID, string intermediate_Path, int ploidy)
{
    this->gene_List = gene_List;
    cout << "Gene list file path\t: " << gene_List << endl
         << endl;
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

void hap_extract::ingress()
{
    functions function = functions();

    fstream reference;
    reference.open(this->reference_File, ios::in);

    string full_Reference = "";

    if (reference.is_open())
    {
        cout << "Loading reference file: " << this->reference_File << endl;
        string line;
        // skip header
        getline(reference, line);
        while (getline(reference, line))
        {
            full_Reference.append(line);
        }
        reference.close();
    }

    transform(full_Reference.begin(), full_Reference.end(), full_Reference.begin(), ::toupper);
    // this->reference_size = full_Reference.size();

    char *reference_full;
    reference_full = (char *)malloc((full_Reference.size() + 1) * sizeof(char));
    cudaMallocManaged(&cuda_reference, (full_Reference.size() + 1) * sizeof(char));
    strcpy(reference_full, full_Reference.c_str());
    cudaMemcpy(cuda_reference, reference_full, (full_Reference.size() + 1) * sizeof(char), cudaMemcpyHostToDevice);

    free(reference_full);

    cout << "Reference file loaded" << endl
         << endl;

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

        this->N = samples * ploidy;
        cout << "Number of sequences in " << country.substr(country.find_last_of("/") + 1, country.length()) << " population [ " << samples << " x " << ploidy << " ] (N)\t: " << N << endl;

        cout << endl;

        fstream gene_File;
        gene_File.open(gene_List, ios::in);
        cout << "Processing gene list:" << endl;
        string output_File = output_Path + "/" +
                             country.substr(country.find_last_of("/") + 1, country.length()) + "_" +
                             filesystem::path(gene_List).stem().string() +
                             ".hsum";
        string intermediate_File = intermediate_Path + "/" +
                                   country.substr(country.find_last_of("/") + 1, country.length()) + "_" +
                                   filesystem::path(gene_List).stem().string() +
                                   ".log_hap";
        cout << endl;
        cout << "Writing summary to file\t: " << output_File << endl;

        string FASTA_folder = output_Path + "/" +
                              country.substr(country.find_last_of("/") + 1, country.length()) + "_" +
                              filesystem::path(gene_List).stem().string();

        if (filesystem::exists(FASTA_folder) == 0)
        {
            cout << "Creating FASTA folder: " << FASTA_folder << endl;

            filesystem::create_directory(FASTA_folder);
        }
        else
        {
            cout << "FASTA folder exists: " << FASTA_folder << endl;
        }

        cout << endl;

        if (gene_File.is_open())
        {
            string gene_Combo;

            if (filesystem::exists(output_File) == 0)
            {
                function.createFile(output_File, "Gene_name\tCoordinates\tHaplotype_number\tMutated_positions\tNumber_of_sequences");
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
                cout << "Coordinates\t: Chromosome: " << coordinates[0] << " Start: " << start_Co << " End: " << end_Co << endl;

                vector<string> collect_Segregrating_sites;
                vector<pair<int, int>> pos_INDEX;
                int count_Segs = 0;

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
                                // cout << pos << endl;
                                collect_Segregrating_sites.push_back(line);
                                pos_INDEX.push_back(make_pair(pos, count_Segs));
                                count_Segs++;
                            }
                            else if (pos > end_Co)
                            {
                                break;
                            }
                        }
                        file.close();
                    }
                }

                // GET Haps
                if (collect_Segregrating_sites.size() != 0)
                {
                    vector<string> write_Lines, write_Sequences;
                    hap_extraction(write_Lines, write_Sequences, collect_Segregrating_sites, pos_INDEX, gene_Name, coordinates[0], start_Co, end_Co);

                    string FASTA_File = FASTA_folder + "/" + gene_Name + ".fasta";
                    fstream FASTA_out;
                    FASTA_out.open(FASTA_File, ios::out);

                    if (pop_Out != "NO")
                    {
                        string FASTA_pop_File = FASTA_folder + "/" + gene_Name + "_population.fasta";
                        fstream FASTA_out_pop;
                        FASTA_out_pop.open(FASTA_pop_File, ios::out);

                        for (int hap = 0; hap < write_Lines.size(); hap++)
                        {
                            output << write_Lines[hap] + "\n";
                            FASTA_out << write_Sequences[hap] + "\n";

                            vector<string> write_Split;
                            vector<string> sequence_Split;

                            function.split(write_Split, write_Lines[hap], '\t');
                            function.split(sequence_Split, write_Sequences[hap], '\n');

                            int number_of_Sequences = stoi(write_Split[4]);
                            for (int seq = 0; seq < number_of_Sequences; seq++)
                            {
                                string seq_ID = sequence_Split[0] + "_" + to_string(seq + 1);
                                FASTA_out_pop << seq_ID + "\n" + sequence_Split[1] + "\n";
                            }
                        }
                        FASTA_out_pop.close();
                    }
                    else
                    {
                        for (int hap = 0; hap < write_Lines.size(); hap++)
                        {
                            output << write_Lines[hap] + "\n";
                            FASTA_out << write_Sequences[hap] + "\n";
                        }
                    }
                    output.flush();
                    FASTA_out.close();
                    // REMOVE break
                    // break;
                }
                cout << endl;
                intermediate << gene_Combo << "\n";
                intermediate.flush();
            }
            output.close();
            intermediate.close();
            gene_File.close();
        }

        // REMOVE AFTER TESTING
        // break;
    }
}

__global__ void cuda_hap_Forge_with_alleles(int total_Segs, char *sites, int *index, char *Hap_array, char *REF_all, char *ALT_all)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < total_Segs)
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

        REF_all[tid] = REF;
        ALT_all[tid] = ALT;

        while (column < 9)
        {
            if (sites[i] == '\t')
            {
                column++;
            }
            i++;
        }

        // int grid_Row = tid;
        // int grid_Column = 0;

        int start_Hap = tid;
        int stride = 0;

        while (i < site_End)
        {
            if (sites[i] == '0' || sites[i] == '1')
            {
                char value = sites[i];
                Hap_array[start_Hap + stride] = value;
                stride = stride + total_Segs;
                // grid[grid_Row][grid_Column] = sites[i];
                // grid_Column++;
            }
            i++;
        }

        tid += blockDim.x * gridDim.x;
    }
}

// __global__ void cuda_haplotype_Forge(int N, int total_Segs, char **grid, char *Hap_array)
// {
//     int tid = threadIdx.x + blockIdx.x * blockDim.x;

//     while (tid < N)
//     {
//         int start = tid * total_Segs;

//         for (size_t stride = 0; stride < total_Segs; stride++)
//         {
//             char value = grid[stride][tid];
//             Hap_array[start + stride] = value;
//         }

//         tid += blockDim.x * gridDim.x;
//     }
// }

__global__ void cuda_sequence_Generation(int sequence_Size, int num_of_Segs, int start, char *ref, char *haplotype, int *pos_Allele, int *index_Allele, char *REF, char *ALT, char *sequence_Full)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < sequence_Size)
    {
        int pos = tid + start;

        // check is pos exists in collection

        // binary search
        int top = 0;
        int bottom = num_of_Segs;
        int middle = top + ((bottom - top) / 2);

        char found = 'N';
        int index_Value = -1;

        while (top <= bottom)
        {
            if (pos_Allele[middle] == pos)
            {
                index_Value = index_Allele[middle];
                found = 'Y';
                break;
            }
            else if (pos_Allele[middle] < pos)
            {
                top = middle + 1;
            }
            else
            {
                bottom = middle - 1;
            }
            middle = top + ((bottom - top) / 2);
        }

        char allele = 'N';

        if (found == 'Y')
        {
            char allele_at_pos = haplotype[index_Value];

            if (allele_at_pos == '0')
            {
                allele = REF[index_Value];
            }
            else
            {
                allele = ALT[index_Value];
            }
        }
        else
        {
            // int pos_Query = 525207;
            // int adjust = pos_Query - 1;
            // cout << reference_full[adjust] << endl;
            allele = ref[pos - 1];
        }

        sequence_Full[tid] = allele;

        tid += blockDim.x * gridDim.x;
    }
}

void hap_extract::hap_extraction(vector<string> &write_Lines, vector<string> &write_Sequences, vector<string> &total_Segregrating_sites, vector<pair<int, int>> &pos_INDEX, string gene_Name, string chr, int start_Pos, int end_Pos)
{
    cout << "\nSystem is conducting Haplotype(s) construction" << endl;

    int num_segregrating_Sites = total_Segregrating_sites.size();

    string Seg_sites = "";
    int *site_Index;
    site_Index = (int *)malloc((num_segregrating_Sites + 1) * sizeof(int));
    site_Index[0] = 0;

    int *pos;
    //*cuda_pos;
    pos = (int *)malloc(num_segregrating_Sites * sizeof(int));

    vector<pair<int, int>> pos_INDEX_temp = pos_INDEX;
    sort(pos_INDEX_temp.begin(), pos_INDEX_temp.end());

    int *pos_Allele, *cuda_pos_Allele, *index_Allele, *cuda_index_Allele;
    pos_Allele = (int *)malloc(num_segregrating_Sites * sizeof(int));
    index_Allele = (int *)malloc(num_segregrating_Sites * sizeof(int));

    for (size_t i = 0; i < num_segregrating_Sites; i++)
    {
        Seg_sites.append(total_Segregrating_sites[i]);
        site_Index[i + 1] = site_Index[i] + total_Segregrating_sites[i].size();

        pos[i] = pos_INDEX[i].first;

        pos_Allele[i] = pos_INDEX_temp[i].first;
        index_Allele[i] = pos_INDEX_temp[i].second;
    }

    char *full_Char;
    full_Char = (char *)malloc((Seg_sites.size() + 1) * sizeof(char));
    strcpy(full_Char, Seg_sites.c_str());

    total_Segregrating_sites.clear();
    pos_INDEX.clear();
    pos_INDEX_temp.clear();

    // char **cuda_snp_N_grid;

    // cudaMallocManaged(&cuda_snp_N_grid, this->N * num_segregrating_Sites * sizeof(char));
    // char **tmp = (char **)malloc(num_segregrating_Sites * sizeof(tmp[0]));
    // for (int i = 0; i < num_segregrating_Sites; i++)
    // {
    //     cudaMalloc((void **)&tmp[i], this->N * sizeof(tmp[0][0]));
    // }
    // cudaMemcpy(cuda_snp_N_grid, tmp, num_segregrating_Sites * sizeof(char *), cudaMemcpyHostToDevice);

    char *cuda_full_Char;
    cudaMallocManaged(&cuda_full_Char, (Seg_sites.size() + 1) * sizeof(char));
    int *cuda_site_Index;
    cudaMallocManaged(&cuda_site_Index, (num_segregrating_Sites + 1) * sizeof(int));

    cudaMemcpy(cuda_full_Char, full_Char, (Seg_sites.size() + 1) * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(cuda_site_Index, site_Index, (num_segregrating_Sites + 1) * sizeof(int), cudaMemcpyHostToDevice);

    char *cuda_REF_char, *cuda_ALT_char;
    //*REF_char, *ALT_char,
    cudaMallocManaged(&cuda_REF_char, (num_segregrating_Sites + 1) * sizeof(char));
    cudaMallocManaged(&cuda_ALT_char, (num_segregrating_Sites + 1) * sizeof(char));

    char *Hap_array, *cuda_Hap_array;
    cudaMallocManaged(&cuda_Hap_array, ((this->N * num_segregrating_Sites) + 1) * sizeof(char));

    // ORGANIZE into array and collect Alleles
    cout << "STEP 1 OF 2: Haplotype forging from segregrating sites" << endl;
    cuda_hap_Forge_with_alleles<<<tot_Blocks, tot_ThreadsperBlock>>>(num_segregrating_Sites, cuda_full_Char, cuda_site_Index, cuda_Hap_array, cuda_REF_char, cuda_ALT_char);

    cudaError_t err = cudaGetLastError();

    if (err != cudaSuccess)
    {
        printf("CUDA Error: %s\n", cudaGetErrorString(err));

        // Possibly: exit(-1) if program cannot continue....
    }
    cudaDeviceSynchronize();
    // cout << "STEP 1 OF 2: Haplotype forging from segregrating sites" << endl;
    Hap_array = (char *)malloc(((this->N * num_segregrating_Sites) + 1) * sizeof(char));
    cudaMemcpy(Hap_array, cuda_Hap_array, ((this->N * num_segregrating_Sites) + 1) * sizeof(char), cudaMemcpyDeviceToHost);

    // REF_char = (char *)malloc((num_segregrating_Sites + 1) * sizeof(char));
    // ALT_char = (char *)malloc((num_segregrating_Sites + 1) * sizeof(char));

    // cudaMemcpy(REF_char, cuda_REF_char, (num_segregrating_Sites + 1) * sizeof(char), cudaMemcpyDeviceToHost);
    // cudaMemcpy(ALT_char, cuda_ALT_char, (num_segregrating_Sites + 1) * sizeof(char), cudaMemcpyDeviceToHost);

    cudaFree(cuda_full_Char);
    cudaFree(cuda_site_Index);
    cudaFree(cuda_Hap_array);

    free(full_Char);
    free(site_Index);

    // CONCAT

    // cout << "STEP 2 OF 3: Forging haplotypes" << endl;
    // cuda_haplotype_Forge<<<tot_Blocks, tot_ThreadsperBlock>>>(this->N, num_segregrating_Sites, cuda_snp_N_grid, cuda_Hap_array);

    // // cudaError_t err = cudaGetLastError();

    // if (err != cudaSuccess)
    // {
    //     printf("CUDA Error: %s\n", cudaGetErrorString(err));

    //     // Possibly: exit(-1) if program cannot continue....
    // }
    // cudaDeviceSynchronize();

    string haplotypes(Hap_array);
    vector<string> Haplotypes_All;

    for (int i = 0; i < (num_segregrating_Sites * this->N); i = i + num_segregrating_Sites)
    {
        // cout << ext_Haplotypes.substr(i, num_segregrating_Sites) << endl;
        Haplotypes_All.push_back(haplotypes.substr(i, num_segregrating_Sites));
    }

    // track all found haplotype locations
    vector<int> found;

    // cudaMallocManaged(&cuda_pos, (num_segregrating_Sites * sizeof(int)));
    // cudaMemcpy(cuda_pos, pos, num_segregrating_Sites * sizeof(int), cudaMemcpyHostToDevice);

    cudaMallocManaged(&cuda_pos_Allele, (num_segregrating_Sites * sizeof(int)));
    cudaMallocManaged(&cuda_index_Allele, (num_segregrating_Sites * sizeof(int)));

    cudaMemcpy(cuda_pos_Allele, pos_Allele, num_segregrating_Sites * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(cuda_index_Allele, index_Allele, num_segregrating_Sites * sizeof(int), cudaMemcpyHostToDevice);

    int sequence_Size = end_Pos - start_Pos + 1;
    int HAP_ID = 0;

    cout << "STEP 2 OF 2: Detecting unique haplotypes and synthesizing their sequences" << endl;
    for (size_t query = 0; query < Haplotypes_All.size(); query++)
    {
        if (binary_search(found.begin(), found.end(), query) == false)
        {
            HAP_ID = HAP_ID + 1;

            int query_core_Count = 0;
            query_core_Count++;

            found.push_back(query);

            string query_Hap = Haplotypes_All[query];
            // cout << query_Hap << endl;

            for (size_t subject = query + 1; subject < Haplotypes_All.size(); subject++)
            {
                if (binary_search(found.begin(), found.end(), subject) == false)
                {
                    string subject_Hap = Haplotypes_All[subject];

                    if (query_Hap.compare(subject_Hap) == 0)
                    {
                        query_core_Count++;
                        found.push_back(subject);
                        sort(found.begin(), found.end());
                    }
                }
            }

            // Process this Haplotype
            string write = gene_Name + "\t" +
                           chr + ":" + to_string(start_Pos) + ":" + to_string(end_Pos) + "\t" +
                           to_string(HAP_ID) + "\t";

            string mutation_positions = "NA";

            for (int stride = 0; stride < num_segregrating_Sites; stride++)
            {
                if (query_Hap.at(stride) == '1')
                {
                    if (mutation_positions == "NA")
                    {
                        mutation_positions = to_string(pos[stride]);
                    }
                    else
                    {
                        mutation_positions = mutation_positions + ";" + to_string(pos[stride]);
                    }
                }
            }

            write = write + mutation_positions + "\t" + to_string(query_core_Count);

            write_Lines.push_back(write);

            // Generate sequences;
            string sequence_ID = ">" + gene_Name + "_" + chr + ":" + to_string(start_Pos) + ":" + to_string(end_Pos) + "_" + to_string(HAP_ID) + "\n";

            char *sequence, *cuda_sequence;
            sequence = (char *)malloc((sequence_Size + 1) * sizeof(char));
            cudaMallocManaged(&cuda_sequence, (sequence_Size + 1) * sizeof(char));

            char *haplotye, *cuda_haplotype;
            haplotye = (char *)malloc((num_segregrating_Sites + 1) * sizeof(char));
            cudaMallocManaged(&cuda_haplotype, (num_segregrating_Sites + 1) * sizeof(char));

            strcpy(haplotye, query_Hap.c_str());
            cudaMemcpy(cuda_haplotype, haplotye, (num_segregrating_Sites + 1) * sizeof(char), cudaMemcpyHostToDevice);

            // cuda_sequence_Generation(int sequence_Size, int num_of_Segs, int start, char *ref, char *haplotype, int *pos_Allele, char *index_Allele, char *REF, char *ALT, char *sequence_Full)
            cuda_sequence_Generation<<<tot_Blocks, tot_ThreadsperBlock>>>(sequence_Size, num_segregrating_Sites, start_Pos, cuda_reference, cuda_haplotype, cuda_pos_Allele, cuda_index_Allele, cuda_REF_char, cuda_ALT_char, cuda_sequence);

            if (err != cudaSuccess)
            {
                printf("CUDA Error: %s\n", cudaGetErrorString(err));

                // Possibly: exit(-1) if program cannot continue....
            }
            cudaDeviceSynchronize();

            cudaMemcpy(sequence, cuda_sequence, (sequence_Size + 1) * sizeof(char), cudaMemcpyDeviceToHost);

            string sequence_string = sequence;
            sequence_string = sequence_string.substr(0, sequence_Size);

            string sequence_Full = sequence_ID + sequence_string;
            write_Sequences.push_back(sequence_Full);

            cudaFree(cuda_haplotype);
            cudaFree(cuda_sequence);
            free(haplotye);
            free(sequence);
        }

        if (found.size() == Haplotypes_All.size())
        {
            break;
        }
        else
        {
            sort(found.begin(), found.end());
        }
    }

    cudaFree(cuda_index_Allele);
    cudaFree(cuda_pos_Allele);
    cudaFree(cuda_REF_char);
    cudaFree(cuda_ALT_char);
    // cudaFree(cuda_snp_N_grid);

    free(index_Allele);
    free(pos_Allele);
    free(Hap_array);
    free(pos);
    // free(tmp);
}