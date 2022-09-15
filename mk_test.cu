#include "mk_test.cuh"
#include "functions.cuh"

mk_test::mk_test(string reference_Path, string alignment_Path, string gene_List, string input_Folder, string ouput_Path, int cuda_ID, string intermediate_Path, int ploidy, string genetic_Code, string start_Codons, string stop_Codons, string mode, string ORF_mode)
{
    cout << "Initiating CUDA powered McDonald–Kreitman Neutrality Index (NI) test calculator" << endl
         << endl;
    this->reference_Path = reference_Path;
    this->alignment_Path = alignment_Path;
    this->gene_List = gene_List;
    cout << "Gene list file path: " << gene_List << endl;
    this->input_Folder = input_Folder;
    this->ouput_Path = ouput_Path;
    this->intermediate_Path = intermediate_Path;
    this->ploidy = ploidy;

    this->mode = mode;

    transform(ORF_mode.begin(), ORF_mode.end(), ORF_mode.begin(), ::toupper);
    if (ORF_mode != "NO")
    {
        this->ORF_mode = "YES";
    }

    cout << "Alignment mode: " << this->mode << endl
         << "ORF's known: " << this->ORF_mode << endl
         << endl;

    this->genetic_Code = genetic_Code;
    this->start_Codons = start_Codons;
    this->stop_Codons = stop_Codons;

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

    this->primary_Intermediate_Path = this->intermediate_Path + "/" + filesystem::path(this->gene_List).stem().string();
    if (filesystem::exists(primary_Intermediate_Path) == 0)
    {
        cout << "Creating primary intermediate index folder: " << primary_Intermediate_Path << endl;
        filesystem::create_directory(primary_Intermediate_Path);
    }
    else
    {
        cout << "Primary intermediate index folder exists" << endl;
    }

    if (mode == "GENE")
    {
        string alignments = this->primary_Intermediate_Path + "/alignments";
        if (filesystem::exists(alignments) == 0)
        {
            cout << "Creating temporary alignment index folder: " << alignments << endl;
            filesystem::create_directory(alignments);
        }
        else
        {
            cout << "Temporary alignment index folder exists" << endl;
        }
    }
    cout << endl;
}

void mk_test::ingress()
{
    functions function = functions();

    vector<string> Code_split;
    print_Code(Code_split);
    cout << "Start codon(s): " << this->start_Codons << endl;
    function.split(this->start_Codons_list, this->start_Codons, ',');

    function.split(this->stop_Codons_list, this->stop_Codons, ',');
    string stop_Codon_All = "";
    for (string stop_Codon : this->stop_Codons_list)
    {
        stop_Codon_All.append(stop_Codon);
    }
    char *stop_Codons;
    stop_Codons = (char *)malloc((stop_Codon_All.size() + 1) * sizeof(char));
    strcpy(stop_Codons, stop_Codon_All.c_str());
    this->stop_Codon_size = stop_Codon_All.size();
    cudaMallocManaged(&cuda_stop_Codons, (stop_Codon_All.size() + 1) * sizeof(char));
    cudaMemcpy(cuda_stop_Codons, stop_Codons, (stop_Codon_All.size() + 1) * sizeof(char), cudaMemcpyHostToDevice);
    free(stop_Codons);
    cout << "Stop codon(s) : " << this->stop_Codons << endl;

    cout << endl;

    prepration();
    // REMOVE AFTER TESTING
    // exit(0);
    process_Genetic_code();
    process_MK();
}

void mk_test::process_Genetic_code()
{
    cout << "Indexing genetic code" << endl;
    functions function = functions();

    string string_Gen_code = "";

    vector<string> split_Aminos;
    function.split(split_Aminos, this->genetic_Code, ';');

    for (string amino : split_Aminos)
    {
        vector<string> split_Amino_codon;
        vector<string> codons;
        function.split(split_Amino_codon, amino, '|');
        function.split(codons, split_Amino_codon[1], ',');
        for (string codon : codons)
        {
            string_Gen_code.append(codon);
            string_Gen_code.append(split_Amino_codon[0]);
        }
    }

    char *index_Gen_code;
    index_Gen_code = (char *)malloc((string_Gen_code.size() + 1) * sizeof(char));
    strcpy(index_Gen_code, string_Gen_code.c_str());

    this->size_of_genetic_Code = string_Gen_code.size() + 1;
    this->index_Gen_code = index_Gen_code;

    cout << "Genetic code indexed" << endl
         << endl;
}

void mk_test::process_MK()
{
    string intermediate_Reference = this->primary_Intermediate_Path + "/" + filesystem::path(this->reference_Path).stem().string();
    if (filesystem::exists(intermediate_Reference) == 0)
    {
        cout << "ERROR: Intermediate reference index folder, " << intermediate_Reference << " has not been found at path.\n";
    }
    else
    {
        cout << "Intermediate reference index folder present: " << intermediate_Reference << endl
             << endl;
        // check if folder EXISTS
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

            fstream gene_File;
            gene_File.open(gene_List, ios::in);
            cout << "Processing gene list:" << endl;
            string output_File = ouput_Path + "/" +
                                 country.substr(country.find_last_of("/") + 1, country.length()) + "_" +
                                 filesystem::path(gene_List).stem().string() +
                                 ".mc";
            string intermediate_File = intermediate_Path + "/" +
                                       country.substr(country.find_last_of("/") + 1, country.length()) + "_" +
                                       filesystem::path(gene_List).stem().string() +
                                       ".log_mc";
            cout << endl;
            cout << "Writing to file\t: " << output_File << endl;
            cout << endl;

            if (gene_File.is_open())
            {
                string gene_Combo;

                if (filesystem::exists(output_File) == 0)
                {
                    function.createFile(output_File, "Gene_name\tGene_Coordinates\tORF_Coordinates\tDs\tDn\tPs\tPn\tNeutrality_Index");
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

                    // Codon coordinates data
                    // check if file EXISTS.

                    string write_Ds, write_Dn, write_Ps, write_Pn, write_NI;
                    vector<string> codon_Coordinates;

                    string codon_Index_File_name = intermediate_Reference + "/" + gene_Name + "_" + coordinates[0] + "_" + to_string(start_Co) + "_" + to_string(end_Co) + ".ca";
                    if (filesystem::exists(codon_Index_File_name) == 0)
                    {
                        cout << "ERROR: Codon alignment index file not found at: " << codon_Index_File_name << endl;
                        cout << "ERROR: Skipping gene: " << gene_Name << endl
                             << endl;
                        write_Ds = "NA";
                        write_Dn = "NA";
                        write_Pn = "NA";
                        write_Ps = "NA";
                        write_NI = "NA";

                        for (size_t i = 0; i < 3; i++)
                        {
                            codon_Coordinates.push_back("NA");
                        }
                    }
                    else
                    {
                        fstream get_Codon_coordinates;
                        get_Codon_coordinates.open(codon_Index_File_name, ios::in);
                        // cout << "Codon alignment index file found at " << codon_Index_File_name << endl;
                        string codon_Line_one;
                        getline(get_Codon_coordinates, codon_Line_one);
                        get_Codon_coordinates.close();
                        function.split(codon_Coordinates, codon_Line_one, '\t');
                        int codon_Start = stoi(codon_Coordinates[1]);
                        int codon_Stop = stoi(codon_Coordinates[2]);

                        cout << "Codon coordinates: Chromosome: " << coordinates[0] << " Start: " << codon_Start << " End: " << codon_Stop << endl;

                        // vector<string> collect_Segregrating_sites;
                        // vector<string> collect_Segregrating_POS;
                        vector<pair<int, string>> collect_Segregrating_site_POS;

                        vector<string> file_List;
                        cout << endl;
                        cout << "System is retrieving file(s)" << endl;
                        if (folder_Index.size() > 1)
                        {
                            file_List = function.compound_interpolationSearch(folder_Index, codon_Start, codon_Stop);
                        }
                        else
                        {
                            file_List.push_back(folder_Index[0].second);
                        }
                        cout << "System has retrieved all file(s)" << endl;
                        cout << endl;

                        cout << "System is collecting SNP site(s)" << endl;
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

                                    if (pos >= codon_Start && pos <= codon_Stop)
                                    {
                                        // collect_Segregrating_sites.push_back(line);
                                        // collect_Segregrating_POS.push_back(pos);
                                        collect_Segregrating_site_POS.push_back(make_pair(pos, line));
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
                                        //     collect_Segregrating_sites.push_back(line);
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

                        // cout << "System has collected " << num_segregrating_Sites << " segregrating site(s)" << endl;
                        // cout << endl;
                        int num_segregrating_Sites = 0;
                        // process_ORF(vector<pair<int, string>> &collect_Segregrating_site_POS, int &real_segregrating_Sites, int codon_Start, int codon_Stop, string codon_Index_File_name, int &tot_Dn, int &tot_Ds, int &tot_Pn, int &tot_Ps, float &NI)
                        int tot_Dn, tot_Ds, tot_Pn, tot_Ps;
                        float NI;
                        process_ORF(collect_Segregrating_site_POS, num_segregrating_Sites, codon_Start, codon_Stop, codon_Index_File_name, tot_Dn, tot_Ds, tot_Pn, tot_Ps, NI);
                        // calc mk syn and nonsy in cuda. Cant use MAF data cause we dont know which one is the MA
                        write_Dn = to_string(tot_Dn);
                        write_Ds = to_string(tot_Ds);
                        write_Pn = to_string(tot_Pn);
                        write_Ps = to_string(tot_Ps);
                        if (isnan(NI))
                        {
                            write_NI = "NA_DIV_0";
                        }
                        else
                        {
                            write_NI = to_string(NI);
                        }
                    }
                    //"Gene_name\tGene_Coordinates\tORF_Coordinates\tDs\tDn\tPs\tPn\tNeutrality_Index"
                    output << gene_Name << "\t"
                           << coordinates[0] << ":" << to_string(start_Co) << ":" << to_string(end_Co) << "\t"
                           << coordinates[0] << ":" << codon_Coordinates[1] << ":" << codon_Coordinates[2] << "\t"
                           << write_Ds << "\t"
                           << write_Dn << "\t"
                           << write_Ps << "\t"
                           << write_Pn << "\t"
                           << write_NI << "\n";

                    intermediate << gene_Combo << "\n";
                    output.flush();
                    intermediate.flush();
                }
                output.close();
                intermediate.close();
                gene_File.close();
            }
        }
    }
}

__global__ void cuda_process_SNPS(char *sites, int *index, int tot_Segregrating_sites, int *REF_Count_all, int *ALT_Count_all, char *REF_all, char *ALT_all)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < tot_Segregrating_sites)
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

        while (column < 7)
        {
            if (sites[i] == '\t')
            {
                column++;
            }
            i++;
        }

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

        REF_all[tid] = REF;
        ALT_all[tid] = ALT;

        REF_Count_all[tid] = REF_count;
        ALT_Count_all[tid] = ALT_count;

        tid += blockDim.x * gridDim.x;
    }
}

__global__ void cuda_process_Codons(int codon_Number, int *positions, char *REF, char *Outgroup, char *seg_REF, char *seg_ALT, int SEG_size, int *SEG_positions, int *seg_REF_count, int *seg_ALT_count, int codon_Start, int size_of_alignment_File, int genetic_Code_size, char *index_Genetic_code, int *VALID_or_NOT, int *Ds, int *Dn, int *Ps, int *Pn)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    // GET SEG SITE POSITIONS FROM PREVIOUS CUDA FUNCTION

    while (tid < codon_Number)
    {

        // start with zero and change to 1 if valid
        VALID_or_NOT[tid] = 0;
        Ds[tid] = 0;
        Dn[tid] = 0;
        Ps[tid] = 0;
        Pn[tid] = 0;

        int start_Pos = (tid * 3) + codon_Start;
        // printf("%d\n", start_Pos);

        // binary search
        char found = 'N';

        int top = 0;
        int bottom = size_of_alignment_File - 1;
        int middle = top + ((bottom - top) / 2);
        int pos_Value = -1;

        while (top <= bottom)
        {
            if (positions[middle] == start_Pos)
            {
                pos_Value = middle;
                found = 'Y';
                break;
            }
            else if (positions[middle] < start_Pos)
            {
                top = middle + 1;
            }
            else
            {
                bottom = middle - 1;
            }
            middle = top + ((bottom - top) / 2);
        }

        if (found == 'Y')
        {
            int second_Pos = start_Pos + 1;
            if (positions[pos_Value + 1] == second_Pos)
            {
                int third_Pos = start_Pos + 2;
                if (positions[pos_Value + 2] == third_Pos)
                {
                    // if (start_Pos == 206331091 || start_Pos == 206331100 || start_Pos == 206331118 || start_Pos == 206331190)
                    // {
                    //     printf("%d\n", start_Pos);
                    // }

                    // process codon if all 3 match up
                    char REF_codon_pos_1 = REF[pos_Value];
                    char REF_codon_pos_2 = REF[pos_Value + 1];
                    char REF_codon_pos_3 = REF[pos_Value + 2];

                    char Outgroup_codon_pos_1 = Outgroup[pos_Value];
                    char Outgroup_codon_pos_2 = Outgroup[pos_Value + 1];
                    char Outgroup_codon_pos_3 = Outgroup[pos_Value + 2];

                    int top_SEG = 0;
                    int bottom_SEG = SEG_size - 1;
                    int middle_SEG = top_SEG + ((bottom_SEG - top_SEG) / 2);

                    int SEG_position_location_pos_1 = -1;
                    int SEG_position_location_pos_2 = -1;
                    int SEG_position_location_pos_3 = -1;
                    char seg_Found_pos_1 = 'N';
                    char seg_Found_pos_2 = 'N';
                    char seg_Found_pos_3 = 'N';

                    // char catch_Point = 'N';

                    // first codon position

                    while (top_SEG < bottom_SEG)
                    {
                        if ((SEG_positions[middle_SEG] >= start_Pos) && (SEG_positions[middle_SEG] <= third_Pos))
                        {
                            break;
                        }
                        else if (SEG_positions[middle_SEG] < start_Pos)
                        {
                            top_SEG = middle_SEG + 1;
                        }
                        else
                        {
                            bottom_SEG = middle_SEG - 1;
                        }
                        middle_SEG = top_SEG + ((bottom_SEG - top_SEG) / 2);
                    }

                    // backward
                    for (int i = middle_SEG; i >= 0; i--)
                    {
                        if (SEG_positions[i] == start_Pos)
                        {
                            SEG_position_location_pos_1 = i;
                            seg_Found_pos_1 = 'Y';
                        }
                        else if (SEG_positions[i] == second_Pos)
                        {
                            SEG_position_location_pos_2 = i;
                            seg_Found_pos_2 = 'Y';
                        }
                        else if (SEG_positions[i] == third_Pos)
                        {
                            SEG_position_location_pos_3 = i;
                            seg_Found_pos_3 = 'Y';
                        }

                        if ((SEG_positions[i] < start_Pos) || (seg_Found_pos_1 != 'Y' && seg_Found_pos_2 != 'Y' && seg_Found_pos_3 != 'Y'))
                        {
                            break;
                        }
                    }

                    // forward
                    if (seg_Found_pos_1 != 'Y' || seg_Found_pos_2 != 'Y' || seg_Found_pos_3 != 'Y')
                    {

                        // prevernt redundancy of the search space
                        if (seg_Found_pos_1 == 'Y' && seg_Found_pos_2 == 'N')
                        {
                            middle_SEG = seg_Found_pos_1;
                        }
                        else if (seg_Found_pos_2 == 'Y')
                        {
                            middle_SEG = seg_Found_pos_2;
                        }

                        for (size_t i = middle_SEG; i < SEG_size; i++)
                        {
                            if (SEG_positions[i] == start_Pos)
                            {
                                SEG_position_location_pos_1 = i;
                                seg_Found_pos_1 = 'Y';
                            }
                            else if (SEG_positions[i] == second_Pos)
                            {
                                SEG_position_location_pos_2 = i;
                                seg_Found_pos_2 = 'Y';
                            }
                            else if (SEG_positions[i] == third_Pos)
                            {
                                SEG_position_location_pos_3 = i;
                                seg_Found_pos_3 = 'Y';
                            }

                            if ((SEG_positions[i] > third_Pos) || (seg_Found_pos_1 != 'Y' && seg_Found_pos_2 != 'Y' && seg_Found_pos_3 != 'Y'))
                            {
                                break;
                            }
                        }
                    }

                    // Use first codon position to try and find second and third
                    // if (seg_Found_pos_1 == 'Y')
                    // {
                    //     if (SEG_positions[SEG_position_location_pos_1 + 1] == second_Pos)
                    //     {
                    //         SEG_position_location_pos_2 = SEG_position_location_pos_1 + 1;
                    //         seg_Found_pos_2 = 'Y';
                    //     }

                    //     if (SEG_positions[SEG_position_location_pos_1 + 2] == third_Pos)
                    //     {
                    //         SEG_position_location_pos_3 = SEG_position_location_pos_1 + 2;
                    //         seg_Found_pos_3 = 'Y';
                    //     }
                    //     else if (SEG_positions[SEG_position_location_pos_1 + 1] == third_Pos)
                    //     {
                    //         SEG_position_location_pos_3 = SEG_position_location_pos_1 + 1;
                    //         seg_Found_pos_3 = 'Y';
                    //     }
                    // }

                    // if second is not present then go find second from scratch
                    // if (seg_Found_pos_2 == 'N')
                    // {

                    //     if (start_Pos == 206331100)
                    //     {
                    //         printf("%d\n", second_Pos);
                    //     }

                    //     top_SEG = 0;
                    //     bottom_SEG = SEG_size - 1;
                    //     // top + ((bottom - top) / 2)
                    //     middle_SEG = top_SEG + ((bottom_SEG - top_SEG) / 2);

                    //     while (top_SEG < bottom_SEG)
                    //     {
                    //         if (SEG_positions[middle_SEG] == second_Pos)
                    //         {
                    //             SEG_position_location_pos_2 = middle_SEG;
                    //             seg_Found_pos_2 = 'Y';
                    //             break;
                    //         }
                    //         else if (SEG_positions[middle_SEG] < second_Pos)
                    //         {
                    //             top_SEG = middle_SEG + 1;
                    //         }
                    //         else
                    //         {
                    //             bottom_SEG = middle_SEG - 1;
                    //         }
                    //         middle_SEG = top_SEG + ((bottom_SEG - top_SEG) / 2);
                    //     }
                    // }

                    // // if second is present use it to find third if third is not already found
                    // if ((seg_Found_pos_2 == 'Y') && (seg_Found_pos_3 == 'N'))
                    // {
                    //     if (SEG_positions[SEG_position_location_pos_2 + 1] == third_Pos)
                    //     {
                    //         SEG_position_location_pos_2 = SEG_position_location_pos_2 + 1;
                    //         seg_Found_pos_3 = 'Y';
                    //     }
                    // }

                    // // if neither second nor first can be used to find third , find third from scratch
                    // if (seg_Found_pos_3 == 'N')
                    // {
                    //     top_SEG = 0;
                    //     bottom_SEG = SEG_size - 1;
                    //     middle_SEG = top_SEG + ((bottom_SEG - top_SEG) / 2);

                    //     while (top_SEG < bottom_SEG)
                    //     {
                    //         if (SEG_positions[middle_SEG] == third_Pos)
                    //         {
                    //             SEG_position_location_pos_3 = middle_SEG;
                    //             seg_Found_pos_3 = 'Y';
                    //             break;
                    //         }
                    //         else if (SEG_positions[middle_SEG] < third_Pos)
                    //         {
                    //             top_SEG = middle_SEG + 1;
                    //         }
                    //         else
                    //         {
                    //             bottom_SEG = middle_SEG - 1;
                    //         }
                    //         middle_SEG = top_SEG + ((bottom_SEG - top_SEG) / 2);
                    //     }
                    // }

                    // get the counts needed for NI

                    // get amino acid for REF, Outgroup and Seg
                    char Seg_codon_pos_1;
                    char Seg_codon_pos_2;
                    char Seg_codon_pos_3;

                    if (seg_Found_pos_1 == 'Y')
                    {
                        if (seg_ALT_count[SEG_position_location_pos_1] == 0)
                        {
                            Seg_codon_pos_1 = REF_codon_pos_1;
                        }
                        else
                        {
                            Seg_codon_pos_1 = seg_ALT[SEG_position_location_pos_1];
                        }
                        // printf("REF: %c \t REF_VCF: %c \n", REF_codon_pos_1, seg_REF[SEG_position_location_pos_1]);
                    }
                    else
                    {
                        Seg_codon_pos_1 = REF_codon_pos_1;
                    }

                    if (seg_Found_pos_2 == 'Y')
                    {

                        if (seg_ALT_count[SEG_position_location_pos_2] == 0)
                        {
                            Seg_codon_pos_2 = REF_codon_pos_2;
                        }
                        else
                        {
                            Seg_codon_pos_2 = seg_ALT[SEG_position_location_pos_2];
                        }
                    }
                    else
                    {
                        // if (start_Pos == 206331100)
                        // {
                        //     printf("%d\n", second_Pos);
                        // }
                        Seg_codon_pos_2 = REF_codon_pos_2;
                    }

                    if (seg_Found_pos_3 == 'Y')
                    {
                        if (seg_ALT_count[SEG_position_location_pos_3] == 0)
                        {
                            Seg_codon_pos_3 = REF_codon_pos_3;
                        }
                        else
                        {
                            Seg_codon_pos_3 = seg_ALT[SEG_position_location_pos_3];
                        }
                    }
                    else
                    {
                        Seg_codon_pos_3 = REF_codon_pos_3;
                    }

                    // check for single allele mutation

                    int num_of_mutations = 0;

                    if (REF_codon_pos_1 != Outgroup_codon_pos_1 || REF_codon_pos_1 != Seg_codon_pos_1 || Outgroup_codon_pos_1 != Seg_codon_pos_1)
                    {
                        num_of_mutations = num_of_mutations + 1;
                    }

                    if (REF_codon_pos_2 != Outgroup_codon_pos_2 || REF_codon_pos_2 != Seg_codon_pos_2 || Outgroup_codon_pos_2 != Seg_codon_pos_2)
                    {
                        num_of_mutations = num_of_mutations + 1;
                    }

                    if (REF_codon_pos_3 != Outgroup_codon_pos_3 || REF_codon_pos_3 != Seg_codon_pos_3 || Outgroup_codon_pos_3 != Seg_codon_pos_3)
                    {
                        num_of_mutations = num_of_mutations + 1;
                    }

                    // process only if pos in the codon is mutatated
                    if (num_of_mutations != 0)
                    {

                        VALID_or_NOT[tid] = 1;
                        // add check to ensure that all translations are found

                        // char REF_amino_acid = '0';
                        // char Outgroup_amino_acid = '0';
                        // char Seg_amino_acid = '0';

                        // char REF_found = 'N';
                        // char Outgroup_found = 'N';
                        // char Seg_found = 'N';

                        // for (int i = 0; i < this->genetic_Code_size; i = i + 4)
                        // {
                        //     // cout << this->index_Genetic_code[i]
                        //     //      << this->index_Genetic_code[i + 1]
                        //     //      << this->index_Genetic_code[i + 2]
                        //     //      << "\t" << this->index_Genetic_code[i + 3] << "\n";

                        //     if (REF_codon_pos_1 == index_Genetic_code[i] && REF_codon_pos_2 == index_Genetic_code[i + 1] && REF_codon_pos_3 == index_Genetic_code[i + 2])
                        //     {
                        //         REF_amino_acid = index_Genetic_code[i + 3];
                        //         REF_found = 'Y';
                        //     }

                        //     if (Outgroup_codon_pos_1 == index_Genetic_code[i] && Outgroup_codon_pos_2 == index_Genetic_code[i + 1] && Outgroup_codon_pos_3 == index_Genetic_code[i + 2])
                        //     {
                        //         Outgroup_amino_acid = index_Genetic_code[i + 3];
                        //         Outgroup_found == 'Y';
                        //     }

                        //     if (Seg_codon_pos_1 == index_Genetic_code[i] && Seg_codon_pos_1 == index_Genetic_code[i + 1] && Seg_codon_pos_1 == index_Genetic_code[i + 2])
                        //     {
                        //         Seg_amino_acid = index_Genetic_code[i + 3];
                        //         Seg_found == 'Y'
                        //     }

                        //     if (REF_found = 'Y' && Outgroup_found == 'Y' && Seg_found == 'Y')
                        //     {
                        //         break;
                        //     }
                        // }

                        // count dn ds pn ps

                        // POSITION 1
                        if (REF_codon_pos_1 == Seg_codon_pos_1)
                        {
                            if (REF_codon_pos_1 != Outgroup_codon_pos_1)
                            {
                                // fixed between
                                char REF_amino_acid = '0';
                                char Outgroup_amino_acid = '0';
                                char REF_found = 'N';
                                char Outgroup_found = 'N';

                                for (int i = 0; i < genetic_Code_size; i = i + 4)
                                {
                                    if (REF_codon_pos_1 == index_Genetic_code[i] && REF_codon_pos_2 == index_Genetic_code[i + 1] && REF_codon_pos_3 == index_Genetic_code[i + 2])
                                    {
                                        REF_amino_acid = index_Genetic_code[i + 3];
                                        REF_found = 'Y';
                                    }

                                    if (Outgroup_codon_pos_1 == index_Genetic_code[i] && REF_codon_pos_2 == index_Genetic_code[i + 1] && REF_codon_pos_3 == index_Genetic_code[i + 2])
                                    {
                                        Outgroup_amino_acid = index_Genetic_code[i + 3];
                                        Outgroup_found = 'Y';
                                    }

                                    if (REF_found == 'Y' && Outgroup_found == 'Y')
                                    {
                                        break;
                                    }
                                }
                                if (REF_found == 'Y' && Outgroup_found == 'Y')
                                {
                                    if (REF_amino_acid == Outgroup_amino_acid)
                                    {
                                        // fixed synonymous
                                        Ds[tid] = Ds[tid] + 1;
                                    }
                                    else
                                    {
                                        Dn[tid] = Dn[tid] + 1;
                                    }
                                }
                                else
                                {
                                    if (REF_found == 'N')
                                    {
                                        printf("ERROR: REFERENCE CODON %c %c %c NOT FOUND\n", REF_codon_pos_1, REF_codon_pos_2, REF_codon_pos_3);
                                    }
                                    if (Outgroup_found == 'N')
                                    {
                                        printf("ERROR: OUTGROUP CODON %c %c %c NOT FOUND\n", Outgroup_codon_pos_1, REF_codon_pos_2, REF_codon_pos_3);
                                    }
                                }
                            }
                        }

                        if (REF_codon_pos_1 != Seg_codon_pos_1)
                        {
                            // polymmorphism within
                            char REF_amino_acid = '0';
                            char Seg_amino_acid = '0';
                            char REF_found = 'N';
                            char Seg_found = 'N';

                            for (int i = 0; i < genetic_Code_size; i = i + 4)
                            {
                                if (REF_codon_pos_1 == index_Genetic_code[i] && REF_codon_pos_2 == index_Genetic_code[i + 1] && REF_codon_pos_3 == index_Genetic_code[i + 2])
                                {
                                    REF_amino_acid = index_Genetic_code[i + 3];
                                    REF_found = 'Y';
                                }

                                if (Seg_codon_pos_1 == index_Genetic_code[i] && REF_codon_pos_2 == index_Genetic_code[i + 1] && REF_codon_pos_3 == index_Genetic_code[i + 2])
                                {
                                    Seg_amino_acid = index_Genetic_code[i + 3];
                                    Seg_found = 'Y';
                                }

                                if (REF_found == 'Y' && Seg_found == 'Y')
                                {
                                    break;
                                }
                            }
                            if (REF_found == 'Y' && Seg_found == 'Y')
                            {
                                // printf("%d\n", start_Pos);
                                if (REF_amino_acid == Seg_amino_acid)
                                {
                                    Ps[tid] = Ps[tid] + 1;
                                }
                                else
                                {
                                    Pn[tid] = Pn[tid] + 1;
                                }
                            }
                            else
                            {
                                if (REF_found == 'N')
                                {
                                    printf("ERROR: REFERENCE CODON %c %c %c NOT FOUND\n", REF_codon_pos_1, REF_codon_pos_2, REF_codon_pos_3);
                                }
                                if (Seg_found == 'N')
                                {
                                    printf("ERROR: OUTGROUP CODON %c %c %c NOT FOUND\n", Seg_codon_pos_1, REF_codon_pos_2, REF_codon_pos_3);
                                }
                            }
                        }

                        // POSITION 2
                        if (REF_codon_pos_2 == Seg_codon_pos_2)
                        {
                            if (REF_codon_pos_2 != Outgroup_codon_pos_2)
                            {
                                // fixed between
                                char REF_amino_acid = '0';
                                char Outgroup_amino_acid = '0';
                                char REF_found = 'N';
                                char Outgroup_found = 'N';

                                for (int i = 0; i < genetic_Code_size; i = i + 4)
                                {
                                    if (REF_codon_pos_1 == index_Genetic_code[i] && REF_codon_pos_2 == index_Genetic_code[i + 1] && REF_codon_pos_3 == index_Genetic_code[i + 2])
                                    {
                                        REF_amino_acid = index_Genetic_code[i + 3];
                                        REF_found = 'Y';
                                    }

                                    if (REF_codon_pos_1 == index_Genetic_code[i] && Outgroup_codon_pos_2 == index_Genetic_code[i + 1] && REF_codon_pos_3 == index_Genetic_code[i + 2])
                                    {
                                        Outgroup_amino_acid = index_Genetic_code[i + 3];
                                        Outgroup_found = 'Y';
                                    }

                                    if (REF_found == 'Y' && Outgroup_found == 'Y')
                                    {
                                        break;
                                    }
                                }

                                // printf("REF: %c \t Outgroup: %c \n", REF_found, Outgroup_found);
                                if (REF_found == 'Y' && Outgroup_found == 'Y')
                                {
                                    if (REF_amino_acid == Outgroup_amino_acid)
                                    {
                                        // fixed synonymous
                                        Ds[tid] = Ds[tid] + 1;
                                    }
                                    else
                                    {
                                        Dn[tid] = Dn[tid] + 1;
                                    }
                                }
                                else
                                {
                                    if (REF_found == 'N')
                                    {
                                        printf("ERROR: REFERENCE CODON %c %c %c NOT FOUND\n", REF_codon_pos_1, REF_codon_pos_2, REF_codon_pos_3);
                                    }
                                    if (Outgroup_found == 'N')
                                    {
                                        printf("ERROR: OUTGROUP CODON %c %c %c NOT FOUND\n", REF_codon_pos_1, Outgroup_codon_pos_2, REF_codon_pos_3);
                                    }
                                }
                            }
                        }

                        if (REF_codon_pos_2 != Seg_codon_pos_2)
                        {
                            // polymmorphism within
                            char REF_amino_acid = '0';
                            char Seg_amino_acid = '0';
                            char REF_found = 'N';
                            char Seg_found = 'N';

                            for (int i = 0; i < genetic_Code_size; i = i + 4)
                            {
                                if (REF_codon_pos_1 == index_Genetic_code[i] && REF_codon_pos_2 == index_Genetic_code[i + 1] && REF_codon_pos_3 == index_Genetic_code[i + 2])
                                {
                                    REF_amino_acid = index_Genetic_code[i + 3];
                                    REF_found = 'Y';
                                }

                                if (REF_codon_pos_1 == index_Genetic_code[i] && Seg_codon_pos_2 == index_Genetic_code[i + 1] && REF_codon_pos_3 == index_Genetic_code[i + 2])
                                {
                                    Seg_amino_acid = index_Genetic_code[i + 3];
                                    Seg_found = 'Y';
                                }

                                if (REF_found == 'Y' && Seg_found == 'Y')
                                {
                                    break;
                                }
                            }
                            if (REF_found == 'Y' && Seg_found == 'Y')
                            {
                                // printf("%d\n", start_Pos);
                                if (REF_amino_acid == Seg_amino_acid)
                                {
                                    Ps[tid] = Ps[tid] + 1;
                                }
                                else
                                {
                                    Pn[tid] = Pn[tid] + 1;
                                }
                            }
                            else
                            {
                                if (REF_found == 'N')
                                {
                                    printf("ERROR: REFERENCE CODON %c %c %c NOT FOUND\n", REF_codon_pos_1, REF_codon_pos_2, REF_codon_pos_3);
                                }
                                if (Seg_found == 'N')
                                {
                                    printf("ERROR: SEG CODON %c %c %c NOT FOUND\n", REF_codon_pos_1, Seg_codon_pos_2, REF_codon_pos_3);
                                }
                            }
                        }

                        // POSITION 3
                        if (REF_codon_pos_3 == Seg_codon_pos_3)
                        {
                            if (REF_codon_pos_3 != Outgroup_codon_pos_3)
                            {
                                // fixed between
                                char REF_amino_acid = '0';
                                char Outgroup_amino_acid = '0';
                                char REF_found = 'N';
                                char Outgroup_found = 'N';

                                for (int i = 0; i < genetic_Code_size; i = i + 4)
                                {
                                    if (REF_codon_pos_1 == index_Genetic_code[i] && REF_codon_pos_2 == index_Genetic_code[i + 1] && REF_codon_pos_3 == index_Genetic_code[i + 2])
                                    {
                                        REF_amino_acid = index_Genetic_code[i + 3];
                                        REF_found = 'Y';
                                    }

                                    if (REF_codon_pos_1 == index_Genetic_code[i] && REF_codon_pos_2 == index_Genetic_code[i + 1] && Outgroup_codon_pos_3 == index_Genetic_code[i + 2])
                                    {
                                        Outgroup_amino_acid = index_Genetic_code[i + 3];
                                        Outgroup_found = 'Y';
                                    }

                                    if (REF_found == 'Y' && Outgroup_found == 'Y')
                                    {
                                        break;
                                    }
                                }

                                // printf(" %c %c %c \n", REF_codon_pos_1, REF_codon_pos_2,REF_codon_pos_3);
                                // printf(" %c %c %c : %c \n", index_Genetic_code[0], index_Genetic_code[1],index_Genetic_code[2],index_Genetic_code[3]);
                                // printf("REF = %c \t Out = %c \n", REF_amino_acid, Outgroup_amino_acid);
                                // printf("ALT = %c\n",Outgroup_amino_acid);
                                if (REF_found == 'Y' && Outgroup_found == 'Y')
                                {
                                    if (REF_amino_acid == Outgroup_amino_acid)
                                    {
                                        // fixed synonymous
                                        // printf("REF = %c \t Out = %c \n", REF_amino_acid, Outgroup_amino_acid);
                                        Ds[tid] = Ds[tid] + 1;
                                    }
                                    else
                                    {
                                        Dn[tid] = Dn[tid] + 1;
                                    }
                                }
                                else
                                {
                                    if (REF_found == 'N')
                                    {
                                        printf("ERROR: REFERENCE CODON %c %c %c NOT FOUND\n", REF_codon_pos_1, REF_codon_pos_2, REF_codon_pos_3);
                                    }
                                    if (Outgroup_found == 'N')
                                    {
                                        printf("ERROR: OUTGROUP CODON %c %c %c NOT FOUND\n", REF_codon_pos_1, REF_codon_pos_2, Outgroup_codon_pos_3);
                                    }
                                }
                            }
                        }

                        if (REF_codon_pos_3 != Seg_codon_pos_3)
                        {
                            // polymmorphism within
                            char REF_amino_acid = '0';
                            char Seg_amino_acid = '0';
                            char REF_found = 'N';
                            char Seg_found = 'N';

                            for (int i = 0; i < genetic_Code_size; i = i + 4)
                            {
                                if (REF_codon_pos_1 == index_Genetic_code[i] && REF_codon_pos_2 == index_Genetic_code[i + 1] && REF_codon_pos_3 == index_Genetic_code[i + 2])
                                {
                                    REF_amino_acid = index_Genetic_code[i + 3];
                                    REF_found = 'Y';
                                }

                                if (REF_codon_pos_1 == index_Genetic_code[i] && REF_codon_pos_2 == index_Genetic_code[i + 1] && Seg_codon_pos_3 == index_Genetic_code[i + 2])
                                {
                                    Seg_amino_acid = index_Genetic_code[i + 3];
                                    Seg_found = 'Y';
                                }

                                if (REF_found == 'Y' && Seg_found == 'Y')
                                {
                                    break;
                                }
                            }

                            // check if both are found
                            if (REF_found == 'Y' && Seg_found == 'Y')
                            {
                                // printf("%d\n", start_Pos);
                                //  printf("FOUND");
                                if (REF_amino_acid == Seg_amino_acid)
                                {
                                    Ps[tid] = Ps[tid] + 1;
                                }
                                else
                                {
                                    Pn[tid] = Pn[tid] + 1;
                                }
                            }
                            else
                            {
                                if (REF_found == 'N')
                                {
                                    printf("ERROR: REFERENCE CODON %c %c %c NOT FOUND\n", REF_codon_pos_1, REF_codon_pos_2, REF_codon_pos_3);
                                }
                                if (Seg_found == 'N')
                                {
                                    printf("ERROR: SEG CODON %c %c %c NOT FOUND\n", REF_codon_pos_1, REF_codon_pos_2, Seg_codon_pos_3);
                                }
                            }
                        }

                        // all 3 dn ds over
                    }
                }
            }
        }
        tid += blockDim.x * gridDim.x;
    }
}

void mk_test::process_ORF(vector<pair<int, string>> &collect_Segregrating_site_POS, int &real_segregrating_Sites, int codon_Start, int codon_Stop, string codon_Index_File_name, int &tot_Dn, int &tot_Ds, int &tot_Pn, int &tot_Ps, float &NI)
{
    cout << "\nCalculating Neutrality Index (NI): " << endl;

    int num_segregrating_Sites = collect_Segregrating_site_POS.size();
    cout << endl
         << "STEP 1 OF 2: Processing " << num_segregrating_Sites << " SNP sites " << endl;
    sort(collect_Segregrating_site_POS.begin(), collect_Segregrating_site_POS.end());
    // map codons from genetic code;
    string Seg_sites = "";
    int site_Index[num_segregrating_Sites + 1];
    site_Index[0] = 0;

    for (size_t i = 0; i < num_segregrating_Sites; i++)
    {
        Seg_sites.append(collect_Segregrating_site_POS[i].second);
        site_Index[i + 1] = site_Index[i] + collect_Segregrating_site_POS[i].second.size();
    }

    char *full_Char;
    full_Char = (char *)malloc((Seg_sites.size() + 1) * sizeof(char));
    strcpy(full_Char, Seg_sites.c_str());
    collect_Segregrating_site_POS.clear();

    char *cuda_full_Char;
    cudaMallocManaged(&cuda_full_Char, (Seg_sites.size() + 1) * sizeof(char));
    int *cuda_site_Index;
    cudaMallocManaged(&cuda_site_Index, (num_segregrating_Sites + 1) * sizeof(int));
    cudaMemcpy(cuda_full_Char, full_Char, (Seg_sites.size() + 1) * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(cuda_site_Index, site_Index, (num_segregrating_Sites + 1) * sizeof(int), cudaMemcpyHostToDevice);

    // int *REF_Count, *ALT_Count;
    int *cuda_REF_Count, *cuda_ALT_Count;
    cudaMallocManaged(&cuda_REF_Count, num_segregrating_Sites * sizeof(int));
    cudaMallocManaged(&cuda_ALT_Count, num_segregrating_Sites * sizeof(int));
    // REF_Count = (int *)malloc(num_segregrating_Sites * sizeof(int));
    // ALT_Count = (int *)malloc(num_segregrating_Sites * sizeof(int));

    // char *REF, *ALT;
    char *cuda_REF, *cuda_ALT;
    cudaMallocManaged(&cuda_REF, num_segregrating_Sites * sizeof(char));
    cudaMallocManaged(&cuda_ALT, num_segregrating_Sites * sizeof(char));
    // REF = (char *)malloc(num_segregrating_Sites * sizeof(char));
    // ALT = (char *)malloc(num_segregrating_Sites * sizeof(char));

    cuda_process_SNPS<<<tot_Blocks, tot_ThreadsperBlock>>>(cuda_full_Char, cuda_site_Index, num_segregrating_Sites, cuda_REF_Count, cuda_ALT_Count, cuda_REF, cuda_ALT);
    cudaDeviceSynchronize();

    // cudaMemcpy(REF_Count, cuda_REF_Count, num_segregrating_Sites * sizeof(int), cudaMemcpyDeviceToHost);
    // cudaMemcpy(ALT_Count, cuda_ALT_Count, num_segregrating_Sites * sizeof(int), cudaMemcpyDeviceToHost);

    // cudaMemcpy(REF, cuda_REF, num_segregrating_Sites * sizeof(char), cudaMemcpyDeviceToHost);
    // cudaMemcpy(ALT, cuda_ALT, num_segregrating_Sites * sizeof(char), cudaMemcpyDeviceToHost);

    free(full_Char);
    cudaFree(cuda_full_Char);
    cudaFree(cuda_site_Index);

    cout << "             Completed processing " << num_segregrating_Sites << " SNP sites " << endl
         << endl;

    int num_of_Codons = (codon_Stop - codon_Start + 1) / 3;

    int *SEG_positions;
    SEG_positions = (int *)malloc(num_segregrating_Sites * sizeof(int));
    for (size_t i = 0; i < num_segregrating_Sites; i++)
    {
        SEG_positions[i] = collect_Segregrating_site_POS[i].first;
        // cout << SEG_positions[i] << endl;
    }

    cout << "STEP 2 OF 2: Processing " << num_of_Codons << " codons " << endl;

    vector<int> positions;
    vector<char> REF_alleles;
    vector<char> OUTGROUP_alleles;

    functions function = functions();

    cout << "             Loading codon alignemnt index file: " << codon_Index_File_name << endl;
    fstream codon_File;
    codon_File.open(codon_Index_File_name, ios::in);
    string line;
    // skip first line
    getline(codon_File, line);
    while (getline(codon_File, line))
    {
        vector<string> split_Line;
        function.split(split_Line, line, '\t');
        positions.push_back(stoi(split_Line[0]));
        REF_alleles.push_back(split_Line[1].at(0));
        OUTGROUP_alleles.push_back(split_Line[2].at(0));
    }

    int size_of_alignment_File = positions.size();

    int *positions_ARRAY, *cuda_positions_ARRAY;
    char *REF_array, *cuda_REF_array, *Outroup_array, *cuda_Outroup_array;
    positions_ARRAY = (int *)malloc(size_of_alignment_File * sizeof(int));
    REF_array = (char *)malloc(size_of_alignment_File * sizeof(char));
    Outroup_array = (char *)malloc(size_of_alignment_File * sizeof(char));

    for (size_t i = 0; i < size_of_alignment_File; i++)
    {
        positions_ARRAY[i] = positions[i];
        REF_array[i] = REF_alleles[i];
        Outroup_array[i] = OUTGROUP_alleles[i];
    }

    positions.clear();
    REF_alleles.clear();
    OUTGROUP_alleles.clear();

    cout << "             Codon alignemnt index file loaded" << endl;
    cout << "             Priming GPU" << endl;
    cudaMallocManaged(&cuda_positions_ARRAY, size_of_alignment_File * sizeof(int));
    cudaMemcpy(cuda_positions_ARRAY, positions_ARRAY, size_of_alignment_File * sizeof(int), cudaMemcpyHostToDevice);
    cudaMallocManaged(&cuda_REF_array, size_of_alignment_File * sizeof(char));
    cudaMemcpy(cuda_REF_array, REF_array, size_of_alignment_File * sizeof(char), cudaMemcpyHostToDevice);
    cudaMallocManaged(&cuda_Outroup_array, size_of_alignment_File * sizeof(char));
    cudaMemcpy(cuda_Outroup_array, Outroup_array, size_of_alignment_File * sizeof(char), cudaMemcpyHostToDevice);

    char *cuda_index_Gen_code;
    cudaMallocManaged(&cuda_index_Gen_code, size_of_genetic_Code * sizeof(char));
    cudaMemcpy(cuda_index_Gen_code, index_Gen_code, size_of_genetic_Code * sizeof(char), cudaMemcpyHostToDevice);

    int *cuda_Seg_positions; // int *SEG_positions;
    cudaMallocManaged(&cuda_Seg_positions, num_segregrating_Sites * sizeof(int));
    cudaMemcpy(cuda_Seg_positions, SEG_positions, num_segregrating_Sites * sizeof(int), cudaMemcpyHostToDevice);

    int *VALID_or_NOT, *cuda_VALID_or_NOT;
    int *Ds, *Dn, *Ps, *Pn, *cuda_Ds, *cuda_Dn, *cuda_Ps, *cuda_Pn;

    cudaMallocManaged(&cuda_VALID_or_NOT, num_of_Codons * sizeof(int));
    VALID_or_NOT = (int *)malloc(num_of_Codons * sizeof(int));

    cudaMallocManaged(&cuda_Dn, num_of_Codons * sizeof(int));
    Dn = (int *)malloc(num_of_Codons * sizeof(int));

    cudaMallocManaged(&cuda_Ds, num_of_Codons * sizeof(int));
    Ds = (int *)malloc(num_of_Codons * sizeof(int));

    cudaMallocManaged(&cuda_Pn, num_of_Codons * sizeof(int));
    Pn = (int *)malloc(num_of_Codons * sizeof(int));

    cudaMallocManaged(&cuda_Ps, num_of_Codons * sizeof(int));
    Ps = (int *)malloc(num_of_Codons * sizeof(int));

    // Fix seg sites only catch cause Allele freq does not matter.
    // GPU load here
    // cuda_process_Codons(int codon_Number, int *positions, char *REF, char *Outgroup, char *seg_REF, char *seg_ALT, int SEG_size, int *SEG_positions, int *seg_REF_count, int *seg_ALT_count, int codon_Start, int size_of_alignment_File, int genetic_Code_size, char *index_Genetic_code, int *VALID_or_NOT, int *Ds, int *Dn, int *Ps, int *Pn)

    cout << "             Launching GPU" << endl;
    cuda_process_Codons<<<tot_Blocks, tot_ThreadsperBlock>>>(num_of_Codons, cuda_positions_ARRAY, cuda_REF_array, cuda_Outroup_array, cuda_REF, cuda_ALT, num_segregrating_Sites, cuda_Seg_positions, cuda_REF_Count, cuda_ALT_Count, codon_Start, size_of_alignment_File, size_of_genetic_Code, cuda_index_Gen_code, cuda_VALID_or_NOT, cuda_Ds, cuda_Dn, cuda_Ps, cuda_Pn);
    cudaDeviceSynchronize();

    cudaError_t err = cudaGetLastError();

    if (err != cudaSuccess)
    {
        printf("CUDA Error: %s\n", cudaGetErrorString(err));

        // Possibly: exit(-1) if program cannot continue....
    }

    cudaMemcpy(VALID_or_NOT, cuda_VALID_or_NOT, num_of_Codons * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(Dn, cuda_Dn, num_of_Codons * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(Ds, cuda_Ds, num_of_Codons * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(Pn, cuda_Pn, num_of_Codons * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(Ps, cuda_Ps, num_of_Codons * sizeof(int), cudaMemcpyDeviceToHost);

    cudaFree(cuda_Seg_positions);
    cudaFree(cuda_REF_Count);
    cudaFree(cuda_ALT_Count);
    cudaFree(cuda_REF);
    cudaFree(cuda_ALT);

    cudaFree(cuda_positions_ARRAY);
    cudaFree(cuda_REF_array);
    cudaFree(cuda_Outroup_array);
    cudaFree(cuda_index_Gen_code);

    cudaFree(cuda_VALID_or_NOT);
    cudaFree(cuda_Dn);
    cudaFree(cuda_Ds);
    cudaFree(cuda_Pn);
    cudaFree(cuda_Ps);

    cout << "             GPU launched" << endl;

    free(SEG_positions);
    free(positions_ARRAY);
    free(REF_array);
    free(Outroup_array);

    cout << "             Completed processing " << num_of_Codons << " codons" << endl;
    cout << endl;

    // process dn ds pn ps
    tot_Dn = 0;
    tot_Ds = 0;
    tot_Pn = 0;
    tot_Ps = 0;

    for (size_t i = 0; i < num_of_Codons; i++)
    {
        if (VALID_or_NOT[i] == 1)
        {
            tot_Dn = tot_Dn + Dn[i];
            tot_Ds = tot_Ds + Ds[i];
            tot_Pn = tot_Pn + Pn[i];
            tot_Ps = tot_Ps + Ps[i];
        }
    }

    cout << "Total Dn: " << tot_Dn << "\t Total Ds: " << tot_Ds << endl;
    cout << "Total Pn: " << tot_Pn << "\t Total Ps: " << tot_Ps << endl;

    float numerator = (float)tot_Pn / (float)tot_Ps;
    float denominator = (float)tot_Dn / (float)tot_Ds;

    NI = numerator / denominator;
    if (isnan(NI))
    {
        cout << "\nNeutrality Index (NI): "
             << "NA_DIV_0" << endl;
    }
    else
    {
        cout << "\nNeutrality Index (NI): " << NI << endl;
    }

    // after that clear. CHECK ALL VARIABLES ARE CLARED
    free(VALID_or_NOT);
    free(Dn);
    free(Ds);
    free(Ps);
    free(Pn);

    cout << endl;
}

void mk_test::prepration()
{
    set<string> log = get_Log();
    // log_Write("reference_mapping_complete");
    // log_Write("alignment_prep_complete");

    if (log.size() > 0)
    {
        if (mode == "CHROM")
        {
            if (log.find("reference_mapping_complete") != log.end())
            {
                cout << "Reference mapping of codons has already been completed" << endl;
            }
            else if (log.find("alignment_prep_complete") != log.end())
            {
                // vector<pair<int, int>> TEMP_file_index = index_alignment_Folder();
                cout << "Alignment file preperation has already been completed" << endl;
                reference_Prep(index_alignment_Folder());
            }
            else
            {
                reference_Prep(alignment_Prep());
                cout << "Chromosome wide reference mapping of codons complete" << endl;
            }
        }
        else
        {
            if (log.find("reference_mapping_per_GENE_complete") != log.end())
            {
                cout << "Per gene reference mapping of codons has already been completed" << endl;
            }
            else
            {
                reference_Prep();
                cout << "Per gene reference mapping of codons complete" << endl;
            }
        }
    }
    else
    {
        if (mode == "CHROM")
        {
            reference_Prep(alignment_Prep());
            cout << "Chromosome wide reference mapping of codons complete" << endl;
        }
        else
        {
            reference_Prep();
            cout << "Per gene reference mapping of codons complete" << endl;
        }
    }
}

void mk_test::reference_Prep()
{
    // GENE VERSION
    cout << "Preparing reference file: " << this->reference_Path << endl
         << endl;

    // string temp_index_Folder = this->primary_Intermediate_Path + "/" + filesystem::path(this->alignment_Path).stem().string();

    functions function = functions();

    string intermediate_Reference = this->primary_Intermediate_Path + "/" + filesystem::path(this->reference_Path).stem().string();
    if (filesystem::exists(intermediate_Reference) == 0)
    {
        cout << "Creating intermediate reference index folder: " << intermediate_Reference << endl;
        filesystem::create_directory(intermediate_Reference);
    }
    else
    {
        cout << "Intermediate reference index folder exists" << endl;
    }

    fstream reference;
    reference.open(this->reference_Path, ios::in);

    string full_Reference = "";

    if (reference.is_open())
    {
        cout << endl
             << "Loading reference file" << endl;
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
    char *reference_full;
    reference_full = (char *)malloc((full_Reference.size() + 1) * sizeof(char));
    cudaMallocManaged(&cuda_reference, (full_Reference.size() + 1) * sizeof(char));
    this->reference_size = full_Reference.size();
    strcpy(reference_full, full_Reference.c_str());
    cudaMemcpy(cuda_reference, reference_full, (full_Reference.size() + 1) * sizeof(char), cudaMemcpyHostToDevice);

    free(reference_full);

    cout << "Reference file loaded" << endl;

    // INDEX reference genome according to gene list ORF and codons.
    fstream gene_File;
    gene_File.open(this->gene_List, ios::in);

    if (gene_File.is_open())
    {
        cout << endl
             << "Indexing reference and alignments according to gene list: " << this->gene_List << endl
             << endl;

        string gene_Combo;

        while (getline(gene_File, gene_Combo))
        {
            vector<string> split_Data;
            function.split(split_Data, gene_Combo, '\t');

            string alignment_File = split_Data[2];
            string combination = split_Data[1];
            string gene_Name = split_Data[0];
            cout << "Gene name\t: " << gene_Name << endl;
            string temp_index_Folder = this->primary_Intermediate_Path + "/alignments/" + gene_Name;
            cout << "Alignment File\t: " << alignment_File << endl;
            vector<string> coordinates;
            function.split(coordinates, split_Data[1], ':');
            // POS - 1 Since in FASTA 1st position is 1, but in C++ its 0
            int start_Co = stoi(coordinates[1]) - 1;
            int end_Co = stoi(coordinates[2]);
            cout << "Coordinates\t: Chromosome: " << coordinates[0] << " Start: " << start_Co + 1 << " End: " << end_Co << endl;

            // replace(combination.begin(), combination.end(), ':', '_');
            // string File = intermediate_Reference + "/" + split_Data[0] + "_" + combination + ".ca";
            string file_Name = intermediate_Reference + "/" + gene_Name + "_" + coordinates[0] + "_" + to_string(start_Co + 1) + "_" + to_string(end_Co) + ".ca";

            if (filesystem::exists(file_Name) == 0)
            {
                cout << endl;

                int found = 0;
                int ORF_start;
                int ORF_stop;

                if (this->ORF_mode == "NO")
                {
                    // Automatically search for ORF
                    cout << "Initiating Open Reading Frame (ORF) search:" << endl;

                    vector<int> start_ORFs;

                    for (int i = start_Co; i < end_Co; i++)
                    {
                        if ((i + 3) < end_Co)
                        {
                            string check_Current{full_Reference.at(i), full_Reference.at(i + 1), full_Reference.at(i + 2)};
                            for (string check : this->start_Codons_list)
                            {
                                if (check_Current == check)
                                {
                                    start_ORFs.push_back(i);
                                    break;
                                }
                            }
                        }
                        else
                        {
                            break;
                        }
                    }

                    if (start_ORFs.size() > 0)
                    {
                        ORF_search(start_ORFs, found, ORF_start, ORF_stop, end_Co);
                    }
                }
                else
                {
                    // ORF is known and gene config is ORF config
                    found = 1;
                    ORF_start = start_Co;
                    ORF_stop = end_Co - 3;
                }
                if (found == 1)
                {
                    cout << "Start codon location\t: " << setfill('0') << setw(to_string(ORF_stop + 1).length()) << to_string(ORF_start + 1)
                         << "\t Start Codon: " << full_Reference.at(ORF_start) << full_Reference.at(ORF_start + 1) << full_Reference.at(ORF_start + 2) << endl;
                    cout << "Stop codon location\t: " << ORF_stop + 1
                         << "\t Stop Codon: " << full_Reference.at(ORF_stop) << full_Reference.at(ORF_stop + 1) << full_Reference.at(ORF_stop + 2) << endl;
                    cout << endl;

                    // MODIFIED FOR PER GENE

                    // int test_Pop = start_Co - 1 + 700;
                    // int test_Pop = 896666;
                    // cout << "Reference at 896666: " << full_Reference.at(test_Pop - 1) << endl;
                    // cout << "Reference at 896667: " << full_Reference.at(test_Pop - 0) << endl;
                    // cout << "Reference at 896668: " << full_Reference.at(test_Pop + 1) << endl;
                    // cout << "Reference at 896669: " << full_Reference.at(test_Pop + 2) << endl;
                    // cout << "Reference at 896670: " << full_Reference.at(test_Pop + 3) << endl;
                    // for (size_t i = 0; i < TEMP_file_index.size(); i++)
                    // {
                    //     cout << TEMP_file_index[i].first << endl;
                    // }

                    cout << "Generating target ORF alignment:" << endl;
                    codon_Alignment_print(alignment_Prep(alignment_File, full_Reference, start_Co, temp_index_Folder), ORF_start + 1, ORF_stop + 3, temp_index_Folder, intermediate_Reference, file_Name);
                    // PURGE THE TEMP FOLDER
                }
                else
                {
                    cout << "ORF NOT FOUND" << endl;
                }
            }
            cout << endl;
        }
        gene_File.close();
    }
    cudaFree(cuda_stop_Codons);
    cudaFree(cuda_reference);

    cout << "Purging temporary alignment index: " << this->primary_Intermediate_Path + "/alignments" << endl
         << endl;
    filesystem::remove_all(this->primary_Intermediate_Path + "/alignments");

    // REMOVE AFTER TEST
    log_Write("reference_mapping_per_GENE_complete");
}

vector<string> mk_test::alignment_Prep(string Gene_alignment_Path, string &full_Reference, int start_Co, string temp_index_Folder)
{
    functions function = functions();

    vector<string> file_List;
    vector<pair<int, int>> TEMP_file_index;

    cout << "STEP 1 OF 3: Processing alignment file" << endl;
    fstream alignment_File;
    alignment_File.open(Gene_alignment_Path, ios::in);

    cout << "STEP 2 OF 3: Creating temporary index" << endl;
    // string temp_index_Folder = this->intermediate_Path + "/" + gene_Name;
    if (filesystem::exists(temp_index_Folder) != 0)
    {
        filesystem::remove_all(temp_index_Folder);
    }
    filesystem::create_directory(temp_index_Folder);

    if (alignment_File.is_open())
    {
        string line;

        while (getline(alignment_File, line))
        {
            if (line.length() > 0)
            {
                if (line.at(0) == '>')
                {
                    break;
                }
            }
        }

        while (getline(alignment_File, line))
        {
            if (line.find("Query") != string::npos)
            {
                string query = line;
                getline(alignment_File, line);
                string subject;
                getline(alignment_File, subject);

                // cout << query << endl
                //      << subject << endl;
                // cout << endl;

                vector<string> split_Query;
                function.split_space(split_Query, query, "  ");
                transform(split_Query[2].begin(), split_Query[2].end(), split_Query[2].begin(), ::toupper);

                vector<string> split_Subject;
                function.split_space(split_Subject, subject, "  ");
                transform(split_Subject[2].begin(), split_Subject[2].end(), split_Subject[2].begin(), ::toupper);

                // cout << split_Query[1] << "\t" << split_Query[2] << endl
                //      << split_Subject[1] << "\t" << split_Subject[2] << endl;
                // cout << endl;

                // int test_Pop = start_Co - 1 + 700;
                int ref_Pos = start_Co - 1 + stoi(split_Query[1]);
                int iterate_Value = ref_Pos + 1;

                string temp_File_name = temp_index_Folder + "/" + to_string(iterate_Value) + ".temp";
                fstream write_Temp;
                write_Temp.open(temp_File_name, ios::out);

                for (int i = 0; i < split_Query[2].length(); i++)
                {
                    if (split_Query[2].at(i) != '-')
                    {
                        if (split_Query[2].at(i) == 'A' || split_Query[2].at(i) == 'T' || split_Query[2].at(i) == 'G' || split_Query[2].at(i) == 'C')
                        {
                            if (split_Subject[2].at(i) == 'A' || split_Subject[2].at(i) == 'T' || split_Subject[2].at(i) == 'G' || split_Subject[2].at(i) == 'C')
                            {
                                // convert char to string
                                string ref(1, split_Query[2].at(i));
                                string query(1, split_Subject[2].at(i));
                                write_Temp << to_string(iterate_Value) << "\t" << ref << "\t" << query << "\n";
                            }
                        }
                        iterate_Value++;
                    }
                }
                write_Temp.close();
                string new_Name = temp_index_Folder + "/" + to_string(ref_Pos + 1) + "_" + to_string(iterate_Value - 1) + ".temp";
                rename(temp_File_name.c_str(), new_Name.c_str());
                file_List.push_back(to_string(ref_Pos + 1) + "_" + to_string(iterate_Value - 1) + ".temp");
                // TEMP_file_index.push_back(make_pair(ref_Pos + 1, iterate_Value - 1));
            }
        }
        alignment_File.close();
    }

    // sort(TEMP_file_index.begin(), TEMP_file_index.end());
    cout << "STEP 3 OF 3: Completed alignment file processing" << endl;
    cout << endl;

    return file_List;
}

vector<pair<int, int>> mk_test::index_alignment_Folder()
{
    functions function = functions();
    string temp_index_Folder = this->primary_Intermediate_Path + "/" + filesystem::path(this->alignment_Path).stem().string();

    vector<pair<int, int>> TEMP_file_index;

    for (const auto &entry : filesystem::directory_iterator(temp_index_Folder))
    {
        string file_Name = entry.path().string();
        int index_of_slash = file_Name.find_last_of('/');
        int index_of_dot = file_Name.find_last_of('.');
        file_Name = file_Name.substr(index_of_slash + 1, index_of_dot - index_of_slash - 1);
        // cout << file_Name << endl;
        vector<string> split_Data;
        function.split(split_Data, file_Name, '_');
        TEMP_file_index.push_back(make_pair(stoi(split_Data[0]), stoi(split_Data[1])));
    }

    sort(TEMP_file_index.begin(), TEMP_file_index.end());
    return TEMP_file_index;
}

set<string> mk_test::get_Log()
{
    functions function = functions();
    set<string> log;
    string prep_Log_path = this->intermediate_Path + "/" + filesystem::path(this->gene_List).stem().string() + ".log_mk_prep";

    if (filesystem::exists(prep_Log_path) != 0)
    {
        fstream prep_Log;
        prep_Log.open(prep_Log_path, ios::in);
        if (prep_Log.is_open())
        {
            string line;
            while (getline(prep_Log, line))
            {
                log.insert(line);
            }
            prep_Log.close();
        }
    }
    else
    {
        function.createFile(prep_Log_path);
    }

    return log;
}

void mk_test::log_Write(string line)
{
    string prep_Log = this->intermediate_Path + "/" + filesystem::path(this->gene_List).stem().string() + ".log_mk_prep";
    fstream log;
    log.open(prep_Log, ios::app);
    log << line << "\n";
    log.close();
}

vector<pair<int, int>> mk_test::alignment_Prep()
{
    cout << "Preparing alignment file: " << this->alignment_Path << endl
         << endl;

    functions function = functions();

    fstream alignment_File;
    alignment_File.open(this->alignment_Path, ios::in);

    vector<pair<int, string>> sort_Alinged_Full;

    if (alignment_File.is_open())
    {
        cout << "STEP 1 OF 5: Reading and processing alignment file" << endl;
        string line;

        int ref_Check = 0;
        int ref_Line_start;
        // start ref_seq query_seq
        string ref_Line;

        while (getline(alignment_File, line))
        {
            if (line.length() > 0)
            {
                if (line.at(0) == 's')
                {
                    vector<string> split_Line;
                    function.split_space(split_Line, line, " ");
                    transform(split_Line[6].begin(), split_Line[6].end(), split_Line[6].begin(), ::toupper);

                    if (ref_Check == 0)
                    {
                        // reference line
                        // cout << split_Line[1] << "\t";
                        ref_Line_start = stoi(split_Line[2]);
                        // int ref_Line_stop = ref_Line_start + stoi(split_Line[3]);
                        ref_Line = split_Line[2] + "\t" + split_Line[6];
                        ref_Check++;
                    }
                    else
                    {
                        // query line
                        // cout << split_Line[1] << endl;
                        ref_Line = ref_Line + "\t" + split_Line[6];
                        sort_Alinged_Full.push_back(make_pair(ref_Line_start, ref_Line));
                        ref_Line = "";
                        ref_Line_start = 0;
                        ref_Check = 0;
                    }
                }
            }
        }
        alignment_File.close();
    }

    cout << "STEP 2 OF 5: Sorting alignment file" << endl;
    sort(sort_Alinged_Full.begin(), sort_Alinged_Full.end());

    vector<pair<int, string>> sort_ref_query;
    set<char> bases;

    string temp_index_Folder = this->primary_Intermediate_Path + "/" + filesystem::path(this->alignment_Path).stem().string();
    cout << "STEP 3 OF 5: Creating temporary alignment index folder: " << temp_index_Folder << endl;

    if (filesystem::exists(temp_index_Folder) != 0)
    {
        cout << "             Folder already exists: PURGING Folder" << endl;
        filesystem::remove_all(temp_index_Folder);
    }

    filesystem::create_directory(temp_index_Folder);
    cout << "STEP 4 OF 5: Temporary alinment folder created" << endl;

    cout << "STEP 5 OF 5: Writing: Base wise mapping of reference to alignment (please wait)" << endl;
    vector<pair<int, int>> TEMP_file_index;

    for (size_t i = 0; i < sort_Alinged_Full.size(); i++)
    {
        // start ref_seq query_seq
        vector<string> split_Line;
        function.split(split_Line, sort_Alinged_Full[i].second, '\t');

        // maf start pos is a zero based number. Therefore it must be incremented by 1
        int ref_start_Pos = stoi(split_Line[0]) + 1;
        string temp_Index = temp_index_Folder + "/" + to_string(ref_start_Pos) + ".temp";
        fstream temp_Write;
        temp_Write.open(temp_Index, ios::out);

        int ref_POS = 0;
        int last_file_Pos = 0;
        for (size_t iterate = 0; iterate < split_Line[1].length(); iterate++)
        {
            // write to hard disk
            if (split_Line[1].at(iterate) != '-')
            {
                if (split_Line[1].at(iterate) == 'A' || split_Line[1].at(iterate) == 'T' || split_Line[1].at(iterate) == 'G' || split_Line[1].at(iterate) == 'C')
                {
                    if (split_Line[2].at(iterate) == 'A' || split_Line[2].at(iterate) == 'T' || split_Line[2].at(iterate) == 'G' || split_Line[2].at(iterate) == 'C')
                    {
                        string ref(1, split_Line[1].at(iterate));
                        string query(1, split_Line[2].at(iterate));
                        string write = to_string(ref_start_Pos + ref_POS) + "\t" + ref + "\t" + query;
                        temp_Write << write << "\n";
                        last_file_Pos = ref_start_Pos + ref_POS;
                    }
                }
                ref_POS++;
            }
        }
        temp_Write.close();
        string new_name = temp_index_Folder + "/" + to_string(ref_start_Pos) + "_" + to_string(last_file_Pos) + ".temp";
        rename(temp_Index.c_str(), new_name.c_str());
        TEMP_file_index.push_back(make_pair(ref_start_Pos, last_file_Pos));
    }

    cout << endl
         << "Completed alignment file preperation" << endl
         << endl;

    log_Write("alignment_prep_complete");
    return TEMP_file_index;
}

void mk_test::reference_Prep(vector<pair<int, int>> TEMP_file_index)
{
    cout << "Preparing reference file: " << this->reference_Path << endl
         << endl;

    string temp_index_Folder = this->primary_Intermediate_Path + "/" + filesystem::path(this->alignment_Path).stem().string();

    functions function = functions();

    string intermediate_Reference = this->primary_Intermediate_Path + "/" + filesystem::path(this->reference_Path).stem().string();
    if (filesystem::exists(intermediate_Reference) == 0)
    {
        cout << "Creating intermediate reference index folder: " << intermediate_Reference << endl;
        filesystem::create_directory(intermediate_Reference);
    }
    else
    {
        cout << "Intermediate reference index folder exists" << endl;
    }

    fstream reference;
    reference.open(this->reference_Path, ios::in);

    string full_Reference = "";

    if (reference.is_open())
    {
        cout << endl
             << "Loading reference file" << endl;
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
    char *reference_full;
    reference_full = (char *)malloc((full_Reference.size() + 1) * sizeof(char));
    cudaMallocManaged(&cuda_reference, (full_Reference.size() + 1) * sizeof(char));
    this->reference_size = full_Reference.size();
    strcpy(reference_full, full_Reference.c_str());
    cudaMemcpy(cuda_reference, reference_full, (full_Reference.size() + 1) * sizeof(char), cudaMemcpyHostToDevice);

    free(reference_full);

    cout << "Reference file loaded" << endl;
    // cout << full_Reference.at(565322 - 1) << endl;

    // INDEX reference genome according to gene list ORF and codons.
    fstream gene_File;
    gene_File.open(this->gene_List, ios::in);

    if (gene_File.is_open())
    {
        cout << endl
             << "Indexing reference according to gene list: " << this->gene_List << endl
             << endl;

        string gene_Combo;

        while (getline(gene_File, gene_Combo))
        {
            vector<string> split_Data;
            function.split(split_Data, gene_Combo, '\t');

            string combination = split_Data[1];
            string gene_Name = split_Data[0];
            cout << "Gene name\t: " << gene_Name << endl;
            vector<string> coordinates;
            function.split(coordinates, split_Data[1], ':');
            // POS - 1 Since in FASTA 1st position is 1, but in C++ its 0
            int start_Co = stoi(coordinates[1]) - 1;
            int end_Co = stoi(coordinates[2]);
            cout << "Coordinates\t: Chromosome: " << coordinates[0] << " Start: " << start_Co + 1 << " End: " << end_Co << endl;

            // replace(combination.begin(), combination.end(), ':', '_');
            // string File = intermediate_Reference + "/" + split_Data[0] + "_" + combination + ".ca";
            string file_Name = intermediate_Reference + "/" + gene_Name + "_" + coordinates[0] + "_" + to_string(start_Co + 1) + "_" + to_string(end_Co) + ".ca";

            if (filesystem::exists(file_Name) == 0)
            {
                cout << endl;
                int found = 0;
                int ORF_start;
                int ORF_stop;

                if (this->ORF_mode == "NO")
                {

                    cout << "Initiating Open Reading Frame (ORF) search:" << endl;

                    vector<int> start_ORFs;

                    for (int i = start_Co; i < end_Co; i++)
                    {
                        if ((i + 3) < end_Co)
                        {
                            string check_Current{full_Reference.at(i), full_Reference.at(i + 1), full_Reference.at(i + 2)};
                            for (string check : this->start_Codons_list)
                            {
                                if (check_Current == check)
                                {
                                    // cout << check << endl;
                                    // ORF_start = i;
                                    // start_ORF_end = i + 3;
                                    // start_Found = 1;
                                    start_ORFs.push_back(i);
                                    break;
                                }
                            }
                        }
                        else
                        {
                            break;
                        }
                        // if (start_Found == 1)
                        // {
                        //     break;
                        // }
                    }

                    // if (start_Found == 1)
                    // {
                    //     for (int i = start_ORF_end; i < end_Co; i = i + 3)
                    //     {
                    //         if ((i + 3) < end_Co)
                    //         {
                    //             string check_Current{full_Reference.at(i), full_Reference.at(i + 1), full_Reference.at(i + 2)};
                    //             for (string check : this->stop_Codons_list)
                    //             {
                    //                 if (check_Current == check)
                    //                 {
                    //                     ORF_stop = i;
                    //                     stop_Found = 1;
                    //                     break;
                    //                 }
                    //             }
                    //         }
                    //         else
                    //         {
                    //             break;
                    //         }
                    //     }
                    // }

                    // for (int i = start_Co; i < end_Co; i++)
                    // {
                    //     if ((i + 3) < end_Co)
                    //     {
                    //         string check_Current{full_Reference.at(i), full_Reference.at(i + 1), full_Reference.at(i + 2)};
                    //         if (start_Found == 0)
                    //         {
                    //             for (string check : this->start_Codons_list)
                    //             {
                    //                 if (check_Current == check)
                    //                 {
                    //                     //cout << check << endl;
                    //                     ORF_start = i;
                    //                     start_Found = 1;
                    //                     break;
                    //                 }
                    //             }
                    //         }
                    //         else
                    //         {
                    //             //look for stop codons only if start codon has been found
                    //             for (string check : this->stop_Codons_list)
                    //             {
                    //                 if (check_Current == check)
                    //                 {
                    //                     ORF_stop = i;
                    //                     stop_Found = 1;
                    //                     break;
                    //                 }
                    //             }
                    //         }
                    //     }
                    //     else
                    //     {
                    //         break;
                    //     }
                    // }

                    // int found = 0;
                    //  int start_Found = 0;
                    //   int start_ORF_end;
                    //  int stop_Found = 0;

                    // int ORF_start;
                    // int ORF_stop;

                    if (start_ORFs.size() > 0)
                    {
                        ORF_search(start_ORFs, found, ORF_start, ORF_stop, end_Co);
                    }
                    // cout << endl;
                }
                else
                {
                    found = 1;
                    ORF_start = start_Co;
                    ORF_stop = end_Co - 3;
                }
                if (found == 1)
                {
                    cout << "Start codon location\t: " << setfill('0') << setw(to_string(ORF_stop + 1).length()) << to_string(ORF_start + 1)
                         << "\t Start Codon: " << full_Reference.at(ORF_start) << full_Reference.at(ORF_start + 1) << full_Reference.at(ORF_start + 2) << endl;
                    cout << "Stop codon location\t: " << ORF_stop + 1
                         << "\t Stop Codon: " << full_Reference.at(ORF_stop) << full_Reference.at(ORF_stop + 1) << full_Reference.at(ORF_stop + 2) << endl;
                    cout << endl;
                    cout << "Fetching target ORF alignment:" << endl;
                    // vector<string> file_List = compound_Interpolation_folder(TEMP_file_index, ORF_start + 1, ORF_stop + 3);
                    cout << "Collecting alignment files" << endl;
                    codon_Alignment_print(compound_Interpolation_folder(TEMP_file_index, ORF_start + 1, ORF_stop + 3), ORF_start + 1, ORF_stop + 3, temp_index_Folder, intermediate_Reference, file_Name);
                    // align regions from alignment file
                }
                else
                {
                    cout << "ORF NOT FOUND" << endl;
                }
                // codon split
            }
            cout << endl;
        }
        gene_File.close();
    }
    cout << "Purging temporary alignment index: " << temp_index_Folder << endl
         << endl;
    filesystem::remove_all(temp_index_Folder);
    cudaFree(cuda_stop_Codons);
    cudaFree(cuda_reference);

    log_Write("reference_mapping_complete");
}

__global__ void cuda_ORF_search(int ORF_nums, int *start_ORFs, char *reference_Full, int reference_Length, char *stop_Codons, int stop_Codons_Length, int gene_End, int *VALID_or_NOT, int *end_ORFs)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < ORF_nums)
    {
        int end_Search = gene_End;

        if (reference_Length < gene_End)
        {
            end_Search = reference_Length;
        }

        int start_Search = start_ORFs[tid] + 3;

        int found = 0;
        int stop_Pos;

        for (int i = start_Search; i < end_Search; i = i + 3)
        {
            if ((i + 3) < end_Search)
            {
                char CODON_pos_1 = reference_Full[i];
                char CODON_pos_2 = reference_Full[i + 1];
                char CODON_pos_3 = reference_Full[i + 2];

                // search stops
                for (int stops = 0; stops < stop_Codons_Length; stops = stops + 3)
                {
                    char STOP_pos_1 = stop_Codons[stops];
                    char STOP_pos_2 = stop_Codons[stops + 1];
                    char STOP_pos_3 = stop_Codons[stops + 2];

                    if ((CODON_pos_1 == STOP_pos_1) && (CODON_pos_2 == STOP_pos_2) && (CODON_pos_3 == STOP_pos_3))
                    {
                        found = 1;
                        stop_Pos = i;
                        break;
                    }
                }
            }
            else
            {
                break;
            }

            if (found == 1)
            {
                break;
            }
        }

        if (found == 1)
        {
            VALID_or_NOT[tid] = 1;
            end_ORFs[tid] = stop_Pos;
        }
        else
        {
            VALID_or_NOT[tid] = 0;
        }

        tid += blockDim.x * gridDim.x;
    }
}

void mk_test::ORF_search(vector<int> start_ORFs, int &found, int &ORF_start, int &ORF_stop, int gene_End)
{
    int potential_ORFs = start_ORFs.size();

    int *start_ORF_array, *cuda_start_ORF_array;
    start_ORF_array = (int *)malloc(potential_ORFs * sizeof(int));
    cudaMallocManaged(&cuda_start_ORF_array, potential_ORFs * sizeof(int));

    int *VALID_or_NOT, *cuda_VALID_or_NOT, *end_ORFs, *cuda_end_ORFs;
    cudaMallocManaged(&cuda_VALID_or_NOT, potential_ORFs * sizeof(int));
    cudaMallocManaged(&cuda_end_ORFs, potential_ORFs * sizeof(int));
    VALID_or_NOT = (int *)malloc(potential_ORFs * sizeof(int));
    end_ORFs = (int *)malloc(potential_ORFs * sizeof(int));

    for (size_t i = 0; i < potential_ORFs; i++)
    {
        start_ORF_array[i] = start_ORFs[i];
    }

    cudaMemcpy(cuda_start_ORF_array, start_ORF_array, potential_ORFs * sizeof(int), cudaMemcpyHostToDevice);

    // cuda_ORF_search(int ORF_nums, int *start_ORFs, char *reference_Full, int reference_Length, char *stop_Codons, int stop_Codons_Length, int gene_End, int *VALID_or_NOT, int *end_ORFs)
    cuda_ORF_search<<<tot_Blocks, tot_ThreadsperBlock>>>(potential_ORFs, cuda_start_ORF_array, cuda_reference, reference_size, cuda_stop_Codons, stop_Codon_size, gene_End, cuda_VALID_or_NOT, cuda_end_ORFs);
    cudaDeviceSynchronize();

    cudaMemcpy(VALID_or_NOT, cuda_VALID_or_NOT, potential_ORFs * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(end_ORFs, cuda_end_ORFs, potential_ORFs * sizeof(int), cudaMemcpyDeviceToHost);

    int max = -1;
    int ORF_pos = -1;
    for (size_t i = 1; i < potential_ORFs; i++)
    {
        if (VALID_or_NOT[i] == 1)
        {
            // cout << "Potential start codon location\t: " << setfill('0') << setw(to_string(start_ORFs[i] + 1).length()) << start_ORFs[i] + 1
            //<< "\t Start Codon\t: " << full_Reference.at(start_ORFs[i]) << full_Reference.at(start_ORFs[i] + 1) << full_Reference.at(start_ORFs[i] + 2) << endl;
            // cout << "Potential stop codon location\t: " << end_ORFs[i] + 1
            //<< "\t Stop Codon\t: " << full_Reference.at(end_ORFs[i]) << full_Reference.at(end_ORFs[i] + 1) << full_Reference.at(end_ORFs[i] + 2) << endl;
            // print potential start and stop codons with coordinates
            int diff = end_ORFs[i] - start_ORFs[i];
            if (diff > max)
            {
                max = diff;
                ORF_pos = i;
            }
            // cout << endl;
        }
    }
    cout << potential_ORFs << " potential ORFs considered";
    if (max != -1)
    {
        found = 1;
        ORF_start = start_ORFs[ORF_pos];
        ORF_stop = end_ORFs[ORF_pos];
    }

    free(VALID_or_NOT);
    free(end_ORFs);
    free(start_ORF_array);

    cudaFree(cuda_VALID_or_NOT);
    cudaFree(cuda_end_ORFs);
    cudaFree(cuda_start_ORF_array);
    cout << endl;
}

void mk_test::codon_Alignment_print(vector<string> file_List, int start_Codon, int end_Codon, string &temp_index_Folder, string &intermediate_Reference, string &file_Name)
{
    cout << "Collecting codon alignments" << endl;
    functions function = functions();
    // int = location, string =LOC\t REF \t Query;
    vector<pair<int, string>> write_To_file;
    vector<int> locations;

    for (string file : file_List)
    {
        // cout << file << endl;
        fstream align_Index;
        align_Index.open(temp_index_Folder + "/" + file, ios::in);
        // cout << file << endl;
        if (align_Index.is_open())
        {
            string alignment;
            while (getline(align_Index, alignment))
            {
                vector<string> line_Split;
                function.split(line_Split, alignment, '\t');
                int location = stoi(line_Split[0]);
                // cout << location << endl;
                if (location >= start_Codon && location <= end_Codon)
                {
                    // vector<int>::iterator itr = std::find(locations.begin(), locations.end(), location);
                    if (binary_search(locations.begin(), locations.end(), location))
                    {
                    }
                    else
                    {
                        write_To_file.push_back(make_pair(location, alignment));
                        locations.push_back(location);
                        sort(locations.begin(), locations.end());
                    }
                }
                else if (location > end_Codon)
                {
                    break;
                }
            }
            align_Index.close();
        }
    }

    locations.clear();
    sort(write_To_file.begin(), write_To_file.end());

    cout << "Writing codon alignments" << endl;
    fstream codon_Alignment;
    // string file_Name = intermediate_Reference + "/" + gene_Name + "_" + chromosome + "_" + start + "_" + stop + ".ca";
    codon_Alignment.open(file_Name, ios::out);
    codon_Alignment << "##Codon_Range:\t" << start_Codon << "\t" << end_Codon << "\n";
    for (size_t i = 0; i < write_To_file.size(); i++)
    {
        codon_Alignment << write_To_file[i].second << "\n";
    }
    codon_Alignment.close();
}

vector<string> mk_test::compound_Interpolation_folder(vector<pair<int, int>> folder_Index, int start_Co, int end_Co)
{
    // functions function = functions();
    vector<string> file_List;
    int start = 0;
    int end = folder_Index.size() - 1;

    int low_Value = folder_Index[start].first;
    int high_Value = folder_Index[end].second;

    // while (start <= end && start_Co >= low_Value && start_Co <= high_Value)
    while (start <= end)
    {
        int pos = start + ((double)(end - start) / ((high_Value - low_Value)) * (start_Co - low_Value));
        int low_Value_atpos = folder_Index[pos].first;
        int high_Value_atpos = folder_Index[pos].second;

        // cout << low_Value_atpos << "\t" << high_Value_atpos << endl;
        // cout << start_Co << endl;
        // cout << end_Co << endl;

        if ((start_Co >= low_Value_atpos && start_Co <= high_Value_atpos) || (start_Co <= low_Value_atpos && end_Co <= high_Value_atpos) || (start_Co <= low_Value_atpos && end_Co >= high_Value_atpos))
        {
            // cout << "Caught: " << low_Value_atpos << endl;
            promise<vector<int>> backward;
            promise<vector<int>> forward;

            future<vector<int>> fut_Backward = backward.get_future();
            future<vector<int>> fut_Forward = forward.get_future();

            thread backward_thread{&mk_test::backward_Search, this, ref(backward), pos, folder_Index, start_Co, end_Co};
            thread forward_thread{&mk_test::forward_Search, this, ref(forward), pos, folder_Index, start_Co, end_Co};

            vector<int> backward_get = fut_Backward.get();
            vector<int> forward_get = fut_Forward.get();

            backward_thread.join();
            forward_thread.join();

            // vector<string> backward_get;
            // vector<string> forward_get;
            // function.split(backward_get, backward_get_s, ",");
            // function.split(forward_get, forward_get_s, ",");

            for (int positions = backward_get.size() - 1; positions >= 0; positions--)
            {
                // int positions = stoi(position);
                // cout << folder_Index[positions].first << "_" << folder_Index[positions].second << ".temp" << endl;
                file_List.push_back(to_string(folder_Index[positions].first) + "_" + to_string(folder_Index[positions].second) + ".temp");
            }

            // cout << folder_Index[pos].first << "_" << folder_Index[pos].second << ".temp" << endl;
            file_List.push_back(to_string(folder_Index[pos].first) + "_" + to_string(folder_Index[pos].second) + ".temp");

            for (auto positions : forward_get)
            {
                // int positions = stoi(position);
                // cout << folder_Index[positions].first << "_" << folder_Index[positions].second << ".temp" << endl;
                file_List.push_back(to_string(folder_Index[positions].first) + "_" + to_string(folder_Index[positions].second) + ".temp");
            }

            break;
        }
        else if (start_Co > low_Value_atpos)
        {
            int new_pos = pos;

            do
            {
                new_pos = new_pos + 1;
            } while (new_pos <= start);

            start = new_pos;
        }
        else
        {
            int new_pos = pos;

            do
            {
                new_pos = new_pos - 1;
            } while (new_pos >= end);

            end = new_pos;
        }
        low_Value = folder_Index[start].first;
        high_Value = folder_Index[end].second;
    }
    return file_List;
}

void mk_test::backward_Search(promise<vector<int>> &backward_Found, int pos, vector<pair<int, int>> folder_Index, int start_Co, int end_Co)
{
    vector<int> backward_get;
    // string backward_Get = "";
    pos = pos - 1;

    while (pos >= 0)
    {
        int low_Value_atpos = folder_Index[pos].first;
        int high_Value_atpos = folder_Index[pos].second;

        if (start_Co > high_Value_atpos)
        {
            break;
        }

        if ((start_Co >= low_Value_atpos && start_Co <= high_Value_atpos) || (start_Co <= low_Value_atpos && end_Co <= high_Value_atpos) || (start_Co <= low_Value_atpos && end_Co >= high_Value_atpos))
        {
            backward_get.push_back(pos);
            // backward_Get = backward_Get + to_string(pos) + ",";
        }

        pos = pos - 1;
    }

    // backward_Found.set_value(backward_Get.substr(0, backward_Get.length() - 1));
    backward_Found.set_value(backward_get);
    // backward_Found.(backward_get);
}

void mk_test::forward_Search(promise<vector<int>> &forward_Found, int pos, vector<pair<int, int>> folder_Index, int start_Co, int end_Co)
{
    vector<int> forward_get;
    // string forward_Get = "";
    pos = pos + 1;

    while (pos < folder_Index.size())
    {
        int low_Value_atpos = folder_Index[pos].first;
        int high_Value_atpos = folder_Index[pos].second;

        if (end_Co < low_Value_atpos)
        {
            break;
        }

        if ((start_Co >= low_Value_atpos && start_Co <= high_Value_atpos) || (start_Co <= low_Value_atpos && end_Co <= high_Value_atpos) || (start_Co <= low_Value_atpos && end_Co >= high_Value_atpos))
        {
            forward_get.push_back(pos);
            // forward_Get = forward_Get + to_string(pos) + ",";
        }

        pos = pos + 1;
    }

    // forward_Found.set_value(forward_Get.substr(0, forward_Get.length() - 1));
    forward_Found.set_value(forward_get);
}

void mk_test::print_Code(vector<string> &Code_split)
{
    functions function = functions();

    cout << "User specified Genetic code:" << endl;
    function.split(Code_split, this->genetic_Code, ';');

    int long_Length = 0;
    for (string longest : Code_split)
    {
        if (longest.length() > long_Length)
        {
            long_Length = longest.length();
        }
    }

    int field = 0;
    int column = 0;
    while (field < Code_split.size())
    {
        cout << setfill(' ') << setw(long_Length) << left << Code_split[field];
        cout << "\t";

        column++;
        if (column == 4)
        {
            cout << "\n";
            column = 0;
        }

        field++;
    }

    cout << "\n\n";
}