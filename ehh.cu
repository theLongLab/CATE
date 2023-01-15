#include "functions.cuh"
#include "ehh.cuh"

ehh::ehh(string range_Mode, string file_Mode_path, string fixed_Mode_value, string input_Folder, string output_Path, int cuda_ID, string intermediate_Path, int ploidy, int default_SNP_count, int EHH_CPU_cores, int default_SNP_BP_count)
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
    else if (range_Mode == "SNP")
    {
        this->default_SNP_count = default_SNP_count;
        // this->EHH_cutoff = EHH_cutoff;
        this->CPU_cores = EHH_CPU_cores;

        cout << "CPU cores: " << this->CPU_cores << endl
             << endl;
    }
    else if (range_Mode == "BP")
    {
        this->default_SNP_BP_count = default_SNP_BP_count;
        this->CPU_cores = EHH_CPU_cores;

        cout << "CPU cores: " << this->CPU_cores << endl
             << endl;
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
    else if (this->range_Mode == "SNP")
    {
        cout << "Default number of SNPs collected: " << this->default_SNP_count << endl;
    }
    else if (this->range_Mode == "BP")
    {
        cout << "Default displacement from core SNP (bp): " << this->default_SNP_BP_count << endl;
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
                             filesystem::path(file_Mode_path).stem().string();
        if (range_Mode == "FILE" || range_Mode == "FIXED")
        {
            // output_File = ouput_Path + "/" +
            //               country.substr(country.find_last_of("/") + 1, country.length()) + "_" +
            //               filesystem::path(file_Mode_path).stem().string() +
            //               ".ehh";
            output_File = output_File + ".ehh";

            cout << endl;
            cout << "Writing to file\t: " << output_File << endl;
            cout << endl;
        }
        else
        {
            cout << endl;
            cout << "Writing to folder\t: " << output_File << endl;
            cout << endl;
        }
        string intermediate_File = intermediate_Path + "/" +
                                   country.substr(country.find_last_of("/") + 1, country.length()) + "_" +
                                   filesystem::path(file_Mode_path).stem().string() +
                                   ".log_ehh";

        if (gene_File.is_open())
        {
            string gene_Combo;

            if (filesystem::exists(intermediate_File) == 0)
            {
                if (range_Mode == "FILE" || range_Mode == "FIXED")
                {
                    function.createFile(output_File, "Gene_name\tCore_coordinates\tExtended_coordinates\tCore_Haplotype_Number\tCt\tTotal_Et\tEHH");
                }
                else
                {
                    filesystem::create_directory(output_File);
                }
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
            if (range_Mode == "FILE" || range_Mode == "FIXED")
            {
                output.open(output_File, ios::app);
            }
            intermediate.open(intermediate_File, ios::app);

            while (getline(gene_File, gene_Combo))
            {
                vector<string> split_Data;
                function.split(split_Data, gene_Combo, '\t');
                string gene_Name = split_Data[0];

                vector<string> Core_coordinates;
                function.split(Core_coordinates, split_Data[1], ':');

                if (range_Mode == "SNP" || range_Mode == "BP")
                {
                    cout << "SNP ID\t: " << gene_Name << endl;
                    int position_SNP = stoi(Core_coordinates[1]);

                    cout << "Coordinates:" << endl;
                    cout << "Chromosome: " << Core_coordinates[0] << "\tPosition: " << position_SNP << endl;

                    // create a new search for single position.

                    int segment_Position = function.single_segment_retrieval(position_SNP, folder_Index);

                    if (segment_Position != -1)
                    {

                        string SNP_file_Name = folder_Index[segment_Position].second;
                        // cout << file_Name << endl;

                        fstream SNP_file;
                        SNP_file.open(SNP_file_Name, ios::in);

                        int SNP_Index_in_file = -1;

                        if (range_Mode == "SNP")
                        {
                            vector<string> collect_Segregrating_sites;

                            int bottom_count = 0;
                            int top_count = 0;

                            if (SNP_file.is_open())
                            {
                                string line;
                                getline(SNP_file, line); // skip first header line
                                while (getline(SNP_file, line))
                                {
                                    collect_Segregrating_sites.push_back(line);

                                    if (SNP_Index_in_file == -1)
                                    {
                                        vector<string> positions;
                                        function.split_getPos_ONLY(positions, line, '\t');
                                        int pos = stoi(positions[1]);

                                        if (pos == position_SNP)
                                        {
                                            cout << "SNP found in repository: " << SNP_file_Name << endl;
                                            SNP_Index_in_file = collect_Segregrating_sites.size() - 1;
                                            top_count = SNP_Index_in_file;
                                            // cout << SNP_Position_in_file << endl;
                                        }
                                    }
                                    else
                                    {
                                        bottom_count = bottom_count + 1;
                                        if (bottom_count == this->default_SNP_count)
                                        {
                                            break;
                                        }
                                    }
                                }

                                SNP_file.close();
                            }

                            // cout << top_count << endl
                            //      << bottom_count << endl;

                            if (SNP_Index_in_file != -1)
                            {
                                process_SNP_EHH(position_SNP, folder_Index, segment_Position, SNP_Index_in_file, collect_Segregrating_sites, N_float, top_count, bottom_count);

                                // Write to file
                                string results_File = output_File + "/" +
                                                      gene_Name + "_" + Core_coordinates[0] + "_" + Core_coordinates[1] + "_" + to_string(this->default_SNP_count) + "_SNP.ehh";
                                function.createFile(results_File, "Displacement\tSNP_position\tEHH_0\tEHH_1");

                                cout << "Writing results to file: " << results_File << endl
                                     << endl;

                                output.open(results_File, ios::app);

                                for (int line = 0; line < positions_Collect.size(); line++)
                                {
                                    output << positions_Collect[line] << "\t";
                                    if (line < EHH_0_up.size())
                                    {
                                        output << EHH_0_up[EHH_0_up.size() - 1 - line] << "\t" << EHH_1_up[EHH_1_up.size() - 1 - line] << "\n";
                                    }
                                    else
                                    {
                                        output << EHH_0_down[line - EHH_0_up.size()] << "\t" << EHH_1_down[line - EHH_0_up.size()] << "\n";
                                    }
                                }

                                output.close();

                                // clear global vectors
                                positions_Collect.clear();
                                EHH_0_up.clear();
                                EHH_1_up.clear();
                                EHH_0_down.clear();
                                EHH_1_down.clear();
                            }
                            else
                            {
                                cout << "ERROR SNP NOT FOUND. PLEASE CHECK THE REPOSITORY" << endl;
                            }
                        }
                        else
                        {
                            vector<string> folder_index_Positions;

                            function.split(folder_index_Positions, folder_Index[0].first, '_');
                            int low_Value = stoi(folder_index_Positions[0]);

                            int min_Limit = position_SNP - this->default_SNP_BP_count;
                            // cout << min_Limit << "\t" << low_Value << endl;

                            if (min_Limit < low_Value)
                            {
                                min_Limit = low_Value;
                            }

                            function.split(folder_index_Positions, folder_Index[folder_Index.size() - 1].first, '_');
                            int high_Value = stoi(folder_index_Positions[1]);

                            int max_Limit = position_SNP + this->default_SNP_BP_count;

                            // cout << max_Limit << "\t" << high_Value << endl;

                            if (max_Limit > high_Value)
                            {
                                max_Limit = high_Value;
                            }

                            // cout << min_Limit << "\t" << low_Value << endl;
                            // cout << max_Limit << "\t" << high_Value << endl;

                            vector<string> collect_SNP_file;
                            vector<string> SNP_positions;

                            if (SNP_file.is_open())
                            {
                                string line;
                                getline(SNP_file, line); // skip first header line
                                while (getline(SNP_file, line))
                                {
                                    vector<string> positions;
                                    function.split_getPos_ONLY(positions, line, '\t');
                                    int pos = stoi(positions[1]);

                                    if (pos >= min_Limit && pos <= max_Limit)
                                    {
                                        collect_SNP_file.push_back(line);
                                        int displacement = pos - position_SNP;
                                        SNP_positions.push_back(to_string(displacement) + "\t" + to_string(pos));

                                        if (pos == position_SNP)
                                        {
                                            cout << "SNP found in repository: " << SNP_file_Name << endl;
                                            SNP_Index_in_file = collect_SNP_file.size() - 1;
                                        }
                                    }
                                    else if (pos > max_Limit)
                                    {
                                        break;
                                    }
                                }
                                SNP_file.close();
                            }

                            if (SNP_Index_in_file != -1)
                            {
                                // Collet remaining files
                                vector<string> file_List;

                                cout << "System is retrieving file(s)" << endl;
                                if (folder_Index.size() > 1)
                                {
                                    file_List = function.compound_interpolationSearch_ordered(folder_Index, min_Limit, max_Limit);
                                }
                                else
                                {
                                    file_List.push_back(folder_Index[0].second);
                                }
                                cout << "System has retrieved all file(s)" << endl;

                                vector<string> seg_Sites_ALL;
                                vector<string> seg_Sites_ALL_Positions;
                                int SNP_Index_in_FULL = -1;

                                //cout << file_List.size();

                                for (string files : file_List)
                                {
                                    // cout << files << "\t" << SNP_file_Name << endl;
                                    if (files == SNP_file_Name)
                                    {
                                        // Fill the SNP file data from what was collected
                                        // get where position is NOW
                                        SNP_Index_in_FULL = SNP_Index_in_file + seg_Sites_ALL.size();

                                        for (int lines = 0; lines < collect_SNP_file.size(); lines++)
                                        {
                                            seg_Sites_ALL.push_back(collect_SNP_file[lines]);
                                            seg_Sites_ALL_Positions.push_back(SNP_positions[lines]);
                                        }

                                        collect_SNP_file.clear();
                                        SNP_positions.clear();
                                    }
                                    else
                                    {
                                        fstream file;
                                        file.open(files, ios::in);
                                        if (file.is_open())
                                        {
                                            string line;
                                            getline(file, line); // skip header
                                            while (getline(file, line))
                                            {
                                                vector<string> positions;
                                                function.split_getPos_ONLY(positions, line, '\t');
                                                int pos = stoi(positions[1]);

                                                if (pos >= min_Limit && pos <= max_Limit)
                                                {
                                                    seg_Sites_ALL.push_back(line);
                                                    int displacement = pos - position_SNP;
                                                    seg_Sites_ALL_Positions.push_back(to_string(displacement) + "\t" + to_string(pos));
                                                }
                                                else if (pos > max_Limit)
                                                {
                                                    break;
                                                }
                                            }
                                            file.close();
                                        }
                                    }
                                }
                                // PROCESS EHH_SNP_Distance
                                // cout << seg_Sites_ALL_Positions[SNP_Index_in_FULL] << endl;
                                // cout << SNP_Index_in_FULL << "\t" << seg_Sites_ALL_Positions.size();
                                process_SNP_EHH_BP(seg_Sites_ALL, SNP_Index_in_FULL, N, SNP_Index_in_FULL, (seg_Sites_ALL.size() - SNP_Index_in_FULL - 1));

                                // print the results
                                string results_File = output_File + "/" +
                                                      gene_Name + "_" + Core_coordinates[0] + "_" + Core_coordinates[1] + "_" + to_string(this->default_SNP_BP_count) + "_BP.ehh";
                                function.createFile(results_File, "Displacement\tSNP_position\tEHH_0\tEHH_1");

                                cout << "\nWriting results to file: " << results_File << endl
                                     << endl;

                                output.open(results_File, ios::app);

                                for (int line = 0; line < seg_Sites_ALL_Positions.size(); line++)
                                {
                                    output << seg_Sites_ALL_Positions[line] << "\t";
                                    if (line < EHH_0_up.size())
                                    {
                                        output << EHH_0_up[EHH_0_up.size() - 1 - line] << "\t" << EHH_1_up[EHH_1_up.size() - 1 - line] << "\n";
                                    }
                                    else
                                    {
                                        output << EHH_0_down[line - EHH_0_up.size()] << "\t" << EHH_1_down[line - EHH_0_up.size()] << "\n";
                                    }
                                }

                                output.close();

                                // clear global vectors
                                // positions_Collect.clear();

                                EHH_0_up.clear();
                                EHH_1_up.clear();
                                EHH_0_down.clear();
                                EHH_1_down.clear();

                                // clear global arrays
                            }
                            else
                            {
                                cout << "ERROR SNP NOT FOUND. PLEASE CHECK THE REPOSITORY" << endl;
                            }
                        }
                    }
                    else
                    {
                        cout << "ERROR SNP NOT FOUND. PLEASE CHECK THE REPOSITORY" << endl;
                    }
                }
                else
                {
                    cout << "Gene name\t: " << gene_Name << endl;
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
                            function.split(Ext_coordinates, split_Data[2], ':');
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

                        cout << "System is collecting segregating site(s)" << endl;

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
                                    function.split_getPos_ONLY(positions, line, '\t');
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
                            // cout << core_Haplotype_Collection[core] << endl;
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

__global__ void cuda_SNP_grid_0_1_BP(int total_Segs, char *sites, int *index, char *Hap_array, int core_SNP_index, int *core_SNP_alleles, int *SNP_counts)
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
            // if (core_SNP_index == tid)
            // {
            //     printf("%c\n", sites[i]);
            // }
            i++;
        }

        int start_Hap = tid;
        int stride = 0;

        // int grid_Row = tid;
        int grid_Column = 0;

        if (core_SNP_index == tid)
        {
            int Hap_0 = 0;
            int Hap_1 = 0;
            while (i < site_End)
            {
                if (sites[i] == '0' || sites[i] == '1')
                {
                    char value = sites[i];
                    Hap_array[start_Hap + stride] = value;
                    stride = stride + total_Segs;
                    if (sites[i] == '0')
                    {
                        core_SNP_alleles[grid_Column] = 0;
                        Hap_0 = Hap_0 + 1;
                    }
                    else
                    {
                        core_SNP_alleles[grid_Column] = 1;
                        Hap_1 = Hap_1 + 1;
                    }
                    grid_Column++;
                }
                i++;
            }
            SNP_counts[0] = Hap_0;
            SNP_counts[1] = Hap_1;
        }
        else
        {
            while (i < site_End)
            {
                if (sites[i] == '0' || sites[i] == '1')
                {
                    char value = sites[i];
                    Hap_array[start_Hap + stride] = value;
                    stride = stride + total_Segs;
                }
                i++;
            }
        }
        tid += blockDim.x * gridDim.x;
    }
}

__global__ void cuda_SNP_grid_0_1(int total_Segs, char *sites, int *index, char *Hap_array, int core_SNP_index, int *core_SNP_alleles, int *SNP_counts, int *cuda_pos_start_Index, int *cuda_pos_end_Index)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < total_Segs)
    {
        int column = 0;
        int site_Start = index[tid];
        int site_End = index[tid + 1];

        int i = site_Start;

        while (column < 1)
        {
            if (sites[i] == '\t')
            {
                column++;
            }
            i++;
        }

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

        while (column < 9)
        {
            if (sites[i] == '\t')
            {
                column++;
            }
            // if (core_SNP_index == tid)
            // {
            //     printf("%c\n", sites[i]);
            // }
            i++;
        }

        int start_Hap = tid;
        int stride = 0;

        // int grid_Row = tid;
        int grid_Column = 0;

        if (core_SNP_index == tid)
        {
            int Hap_0 = 0;
            int Hap_1 = 0;
            while (i < site_End)
            {
                if (sites[i] == '0' || sites[i] == '1')
                {
                    char value = sites[i];
                    Hap_array[start_Hap + stride] = value;
                    stride = stride + total_Segs;
                    if (sites[i] == '0')
                    {
                        core_SNP_alleles[grid_Column] = 0;
                        Hap_0 = Hap_0 + 1;
                    }
                    else
                    {
                        core_SNP_alleles[grid_Column] = 1;
                        Hap_1 = Hap_1 + 1;
                    }
                    grid_Column++;
                }
                i++;
            }
            SNP_counts[0] = Hap_0;
            SNP_counts[1] = Hap_1;
        }
        else
        {
            while (i < site_End)
            {
                if (sites[i] == '0' || sites[i] == '1')
                {
                    char value = sites[i];
                    Hap_array[start_Hap + stride] = value;
                    stride = stride + total_Segs;
                }
                i++;
            }
        }
        tid += blockDim.x * gridDim.x;
    }
}

__global__ void cuda_EHH_up_0(int total_Segs_UP, char **grid, int core_SNP_index, int *core_SNP_alleles, int N, int **Indexes_found_Zero, int zero_Count, float *EHH_Zero, int combo_Zero)
{
    // ! TEST CODE, NOT USED IN THE ACTUAL PROGRAMME
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < total_Segs_UP)
    {
        int augment = tid + 1;
        int check_until_row = core_SNP_index - augment;

        // printf("%d\n", check_until_row);

        int found_column = 0;
        int found_row = tid;

        int sum_0_Combos = 0;

        for (int query_Hap = 0; query_Hap < N; query_Hap++)
        {
            if (found_column < zero_Count)
            {
                if (core_SNP_alleles[query_Hap] == 0)
                {
                    int found_YES_NO_query = 0;
                    for (size_t check_Found = 0; check_Found < found_column; check_Found++)
                    {
                        if (Indexes_found_Zero[found_row][check_Found] == query_Hap)
                        {
                            found_YES_NO_query = 1;
                            break;
                        }
                    }
                    if (found_YES_NO_query == 0)
                    {
                        Indexes_found_Zero[found_row][found_column] = query_Hap;
                        int query_Count = 1;
                        found_column++;
                        for (int subject_Hap = query_Hap + 1; subject_Hap < N; subject_Hap++)
                        {
                            if (core_SNP_alleles[subject_Hap] == 0)
                            {
                                int found_YES_NO_subject = 0;
                                for (size_t check_Found = 0; check_Found < found_column; check_Found++)
                                {
                                    if (Indexes_found_Zero[found_row][check_Found] == subject_Hap)
                                    {
                                        found_YES_NO_subject = 1;
                                        break;
                                    }
                                }
                                if (found_YES_NO_subject == 0)
                                {
                                    int match = 0;
                                    for (int grid_row = core_SNP_index - 1; grid_row >= check_until_row; grid_row--)
                                    {
                                        // printf("%d %d %d\n", tid, grid_row, check_until_row);
                                        if (grid[grid_row][query_Hap] != grid[grid_row][subject_Hap])
                                        {
                                            match = 1;
                                            break;
                                        }
                                    }
                                    if (match == 0)
                                    {
                                        // same hap found
                                        query_Count++;
                                        Indexes_found_Zero[found_row][found_column] = subject_Hap;
                                        found_column++;
                                    }
                                }
                            }
                        }
                        // get_combo
                        int combo = (query_Count * (query_Count - 1)) / 2;
                        sum_0_Combos = sum_0_Combos + combo;
                    }
                }
            }
            else
            {
                break;
            }
        }

        EHH_Zero[tid] = (float)sum_0_Combos / (float)combo_Zero;
        // printf("%d\n", sum_0_Combos);
        //  EHH_Zero[tid] = (float)sum_0_Combos;

        tid += blockDim.x * gridDim.x;
    }
}

void ehh::process_SNP_EHH_BP(vector<string> &total_Segregrating_sites, int SNP_Index_in_FULL, int N, int SNPs_above, int SNPs_below)
{
    functions function = functions();

    cout << "\nSystem is calculating EHH: \n"
         << endl;

    cout << "STEP 1 of 4: Organizing SNPs for GPU" << endl;

    int num_segregrating_Sites = total_Segregrating_sites.size();
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

    int *core_SNP_alleles;
    int *cuda_core_SNP_alleles;
    int *SNP_counts, *cuda_SNP_counts;
    core_SNP_alleles = (int *)malloc(N * sizeof(int));
    SNP_counts = (int *)malloc(2 * sizeof(int));
    cudaMallocManaged(&cuda_SNP_counts, 2 * sizeof(int));
    cudaMallocManaged(&cuda_core_SNP_alleles, N * sizeof(int));

    cudaMemcpy(cuda_full_Char, full_Char, (Seg_sites.size() + 1) * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(cuda_site_Index, site_Index, (num_segregrating_Sites + 1) * sizeof(int), cudaMemcpyHostToDevice);

    /**
     * @param cuda_Hap_array stores the forged Haplotypes for the region under study.
     * @param Hap_array is used by the CPU. Is a COPY of cuda_Hap_array.
     **/
    char *Hap_array, *cuda_Hap_array;
    cudaMallocManaged(&cuda_Hap_array, ((N * num_segregrating_Sites) + 1) * sizeof(char));

    cout << "STEP 2 OF 4: Haplotype forging from segregrating sites" << endl;
    // cuda_SNP_grid_0_1_BP(int total_Segs, char *sites, int *index, char *Hap_array, int core_SNP_index, int *core_SNP_alleles, int *SNP_counts)
    cuda_SNP_grid_0_1_BP<<<tot_Blocks, tot_ThreadsperBlock>>>(num_segregrating_Sites, cuda_full_Char, cuda_site_Index, cuda_Hap_array, SNP_Index_in_FULL, cuda_core_SNP_alleles, cuda_SNP_counts);

    cudaError_t err = cudaGetLastError();

    if (err != cudaSuccess)
    {
        printf("CUDA Error: %s\n", cudaGetErrorString(err));

        // Possibly: exit(-1) if program cannot continue....
    }
    cudaDeviceSynchronize();

    cudaMemcpy(core_SNP_alleles, cuda_core_SNP_alleles, N * sizeof(int), cudaMemcpyDeviceToHost);

    Hap_array = (char *)malloc(((N * num_segregrating_Sites) + 1) * sizeof(char));
    cudaMemcpy(Hap_array, cuda_Hap_array, ((N * num_segregrating_Sites) + 1) * sizeof(char), cudaMemcpyDeviceToHost);

    cudaMemcpy(SNP_counts, cuda_SNP_counts, 2 * sizeof(int), cudaMemcpyDeviceToHost);

    cudaFree(cuda_full_Char);
    cudaFree(cuda_site_Index);
    cudaFree(cuda_core_SNP_alleles);
    cudaFree(cuda_SNP_counts);
    cudaFree(cuda_Hap_array);

    cout << "STEP 3 OF 4: Collecting haplotypes" << endl;

    string haplotypes(Hap_array);
    vector<string> Haplotypes_0, Haplotypes_1;

    int N_count = 0;
    for (int i = 0; i < (num_segregrating_Sites * N); i = i + num_segregrating_Sites)
    {
        // cout << core_SNP_alleles[N_count];
        if (core_SNP_alleles[N_count] == 0)
        {
            Haplotypes_0.push_back(haplotypes.substr(i, num_segregrating_Sites));
        }
        else
        {
            Haplotypes_1.push_back(haplotypes.substr(i, num_segregrating_Sites));
        }

        N_count++;
    }

    cout << "STEP 4 OF 4: Calculating EHH" << endl;

    // cout << Haplotypes_0[0].size() << endl;

    int combo_Zero = (SNP_counts[0] * (SNP_counts[0] - 1)) / 2;
    int combo_One = (SNP_counts[1] * (SNP_counts[1] - 1)) / 2;

    // cout << SNPs_above << endl;
    // cout << SNPs_below << endl;
    // exit(3);

    for (size_t i = 0; i < num_segregrating_Sites; i++)
    {
        if (i < SNPs_above)
        {
            EHH_0_up.push_back("");
            EHH_1_up.push_back("");
        }
        else
        {
            EHH_0_down.push_back("");
            EHH_1_down.push_back("");
        }
    }

    vector<thread> threads_vec;

    cout << "             Calculating upper bound EHH 0" << endl;

    int SNPs_in_thread = SNPs_above / CPU_cores;
    int remainder = SNPs_above % CPU_cores;

    for (size_t i = 0; i < CPU_cores; i++)
    {
        int start = i * SNPs_in_thread;
        threads_vec.push_back(thread{&ehh::calc_EHH_0_1_up, this, Haplotypes_0, SNP_Index_in_FULL, start, start + SNPs_in_thread, combo_Zero, 0});
    }

    for (thread &t : threads_vec)
    {
        if (t.joinable())
        {
            t.join();
        }
    }

    threads_vec.clear();

    // cout << "             Calculating upper bound EHH 0" << endl;

    if (remainder != 0)
    {
        int start = SNPs_above - remainder;
        threads_vec.push_back(thread{&ehh::calc_EHH_0_1_up, this, Haplotypes_0, SNP_Index_in_FULL, start, SNPs_above, combo_Zero, 0});

        for (thread &t : threads_vec)
        {
            if (t.joinable())
            {
                t.join();
            }
        }

        threads_vec.clear();
    }

    cout << "             Calculating upper bound EHH 1" << endl;

    for (size_t i = 0; i < CPU_cores; i++)
    {
        int start = i * SNPs_in_thread;
        threads_vec.push_back(thread{&ehh::calc_EHH_0_1_up, this, Haplotypes_1, SNP_Index_in_FULL, start, start + SNPs_in_thread, combo_One, 1});
    }

    for (thread &t : threads_vec)
    {
        if (t.joinable())
        {
            t.join();
        }
    }

    threads_vec.clear();

    if (remainder != 0)
    {
        int start = SNPs_above - remainder;
        threads_vec.push_back(thread{&ehh::calc_EHH_0_1_up, this, Haplotypes_1, SNP_Index_in_FULL, start, SNPs_above, combo_One, 1});

        for (thread &t : threads_vec)
        {
            if (t.joinable())
            {
                t.join();
            }
        }

        threads_vec.clear();
    }

    cout << "             Calculating lower bound EHH 0" << endl;

    SNPs_in_thread = (SNPs_below + 1) / CPU_cores;
    remainder = (SNPs_below + 1) % CPU_cores;

    for (size_t i = 0; i < CPU_cores; i++)
    {
        int start = i * SNPs_in_thread;
        threads_vec.push_back(thread{&ehh::calc_EHH_0_1_down, this, Haplotypes_0, SNP_Index_in_FULL, start, start + SNPs_in_thread, combo_Zero, 0});
    }

    for (thread &t : threads_vec)
    {
        if (t.joinable())
        {
            t.join();
        }
    }

    threads_vec.clear();

    if (remainder != 0)
    {
        int start = (SNPs_below + 1) - remainder;
        threads_vec.push_back(thread{&ehh::calc_EHH_0_1_down, this, Haplotypes_0, SNP_Index_in_FULL, start, SNPs_below + 1, combo_Zero, 0});

        for (thread &t : threads_vec)
        {
            if (t.joinable())
            {
                t.join();
            }
        }

        threads_vec.clear();
    }

    cout << "             Calculating lower bound EHH 1" << endl;

    for (size_t i = 0; i < CPU_cores; i++)
    {
        int start = i * SNPs_in_thread;
        threads_vec.push_back(thread{&ehh::calc_EHH_0_1_down, this, Haplotypes_1, SNP_Index_in_FULL, start, start + SNPs_in_thread, combo_One, 1});
    }

    for (thread &t : threads_vec)
    {
        if (t.joinable())
        {
            t.join();
        }
    }

    threads_vec.clear();

    if (remainder != 0)
    {
        int start = (SNPs_below + 1) - remainder;
        threads_vec.push_back(thread{&ehh::calc_EHH_0_1_down, this, Haplotypes_1, SNP_Index_in_FULL, start, SNPs_below + 1, combo_One, 1});

        for (thread &t : threads_vec)
        {
            if (t.joinable())
            {
                t.join();
            }
        }

        threads_vec.clear();
    }

    free(full_Char);
    free(Hap_array);
    free(SNP_counts);
}

void ehh::process_SNP_EHH(int position, vector<pair<string, string>> &folder_Index, int segment_Position, int SNP_Index_in_file, vector<string> &segment_Segregrating_sites, float N, int num_top, int num_bottom)
{
    functions function = functions();

    cout << "\nSystem is initiating EHH engine: \n"
         << endl;

    cout << "STEP 1 of 6: Completing SNP collection" << endl;

    vector<vector<string>> top_Segregrating_sites_Total;
    vector<string> bottom_Segregrating_sites;

    int top_collect = num_top;
    int bottom_collect = num_bottom;

    int top_Collected = 0;

    if (top_collect < this->default_SNP_count)
    {
        cout << "             Collecting upper bound SNPs" << endl;
        int top_file_Index = segment_Position;

        while (top_collect < (default_SNP_count + 1))
        {
            top_file_Index = top_file_Index - 1;
            if (top_file_Index >= 0)
            {
                string SNP_file_Name = folder_Index[top_file_Index].second;
                fstream SNP_file;
                SNP_file.open(SNP_file_Name, ios::in);
                if (SNP_file.is_open())
                {
                    vector<string> top_Segregrating_sites;
                    string line;
                    getline(SNP_file, line); // skip first header line
                    while (getline(SNP_file, line))
                    {
                        top_Segregrating_sites.push_back(line);
                    }
                    top_Segregrating_sites_Total.push_back(top_Segregrating_sites);
                    top_Collected = top_Collected + top_Segregrating_sites.size();
                    top_collect = top_collect + top_Segregrating_sites.size();
                    SNP_file.close();
                }
            }
            else
            {
                break;
            }
        }
    }

    if (bottom_collect < this->default_SNP_count)
    {
        cout << "             Collecting lower bound SNPs" << endl;
        int bottom_file_Index = segment_Position;
        while (bottom_collect != default_SNP_count)
        {
            bottom_file_Index = bottom_file_Index + 1;
            if (bottom_file_Index < folder_Index.size())
            {
                string SNP_file_Name = folder_Index[bottom_file_Index].second;
                fstream SNP_file;
                SNP_file.open(SNP_file_Name, ios::in);
                if (SNP_file.is_open())
                {
                    string line;
                    getline(SNP_file, line); // skip first header line
                    while (getline(SNP_file, line))
                    {
                        bottom_Segregrating_sites.push_back(line);
                        bottom_collect = bottom_collect + 1;
                        if (bottom_collect == default_SNP_count)
                        {
                            break;
                        }
                    }
                    SNP_file.close();
                }
            }
            else
            {
                break;
            }
            // bottom_collect = bottom_collect + bottom_Segregrating_sites.size();
        }
    }

    cout << "STEP 2 of 6: Organizing SNPs for GPU" << endl;

    int num_segregrating_Sites;

    if (top_collect > this->default_SNP_count)
    {
        num_segregrating_Sites = this->default_SNP_count;
    }
    else
    {
        num_segregrating_Sites = top_collect;
    }

    num_segregrating_Sites = num_segregrating_Sites + 1 + bottom_collect;

    int stride_track = 0;

    string Seg_sites = "";
    int site_Index[num_segregrating_Sites + 1];
    site_Index[0] = 0;

    // fill top bit
    if (num_top >= this->default_SNP_count)
    {
        for (size_t i = (SNP_Index_in_file - this->default_SNP_count); i < SNP_Index_in_file; i++)
        {
            Seg_sites.append(segment_Segregrating_sites[i]);
            site_Index[stride_track + 1] = site_Index[stride_track] + segment_Segregrating_sites[i].size();
            stride_track++;
        }
    }
    else
    {
        // fill remaining top from collected

        int remaining_get = this->default_SNP_count - num_top;
        if (remaining_get > top_Collected)
        {
            // remaining_get = top_Collected;
            if (top_Collected > 0)
            {
                int files = top_Segregrating_sites_Total.size();
                for (int file_current = files - 1; file_current >= 0; file_current--)
                {
                    vector<string> top_Segregrating_sites = top_Segregrating_sites_Total[file_current];
                    for (int i = 0; i < top_Segregrating_sites.size(); i++)
                    {
                        Seg_sites.append(top_Segregrating_sites[i]);
                        site_Index[stride_track + 1] = site_Index[stride_track] + top_Segregrating_sites[i].size();
                        stride_track++;
                    }
                }
            }
            for (int i = 0; i < SNP_Index_in_file; i++)
            {
                Seg_sites.append(segment_Segregrating_sites[i]);
                site_Index[stride_track + 1] = site_Index[stride_track] + segment_Segregrating_sites[i].size();
                stride_track++;
            }
        }
        else
        {
            int files = top_Segregrating_sites_Total.size();
            int file_Num, line_Num;

            for (int i = 0; i < files; i++)
            {
                int Num_seg_Sites = top_Segregrating_sites_Total[i].size();
                if (remaining_get <= Num_seg_Sites)
                {
                    file_Num = i;
                    line_Num = Num_seg_Sites - remaining_get;
                    break;
                }
                remaining_get = remaining_get - Num_seg_Sites;
            }

            vector<string> partial_top_Segregrating_sites = top_Segregrating_sites_Total[file_Num];

            for (int i = line_Num; i < partial_top_Segregrating_sites.size(); i++)
            {
                Seg_sites.append(partial_top_Segregrating_sites[i]);
                site_Index[stride_track + 1] = site_Index[stride_track] + partial_top_Segregrating_sites[i].size();
                stride_track++;
            }

            for (int file_current = file_Num - 1; file_current >= 0; file_current--)
            {
                vector<string> top_Segregrating_sites = top_Segregrating_sites_Total[file_current];
                for (size_t i = 0; i < top_Segregrating_sites.size(); i++)
                {
                    Seg_sites.append(top_Segregrating_sites[i]);
                    site_Index[stride_track + 1] = site_Index[stride_track] + top_Segregrating_sites[i].size();
                    stride_track++;
                }
            }

            for (int i = 0; i < SNP_Index_in_file; i++)
            {
                Seg_sites.append(segment_Segregrating_sites[i]);
                site_Index[stride_track + 1] = site_Index[stride_track] + segment_Segregrating_sites[i].size();
                stride_track++;
            }
        }
    }

    int SNP_position_in_full = stride_track;
    Seg_sites.append(segment_Segregrating_sites[SNP_Index_in_file]);
    site_Index[stride_track + 1] = site_Index[stride_track] + segment_Segregrating_sites[SNP_Index_in_file].size();
    stride_track++;
    // fill bottom bit

    if (num_bottom > this->default_SNP_count)
    {
        for (int i = SNP_Index_in_file + 1; i < (SNP_Index_in_file + this->default_SNP_count); i++)
        {
            Seg_sites.append(segment_Segregrating_sites[i]);
            site_Index[stride_track + 1] = site_Index[stride_track] + segment_Segregrating_sites[i].size();
            stride_track++;
        }
    }
    else
    {
        for (int i = SNP_Index_in_file + 1; i < segment_Segregrating_sites.size(); i++)
        {
            Seg_sites.append(segment_Segregrating_sites[i]);
            site_Index[stride_track + 1] = site_Index[stride_track] + segment_Segregrating_sites[i].size();
            stride_track++;
        }

        for (int i = 0; i < bottom_Segregrating_sites.size(); i++)
        {
            Seg_sites.append(bottom_Segregrating_sites[i]);
            site_Index[stride_track + 1] = site_Index[stride_track] + bottom_Segregrating_sites[i].size();
            stride_track++;
        }
    }

    char *full_Char;
    full_Char = (char *)malloc((Seg_sites.size() + 1) * sizeof(char));
    strcpy(full_Char, Seg_sites.c_str());

    top_Segregrating_sites_Total.clear();
    segment_Segregrating_sites.clear();
    bottom_Segregrating_sites.clear();

    char *cuda_full_Char;
    cudaMallocManaged(&cuda_full_Char, (Seg_sites.size() + 1) * sizeof(char));
    int *cuda_site_Index;
    cudaMallocManaged(&cuda_site_Index, (num_segregrating_Sites + 1) * sizeof(int));

    int *core_SNP_alleles;
    int *cuda_core_SNP_alleles;
    int *SNP_counts, *cuda_SNP_counts;
    core_SNP_alleles = (int *)malloc(N * sizeof(int));
    SNP_counts = (int *)malloc(2 * sizeof(int));
    cudaMallocManaged(&cuda_SNP_counts, 2 * sizeof(int));
    cudaMallocManaged(&cuda_core_SNP_alleles, N * sizeof(int));

    cudaMemcpy(cuda_full_Char, full_Char, (Seg_sites.size() + 1) * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(cuda_site_Index, site_Index, (num_segregrating_Sites + 1) * sizeof(int), cudaMemcpyHostToDevice);

    /**
     * @param cuda_Hap_array stores the forged Haplotypes for the region under study.
     * @param Hap_array is used by the CPU. Is a COPY of cuda_Hap_array.
     **/
    char *Hap_array, *cuda_Hap_array;
    cudaMallocManaged(&cuda_Hap_array, ((N * num_segregrating_Sites) + 1) * sizeof(char));

    int *pos_start_Index, *pos_end_Index, *cuda_pos_start_Index, *cuda_pos_end_Index;
    cudaMallocManaged(&cuda_pos_start_Index, num_segregrating_Sites * sizeof(int));
    cudaMallocManaged(&cuda_pos_end_Index, num_segregrating_Sites * sizeof(int));
    pos_start_Index = (int *)malloc(num_segregrating_Sites * sizeof(int));
    pos_end_Index = (int *)malloc(num_segregrating_Sites * sizeof(int));

    // char **cuda_snp_N_grid;
    // cudaMallocManaged(&cuda_snp_N_grid, N * num_segregrating_Sites * sizeof(char));
    // char **tmp = (char **)malloc(num_segregrating_Sites * sizeof(tmp[0]));
    // for (int i = 0; i < num_segregrating_Sites; i++)
    // {
    //     cudaMalloc((void **)&tmp[i], N * sizeof(tmp[0][0]));
    // }
    // cudaMemcpy(cuda_snp_N_grid, tmp, num_segregrating_Sites * sizeof(char *), cudaMemcpyHostToDevice);
    // free(tmp);

    cout << "STEP 3 OF 6: Haplotype forging from segregrating sites" << endl;
    cuda_SNP_grid_0_1<<<tot_Blocks, tot_ThreadsperBlock>>>(num_segregrating_Sites, cuda_full_Char, cuda_site_Index, cuda_Hap_array, SNP_position_in_full, cuda_core_SNP_alleles, cuda_SNP_counts, cuda_pos_start_Index, cuda_pos_end_Index);

    cudaError_t err = cudaGetLastError();

    if (err != cudaSuccess)
    {
        printf("CUDA Error: %s\n", cudaGetErrorString(err));

        // Possibly: exit(-1) if program cannot continue....
    }
    cudaDeviceSynchronize();

    cudaMemcpy(core_SNP_alleles, cuda_core_SNP_alleles, N * sizeof(int), cudaMemcpyDeviceToHost);

    Hap_array = (char *)malloc(((N * num_segregrating_Sites) + 1) * sizeof(char));
    cudaMemcpy(Hap_array, cuda_Hap_array, ((N * num_segregrating_Sites) + 1) * sizeof(char), cudaMemcpyDeviceToHost);

    cudaMemcpy(SNP_counts, cuda_SNP_counts, 2 * sizeof(int), cudaMemcpyDeviceToHost);

    cudaMemcpy(pos_start_Index, cuda_pos_start_Index, num_segregrating_Sites * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(pos_end_Index, cuda_pos_end_Index, num_segregrating_Sites * sizeof(int), cudaMemcpyDeviceToHost);

    cudaFree(cuda_full_Char);
    cudaFree(cuda_site_Index);
    cudaFree(cuda_core_SNP_alleles);
    cudaFree(cuda_SNP_counts);
    cudaFree(cuda_pos_start_Index);
    cudaFree(cuda_pos_end_Index);
    cudaFree(cuda_Hap_array);

    // free(site_Index);

    // get positions of 0 and 1 haps (no need conversion).
    // store 0 ex haps in one vector and 1 in another.
    // spawn arrays for all in 0 and in 1 till end.

    cout << "STEP 4 OF 6: Collecting haplotypes" << endl;

    string haplotypes(Hap_array);
    vector<string> Haplotypes_0, Haplotypes_1;

    int SNPs_above = SNP_position_in_full;
    int SNPs_below = num_segregrating_Sites - SNP_position_in_full - 1;

    int N_count = 0;
    for (int i = 0; i < (num_segregrating_Sites * N); i = i + num_segregrating_Sites)
    {
        // cout << core_SNP_alleles[N_count];
        if (core_SNP_alleles[N_count] == 0)
        {
            Haplotypes_0.push_back(haplotypes.substr(i, num_segregrating_Sites));
        }
        else
        {
            Haplotypes_1.push_back(haplotypes.substr(i, num_segregrating_Sites));
        }

        N_count++;

        //  initialise global multithread vectors.

        // if (N_count < SNPs_above)
        // {
        //     EHH_0_up.push_back("-1");
        //     EHH_1_up.push_back("-1");
        // }
        // else
        // {
        //     EHH_0_down.push_back("-1");
        //     EHH_1_down.push_back("-1");
        // }
        // positions_Collect.push_back("-1");
    }
    // cout << "\n\n";
    // cout << Haplotypes_0[0].size() << endl
    //      << SNP_position_in_full << endl;

    // cout << SNPs_above << endl
    //      << SNPs_below << endl;

    cout << "STEP 5 OF 6: Calculating EHH" << endl;

    // cout << "0 " << SNP_counts[0] << endl;
    // cout << "1 " << SNP_counts[1] << " " << Haplotypes_1.size() << endl
    //      << endl;

    int combo_Zero = (SNP_counts[0] * (SNP_counts[0] - 1)) / 2;
    int combo_One = (SNP_counts[1] * (SNP_counts[1] - 1)) / 2;

    // int **cuda_Indexes_found_Zero;
    // cudaMallocManaged(&cuda_Indexes_found_Zero, SNP_counts[0] * SNPs_above * sizeof(int));
    // int **tmp_2 = (int **)malloc(SNPs_above * sizeof(tmp_2[0]));
    // for (int i = 0; i < SNPs_above; i++)
    // {
    //     cudaMalloc((void **)&tmp_2[i], SNP_counts[0] * sizeof(tmp_2[0][0]));
    // }
    // cudaMemcpy(cuda_Indexes_found_Zero, tmp_2, SNPs_above * sizeof(int *), cudaMemcpyHostToDevice);
    // free(tmp_2);

    // float *cuda_EHH_Zero, *EHH_Zero;
    // cudaMallocManaged(&cuda_EHH_Zero, SNPs_above * sizeof(float));
    // EHH_Zero = (float *)malloc(SNPs_above * sizeof(float));

    // cuda_EHH_up_0(int total_Segs_UP, char **grid, int core_SNP_index, int *core_SNP_alleles, int N, int **Indexes_found_Zero, int zero_Count, float *EHH_Zero, int combo_Zero)
    // cout << "STEP 5 OF X: Calculating EHH" << endl;
    // cuda_EHH_up_0<<<tot_Blocks, tot_ThreadsperBlock>>>(SNPs_above, cuda_snp_N_grid, SNP_position_in_full, cuda_core_SNP_alleles, (int)N, cuda_Indexes_found_Zero, SNP_counts[0], cuda_EHH_Zero, combo_Zero);
    // cudaError_t err2 = cudaGetLastError();

    // if (err2 != cudaSuccess)
    // {
    //     printf("CUDA Error: %s\n", cudaGetErrorString(err2));

    //     // Possibly: exit(-1) if program cannot continue....
    // }
    // cudaDeviceSynchronize();

    // cudaMemcpy(EHH_Zero, cuda_EHH_Zero, SNPs_above * sizeof(float), cudaMemcpyDeviceToHost);

    // for (size_t i = 0; i < SNPs_above; i++)
    // {
    //     cout << EHH_Zero[i] << endl;
    // }

    for (size_t i = 0; i < num_segregrating_Sites; i++)
    {
        if (i < SNPs_above)
        {
            EHH_0_up.push_back("");
            EHH_1_up.push_back("");
        }
        else
        {
            EHH_0_down.push_back("-1");
            EHH_1_down.push_back("-1");
        }
        // EHH_values.push_back("1");
        positions_Collect.push_back("-1");
    }

    // multithread

    vector<thread> threads_vec;

    cout << "             Calculating upper bound EHH 0" << endl;

    int SNPs_in_thread = SNPs_above / CPU_cores;
    int remainder = SNPs_above % CPU_cores;

    for (size_t i = 0; i < CPU_cores; i++)
    {
        int start = i * SNPs_in_thread;
        threads_vec.push_back(thread{&ehh::calc_EHH_0_1_up, this, Haplotypes_0, SNP_position_in_full, start, start + SNPs_in_thread, combo_Zero, 0});
    }

    for (thread &t : threads_vec)
    {
        if (t.joinable())
        {
            t.join();
        }
    }

    threads_vec.clear();

    if (remainder != 0)
    {
        int start = SNPs_above - remainder;
        threads_vec.push_back(thread{&ehh::calc_EHH_0_1_up, this, Haplotypes_0, SNP_position_in_full, start, SNPs_above, combo_Zero, 0});

        for (thread &t : threads_vec)
        {
            if (t.joinable())
            {
                t.join();
            }
        }

        threads_vec.clear();
    }

    cout << "             Calculating upper bound EHH 1" << endl;

    for (size_t i = 0; i < CPU_cores; i++)
    {
        int start = i * SNPs_in_thread;
        threads_vec.push_back(thread{&ehh::calc_EHH_0_1_up, this, Haplotypes_1, SNP_position_in_full, start, start + SNPs_in_thread, combo_One, 1});
    }

    for (thread &t : threads_vec)
    {
        if (t.joinable())
        {
            t.join();
        }
    }

    threads_vec.clear();

    if (remainder != 0)
    {
        int start = SNPs_above - remainder;
        threads_vec.push_back(thread{&ehh::calc_EHH_0_1_up, this, Haplotypes_1, SNP_position_in_full, start, SNPs_above, combo_One, 1});

        for (thread &t : threads_vec)
        {
            if (t.joinable())
            {
                t.join();
            }
        }

        threads_vec.clear();
    }

    cout << "             Calculating lower bound EHH 0" << endl;

    SNPs_in_thread = (SNPs_below + 1) / CPU_cores;
    remainder = (SNPs_below + 1) % CPU_cores;

    for (size_t i = 0; i < CPU_cores; i++)
    {
        int start = i * SNPs_in_thread;
        threads_vec.push_back(thread{&ehh::calc_EHH_0_1_down, this, Haplotypes_0, SNP_position_in_full, start, start + SNPs_in_thread, combo_Zero, 0});
    }

    for (thread &t : threads_vec)
    {
        if (t.joinable())
        {
            t.join();
        }
    }

    threads_vec.clear();

    if (remainder != 0)
    {
        int start = (SNPs_below + 1) - remainder;
        threads_vec.push_back(thread{&ehh::calc_EHH_0_1_down, this, Haplotypes_0, SNP_position_in_full, start, SNPs_below + 1, combo_Zero, 0});

        for (thread &t : threads_vec)
        {
            if (t.joinable())
            {
                t.join();
            }
        }

        threads_vec.clear();
    }

    cout << "             Calculating lower bound EHH 1" << endl;

    for (size_t i = 0; i < CPU_cores; i++)
    {
        int start = i * SNPs_in_thread;
        threads_vec.push_back(thread{&ehh::calc_EHH_0_1_down, this, Haplotypes_1, SNP_position_in_full, start, start + SNPs_in_thread, combo_One, 1});
    }

    for (thread &t : threads_vec)
    {
        if (t.joinable())
        {
            t.join();
        }
    }

    threads_vec.clear();

    if (remainder != 0)
    {
        int start = SNPs_below + 1 - remainder;
        threads_vec.push_back(thread{&ehh::calc_EHH_0_1_down, this, Haplotypes_1, SNP_position_in_full, start, SNPs_below + 1, combo_One, 1});

        for (thread &t : threads_vec)
        {
            if (t.joinable())
            {
                t.join();
            }
        }

        threads_vec.clear();
    }

    // exit(3);

    cout << "STEP 6 OF 6: Extracting positions\n"
         << endl;

    SNPs_in_thread = num_segregrating_Sites / CPU_cores;
    remainder = num_segregrating_Sites % CPU_cores;

    for (size_t i = 0; i < CPU_cores; i++)
    {
        int start = i * SNPs_in_thread;
        threads_vec.push_back(thread{&ehh::pos_Construct, this, full_Char, pos_start_Index, pos_end_Index, position, start, start + SNPs_in_thread});
    }

    for (thread &t : threads_vec)
    {
        if (t.joinable())
        {
            t.join();
        }
    }

    threads_vec.clear();

    if (remainder != 0)
    {
        int start = num_segregrating_Sites - CPU_cores + 1;
        threads_vec.push_back(thread{&ehh::pos_Construct, this, full_Char, pos_start_Index, pos_end_Index, position, start, num_segregrating_Sites});

        for (thread &t : threads_vec)
        {
            if (t.joinable())
            {
                t.join();
            }
        }
        threads_vec.clear();
    }

    // for (size_t i = 0; i < num_segregrating_Sites; i++)
    // {
    //     cout << positions_Collect[i] << "\n";
    // }

    free(full_Char);
    free(pos_start_Index);
    free(pos_end_Index);
    free(Hap_array);
    free(SNP_counts);

    // CLEAR 0 1 up down vectors
}

void ehh::calc_EHH_0_1_down(vector<string> Haplotypes, int position_of_core_SNP, int augment_start, int augment_stop, int combo_zero_one, int zero_one)
{
    vector<string> EHH_all;
    vector<int> positions_all;

    for (size_t i = augment_start; i < augment_stop; i++)
    {
        int augment = i;
        int pos, length;

        pos = position_of_core_SNP;
        length = augment + 1;

        int store_pos = augment;

        int found = 0;
        int sum_Extended = 0;

        vector<int> track_Found;

        for (int increment = 0; increment < Haplotypes.size(); increment++)
        {
            track_Found.push_back(-1);
        }

        for (int query = 0; query < Haplotypes.size(); query++)
        {
            // if (binary_search(found.begin(), found.end(), query) == false)
            if (track_Found[query] == -1)
            {
                string query_Hap = Haplotypes[query].substr(pos, length);

                int query_core_Count = 0;
                query_core_Count++;

                track_Found[query] = query;
                found++;
                // found.push_back(query);

                for (int subject = query + 1; subject < Haplotypes.size(); subject++)
                {
                    if (track_Found[subject] == -1)
                    // if (binary_search(found.begin(), found.end(), subject) == false)
                    {
                        string subject_Hap = Haplotypes[subject].substr(pos, length);
                        if (query_Hap.compare(subject_Hap) == 0)
                        {
                            query_core_Count++;

                            track_Found[subject] = subject;
                            found++;
                            // found.push_back(subject);
                            // sort(found.begin(), found.end());
                        }
                    }
                }
                int sum_Partial = (query_core_Count * (query_core_Count - 1)) / 2;
                // cout << store_pos << "\t" << sum_Partial << endl;
                sum_Extended = sum_Extended + sum_Partial;
            }

            if (found == Haplotypes.size())
            {
                break;
            }
            // else
            // {
            //     /**
            //      * The modified vector is now sorted to enable binary search.
            //      **/
            //     sort(found.begin(), found.end());
            // }
        }

        // get division
        float EHH = (float)sum_Extended / (float)combo_zero_one;
        if (isnan(EHH))
        {
            EHH_all.push_back("DIV_0");
        }
        else
        {
            EHH_all.push_back(to_string(EHH));
        }
        positions_all.push_back(store_pos);

        // cout << store_pos << "\t" << EHH << endl;
    }

    unique_lock<shared_mutex> ul(g_mutex);
    if (zero_one == 0)
    {
        for (size_t i = 0; i < positions_all.size(); i++)
        {
            EHH_0_down[positions_all[i]] = EHH_all[i];
        }
    }
    else
    {
        for (size_t i = 0; i < positions_all.size(); i++)
        {
            EHH_1_down[positions_all[i]] = EHH_all[i];
        }
    }
}

void ehh::calc_EHH_0_1_up(vector<string> Haplotypes, int position_of_core_SNP, int augment_start, int augment_stop, int combo_zero_one, int zero_one)
{
    vector<string> EHH_all;
    vector<int> positions_all;
    // cout << augment << endl;
    for (size_t i = augment_start; i < augment_stop; i++)
    {
        // cout << i << " " << position_of_core_SNP << endl
        //      << endl;
        int augment = i + 1;
        int pos, length;

        pos = position_of_core_SNP - augment;
        length = augment + 1;

        int store_pos = augment - 1;

        int found = 0;
        int sum_Extended = 0;

        vector<int> track_Found;

        for (int increment = 0; increment < Haplotypes.size(); increment++)
        {
            track_Found.push_back(-1);
        }

        for (int query = 0; query < Haplotypes.size(); query++)
        {
            // if (binary_search(found.begin(), found.end(), query) == false)
            if (track_Found[query] == -1)
            {
                string query_Hap = Haplotypes[query].substr(pos, length);

                int query_core_Count = 0;
                query_core_Count++;

                track_Found[query] = query;
                found++;
                // found.push_back(query);

                for (int subject = query + 1; subject < Haplotypes.size(); subject++)
                {
                    if (track_Found[subject] == -1)
                    // if (binary_search(found.begin(), found.end(), subject) == false)
                    {
                        string subject_Hap = Haplotypes[subject].substr(pos, length);
                        if (query_Hap.compare(subject_Hap) == 0)
                        {
                            query_core_Count++;

                            track_Found[subject] = subject;
                            found++;
                            // found.push_back(subject);
                            // sort(found.begin(), found.end());
                        }
                    }
                }
                int sum_Partial = (query_core_Count * (query_core_Count - 1)) / 2;
                // cout << store_pos << "\t" << sum_Partial << endl;
                sum_Extended = sum_Extended + sum_Partial;
            }

            if (found == Haplotypes.size())
            {
                break;
            }
            // else
            // {
            //     /**
            //      * The modified vector is now sorted to enable binary search.
            //      **/
            //     sort(found.begin(), found.end());
            // }
        }

        // get division
        float EHH = (float)sum_Extended / (float)combo_zero_one;
        if (isnan(EHH))
        {
            EHH_all.push_back("DIV_0");
        }
        else
        {
            EHH_all.push_back(to_string(EHH));
        }
        positions_all.push_back(store_pos);

        // cout << store_pos << "\t" << EHH << endl;
    }

    unique_lock<shared_mutex> ul(g_mutex);
    if (zero_one == 0)
    {
        for (size_t i = 0; i < positions_all.size(); i++)
        {
            EHH_0_up[positions_all[i]] = EHH_all[i];
        }
    }
    else
    {
        for (size_t i = 0; i < positions_all.size(); i++)
        {
            EHH_1_up[positions_all[i]] = EHH_all[i];
        }
    }
}

void ehh::pos_Construct(char *full_Char, int *pos_start_Index, int *pos_end_Index, int gen_position_of_Core, int start_pos_Array, int stop_pos_Array)
{
    vector<string> positions_data;
    vector<int> positions;

    for (int index = start_pos_Array; index < stop_pos_Array; index++)
    {
        /* code */

        string POS_string = "";

        for (int i = pos_start_Index[index]; i < pos_end_Index[index]; i++)
        {
            /**
             * Extracting the position through concatenation of the data in the Position column of the SNP data.
             **/
            POS_string = POS_string + full_Char[i];
        }

        int displacement = stoi(POS_string) - gen_position_of_Core;

        positions_data.push_back(to_string(displacement) + "\t" + POS_string);
        positions.push_back(index);
    }

    unique_lock<shared_mutex> ul(g_mutex);
    for (size_t i = 0; i < positions.size(); i++)
    {
        positions_Collect[positions[i]] = positions_data[i];
    }
}

// void ehh::calc_EHH_0_1(vector<string> Haplotypes, int position_of_core_SNP, char up_down, int augment, int one_zero, int combo_one_zero)
// {
//     // get pos of snp here itself

//     // int Seg_No;
//     int pos, length;

//     if (up_down == 'D')
//     {
//         // Seg_No = position_of_core_SNP + augment;
//         pos = position_of_core_SNP;
//         length = augment + 1;
//     }
//     else
//     {
//         // Seg_No = position_of_core_SNP - augment;
//         pos = position_of_core_SNP - augment;
//         length = augment + 1;
//     }

//     // do sort for binary sort
//     vector<int> found;
//     int sum_Extended = 0;

//     for (int query = 0; query < Haplotypes.size(); query++)
//     {
//         if (binary_search(found.begin(), found.end(), query) == false)
//         {
//             string query_Hap = Haplotypes[query].substr(pos, length);

//             int query_core_Count = 0;
//             query_core_Count++;

//             found.push_back(query);

//             for (int subject = query + 1; subject < Haplotypes.size(); subject++)
//             {
//                 if (find(found.begin(), found.end(), subject) == found.end())
//                 {
//                     string subject_Hap = Haplotypes[subject].substr(pos, length);
//                     if (query_Hap.compare(subject_Hap) == 0)
//                     {
//                         query_core_Count++;
//                         found.push_back(subject);
//                         // sort(found.begin(), found.end());
//                     }
//                 }
//             }

//             int sum_Partial = (query_core_Count * (query_core_Count - 1)) / 2;
//             sum_Extended = sum_Extended + sum_Partial;
//         }

//         if (found.size() == Haplotypes.size())
//         {
//             break;
//         }
//         else
//         {
//             /**
//              * The modified vector is now sorted to enable binary search.
//              **/
//             sort(found.begin(), found.end());
//         }
//     }

//     // get division
//     float EHH = (float)sum_Extended / (float)combo_one_zero;

//     cout << augment << "\t" << pos << "\t" << EHH << endl;
//     // int position_Store = -1;
//     // if (up_down == 'D')
//     // {
//     //     position_Store = pos + length;
//     // }
//     // else
//     // {
//     //     position_Store = pos;
//     // }

//     // unique_lock<shared_mutex> ul(g_mutex);
//     // EHH_values[position_Store] = EHH;
//     // if (one_zero == 0)
//     // {
//     //     if (up_down == 'D')
//     //     {
//     //         EHH_0_down[augment] = to_string(EHH);
//     //     }
//     //     else
//     //     {
//     //         EHH_0_up[pos] = to_string(EHH);
//     //     }
//     // }
//     // else
//     // {
//     //     if (up_down == 'D')
//     //     {
//     //         EHH_1_down[augment] = to_string(EHH);
//     //     }
//     //     else
//     //     {
//     //         EHH_1_up[pos] = to_string(EHH);
//     //     }
//     // }
// }

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

    // snp_N_grid = (int **)malloc(num_segregrating_Sites * sizeof(int *));**
    // for (int i = 0; i < num_segregrating_Sites; i++)
    // {
    //     snp_N_grid[i] = (int *)malloc(N * sizeof(int));
    // }

    cudaMallocManaged(&cuda_snp_N_grid, (N + 1) * num_segregrating_Sites * sizeof(char));
    char **tmp = (char **)malloc(num_segregrating_Sites * sizeof(tmp[0]));
    for (int i = 0; i < num_segregrating_Sites; i++)
    {
        cudaMalloc((void **)&tmp[i], (N + 1) * sizeof(tmp[0][0]));
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
    cout << "STEP 1 OF 3: Sorting segregating sites into an array" << endl;
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
    int found_Count = 0;

    for (size_t i = 0; i < core_Haplotypes_All.size(); i++)
    {
        found.push_back(-1);
    }

    cout << "STEP 3 OF 3: Detecting unique core haplotypes with their unique extended haplotypes" << endl;
    for (size_t query = 0; query < core_Haplotypes_All.size(); query++)
    {
        // int present_OR_not = 0;

        if (found[query] == -1)
        // if (binary_search(found.begin(), found.end(), query) == false)
        {
            int query_core_Count = 0;
            query_core_Count++;
            found[query] = query;
            found_Count++;
            // vector<int> locations;
            // locations.push_back(query);
            // found.push_back(query);
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
                if (found[subject] == -1)
                // if (binary_search(found.begin(), found.end(), subject) == false)
                {
                    string subject_Hap = core_Haplotypes_All[subject];
                    // cout << query_Hap << "\t" << subject_Hap << endl;
                    if (query_Hap.compare(subject_Hap) == 0)
                    {
                        // cout << query_Hap << "\t" << subject_Hap << endl;
                        // cout << ext_Haplotypes_All[query] << "\t" << ext_Haplotypes_All[subject] << endl;
                        // locations.push_back(subject);
                        query_core_Count++;
                        found[subject] = subject;
                        found_Count++;
                        // found.push_back(subject);
                        //  hap_Count++;

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

                        // sort(found.begin(), found.end());
                    }
                }
            }

            Unique_haplotypes_sum_Partials.push_back(query_core_Count);

            int Ext_sum = 0;
            for (size_t unique_Haps = 0; unique_Haps < ext_Haplotypes_count.size(); unique_Haps++)
            {
                // cout << "Haps: " << unique_Haps + 1 << endl;
                // cout << ext_Haplotypes_count[unique_Haps].second << endl;
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

        if (found_Count == core_Haplotypes_All.size())
        {
            break;
        }
        // else
        // {
        //     sort(found.begin(), found.end());
        // }
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