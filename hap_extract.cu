#include "functions.cuh"
#include "hap_extract.cuh"

hap_extract::hap_extract(string gene_List, string input_Folder, string output_Path, int cuda_ID, string intermediate_Path, int ploidy, string reference_File, string pop_Out)
{
    /**
     * * Constructor Function
     * Assigns passed variables to the classes' private variable.
     **/

    cout << "Initiating CUDA powered Haplotype extractor" << endl
         << endl;

    set_Values(gene_List, input_Folder, output_Path, cuda_ID, intermediate_Path, ploidy);

    this->reference_File = reference_File;

    /**
     * Configures the population sequence print.
     * Conducts an uppercase conversion to prevent user error.
     **/
    transform(pop_Out.begin(), pop_Out.end(), pop_Out.begin(), ::toupper);
    if (pop_Out != "NO")
    {
        this->pop_Out = "YES";
    }
}

void hap_extract::set_Values(string gene_List, string input_Folder, string output_Path, int cuda_ID, string intermediate_Path, int ploidy)
{
    /**
     * This function is used in conjunction with the constructor to set the common private variables.
     * Here the first call to the selected CUDA device occurs.
     **/

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
    /**
     * Execution function.
     **/

    /**
     * Call the "functions" class. Bespoke functions commonly used by CATE.
     **/
    functions function = functions();

    /**
     * Read the reference file and load it into the GPU memory for processing.
     **/
    fstream reference;
    reference.open(this->reference_File, ios::in);

    /**
     * @param full_Reference used to capture the reference sequence into the RAM.
     **/
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

    /**
     * Converts the reference sequence to uppercase to provide uniformity.
     **/
    transform(full_Reference.begin(), full_Reference.end(), full_Reference.begin(), ::toupper);
    // this->reference_size = full_Reference.size();

    /**
     * @param reference_full is used to convert the string into a char pointer that can then be transferred into the GPU memory.
     **/
    char *reference_full;
    reference_full = (char *)malloc((full_Reference.size() + 1) * sizeof(char));
    cudaMallocManaged(&cuda_reference, (full_Reference.size() + 1) * sizeof(char));
    strcpy(reference_full, full_Reference.c_str());
    cudaMemcpy(cuda_reference, reference_full, (full_Reference.size() + 1) * sizeof(char), cudaMemcpyHostToDevice);

    /**
     * RAM is cleared of the reference sequence to prevent redundancy.
     **/
    free(reference_full);

    cout << "Reference file loaded" << endl
         << endl;

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
        vector<pair<string, string>> folder_Index = function.index_Folder(country);
        cout << "Completed indexing folder\t: " << country << endl;

        cout << endl;

        /**
         * The first VCF file is read to obtain information of the sample size.
         * @param samples captures the sample size of the population under study.
         **/
        int samples = function.getN_Split(folder_Index[0].second);
        cout << "Number of samples in " << country.substr(country.find_last_of("/") + 1, country.length()) << " population\t: " << samples << endl;

        this->N = samples * ploidy;
        cout << "Number of sequences in " << country.substr(country.find_last_of("/") + 1, country.length()) << " population [ " << samples << " x " << ploidy << " ] (N)\t: " << N << endl;

        cout << endl;

        /**
         * Initiate the reading of the gene file.
         **/
        fstream gene_File;
        gene_File.open(gene_List, ios::in);
        cout << "Processing gene list:" << endl;

        /**
         * Output file is created for the population in the output folder.
         * @param output_File stores the output file's location.
         * This file contains the summary information of the haplotypes. It is a tab deliminated text file.
         **/
        string output_File = output_Path + "/" +
                             country.substr(country.find_last_of("/") + 1, country.length()) + "_" +
                             filesystem::path(gene_List).stem().string() +
                             ".hsum";

        /**
         * Log file created in the intermediate folder for the population.
         * @param intermediate_File stores the log file's location.
         * ! This helps with the resume function. Automatically resumes from the last completely processed gene in the event of a program crash.
         **/
        string intermediate_File = intermediate_Path + "/" +
                                   country.substr(country.find_last_of("/") + 1, country.length()) + "_" +
                                   filesystem::path(gene_List).stem().string() +
                                   ".log_hap";
        cout << endl;
        cout << "Writing summary to file\t: " << output_File << endl;

        /**
         * FASTA folder created to store the resultant FASTA sequences.
         * Both the unique haplotype sequences and when required the population sequences will be printed here.
         * @param FASTA_folder stores the FASTA folder location.
         **/
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
                function.createFile(output_File, "Gene_name\tCoordinates\tHaplotype_number\tMutated_positions\tNumber_of_sequences");
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
                 * The SNPs (Segregating sites) that fall within the query region are collected from the VCF's.
                 * @param collect_Segregrating_sites vector stores the collected SNPs.
                 * @param pos_INDEX paired vector stores the SNP position and its position in the collect_Segregrating_sites vector,
                 * this helps to sort the SNPs by position (for faster binary searches) as well as keep track of their overall positions later on.
                 * @param count_Segs counts the number of collected SNPs.
                 **/
                vector<string> collect_Segregrating_sites;
                vector<pair<int, int>> pos_INDEX;
                int count_Segs = 0;

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

                /**
                 * Once the required files are found they are read sequentially to get the required SNP data for processing.
                 **/
                cout << "System is collecting segregating site(s)" << endl;
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
                                // cout << pos << endl;

                                /**
                                 * Required information on target SNP data is collected.
                                 **/
                                collect_Segregrating_sites.push_back(line);
                                pos_INDEX.push_back(make_pair(pos, count_Segs));
                                count_Segs++;
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

                // GET Haps

                /**
                 * Once the required data is collected CATE begins the process of Haplotype reconstruction.
                 * First it is ensured that for the query region there was SNP information.
                 **/
                if (collect_Segregrating_sites.size() != 0)
                {
                    /**
                     * @param write_Lines vector stores the output line that will be written to *.hsum file. Each Haplotype discovered will be written into a separate line.
                     * @param write_Sequences stores the sequences of unique haplotypes.
                     **/
                    vector<string> write_Lines, write_Sequences;

                    /**
                     * Initiate haplotype reconstruction.
                     **/
                    hap_extraction(write_Lines, write_Sequences, collect_Segregrating_sites, pos_INDEX, gene_Name, coordinates[0], start_Co, end_Co);

                    string FASTA_File = FASTA_folder + "/" + gene_Name + ".fasta";
                    fstream FASTA_out;
                    FASTA_out.open(FASTA_File, ios::out);

                    if (pop_Out != "NO")
                    {
                        /**
                         * If the user has selected population reconstruction, then,
                         * the entire population of FASTA sequences are reconstructed.
                         * @param FASTA_pop_File stored the population's FASTA sequence file's location information.
                         **/

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

                            /**
                             * @param number_of_Sequences is used to capture the number of copies of that haplotype that are present.
                             **/
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
    /**
     * @param tid is used to get the unique thread ID. In this instance thread ID is used to keep track of the Seg/SNP sites.
     **/
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < total_Segs)
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

        /**
         * @param REF is used capture the reference ALLELE.
         **/

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
         * Converts ASCII simple letters to CAPITAL, to create uniformity.
         **/

        if (sites[i] >= 97)
        {
            REF = sites[i] - 32;
        }
        else
        {
            REF = sites[i];
        }

        /**
         * @param ALT is used capture the alternate ALLELE.
         **/

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

        /**
         * Assigns the captured allelic information to the arrays for storage and post processing.
         **/

        REF_all[tid] = REF;
        ALT_all[tid] = ALT;

        /**
         * Beyond column 9 we grab sample information.
         **/

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
                /**
                 * We construct the haplotype while string along the SNP.
                 **/
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
    /**
     * @param tid is used to get the unique thread ID. In this instance thread ID is used to keep track of the Seg/SNP sites.
     **/
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < sequence_Size)
    {
        /**
         * EACH thread in this GPU call is a position on the FASTA sequence.
         * This enables the entire sequence to be reconstructed at once inside the GPU.
         **/

        /**
         * @param pos is use to capture the base pair position in the reference genome.
         **/

        int pos = tid + start;

        // check is pos exists in collection
        /**
         * We have to first check if the target region is present in the collected VCF information.
         * If it is not present we assume that no polymorphism has occurred at this site and simply take whatever allele is present
         * at that position in the reference genome.
         **/

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

        /**
         * If the allele is found in the search space we take the value from the VCF's SNP information.
         * 0 = REFERENCE ALLELE
         * 1 = ALTERNATE ALLELE
         * Else we take the value from the reference sequence.
         **/
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

        /**
         * Once the allele for that position is found it is assigned to that position of the FASTA sequence.
         **/
        sequence_Full[tid] = allele;

        tid += blockDim.x * gridDim.x;
    }
}

void hap_extract::hap_extraction(vector<string> &write_Lines, vector<string> &write_Sequences, vector<string> &total_Segregrating_sites, vector<pair<int, int>> &pos_INDEX, string gene_Name, string chr, int start_Pos, int end_Pos)
{
    /**
     * Administrative function responsible for haplotype reconstruction.
     * 1. Conversion of SNP strings into char pointers for GPU accessability.
     * 2. Call GPU for Haplotype reconstruction.
     * 3. Detect unique Haplotypes.
     * 4. Call GPU Sequence reconstruction.
     **/

    cout << "\nSystem is conducting Haplotype(s) construction" << endl;

    /**
     * @param num_segregrating_Sites is used to track the number of SNPs collected for the query region.
     * This track is vital for navigating through the data in the GPU. For the data is stored in the form of a 1D array.
     **/
    int num_segregrating_Sites = total_Segregrating_sites.size();

    /**
     * @param Seg_sites is used to stitch the SNP data end to end, before converting it to a char array.
     **/
    string Seg_sites = "";
    /**
     * @param site_Index is used to keep track of the start and ends of each SNP's data.
     **/
    int *site_Index;
    site_Index = (int *)malloc((num_segregrating_Sites + 1) * sizeof(int));
    site_Index[0] = 0;

    /**
     * @param pos used to store the collected SNP's position information.
     **/
    int *pos;
    //*cuda_pos;
    pos = (int *)malloc(num_segregrating_Sites * sizeof(int));

    vector<pair<int, int>> pos_INDEX_temp = pos_INDEX;
    /**
     * SNPs are sorted by position.
     **/
    sort(pos_INDEX_temp.begin(), pos_INDEX_temp.end());

    /**
     * @param pos_Allele stores the position of each SNP.
     * @param cuda_pos_Allele is used by the GPU. Is a COPY of pos_Allele.
     *
     * @param index_Allele stores the index of each SNP (location in the array).
     * @param cuda_index_Allele is used by the GPU. Is a COPY of index_Allele.
     **/
    int *pos_Allele, *cuda_pos_Allele, *index_Allele, *cuda_index_Allele;
    pos_Allele = (int *)malloc(num_segregrating_Sites * sizeof(int));
    index_Allele = (int *)malloc(num_segregrating_Sites * sizeof(int));

    /**
     * Conversion of vector SNP information into a 1D array by concatenating the vector data into a single string.
     **/
    for (size_t i = 0; i < num_segregrating_Sites; i++)
    {
        Seg_sites.append(total_Segregrating_sites[i]);
        site_Index[i + 1] = site_Index[i] + total_Segregrating_sites[i].size();

        pos[i] = pos_INDEX[i].first;

        pos_Allele[i] = pos_INDEX_temp[i].first;
        index_Allele[i] = pos_INDEX_temp[i].second;
    }

    /**
     * Final assignment of concatented string into a 1D char array.
     **/
    char *full_Char;
    full_Char = (char *)malloc((Seg_sites.size() + 1) * sizeof(char));
    strcpy(full_Char, Seg_sites.c_str());

    /**
     * RAM is released to prevent redundancy.
     **/
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
     * Transfer of data to the GPU.
     **/
    cudaMemcpy(cuda_full_Char, full_Char, (Seg_sites.size() + 1) * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(cuda_site_Index, site_Index, (num_segregrating_Sites + 1) * sizeof(int), cudaMemcpyHostToDevice);

    /**
     * @param cuda_REF_char used to capture the Reference allele in the SNP position.
     * @param cuda_ALT_char used to capture the Alternate allele in the SNP position.
     **/
    char *cuda_REF_char, *cuda_ALT_char;
    //*REF_char, *ALT_char,
    cudaMallocManaged(&cuda_REF_char, (num_segregrating_Sites + 1) * sizeof(char));
    cudaMallocManaged(&cuda_ALT_char, (num_segregrating_Sites + 1) * sizeof(char));

    /**
     * @param cuda_Hap_array stores the forged Haplotypes for the region under study.
     * @param Hap_array is used by the CPU. Is a COPY of cuda_Hap_array.
     **/
    char *Hap_array, *cuda_Hap_array;
    cudaMallocManaged(&cuda_Hap_array, ((this->N * num_segregrating_Sites) + 1) * sizeof(char));

    // ORGANIZE into array and collect Alleles
    /**
     * CALL THE GPU.
     * * GPU WILL CONDUCT HAPLOTYPE RECONSTRUCTION AND COLLECT THE ALLELIC INFORMATION REQUIRED FOR POST PROCESSING.
     **/
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

    /**
     * Conversion of haplotype information to strings to enable CPU based comparisons.
     **/

    /**
     * @param haplotypes captures the haplotype array in string format. Allows substr function and easier submission into a vector.
     * @param Haplotypes_All vector collects all individual haplotypes.
     **/
    string haplotypes(Hap_array);
    vector<string> Haplotypes_All;

    for (int i = 0; i < (num_segregrating_Sites * this->N); i = i + num_segregrating_Sites)
    {
        // cout << ext_Haplotypes.substr(i, num_segregrating_Sites) << endl;
        Haplotypes_All.push_back(haplotypes.substr(i, num_segregrating_Sites));
    }

    /**
     * @param found vector is used to track all found haplotype. The position of the haplotype in the Haplotypes_All vector is recorded.
     **/
    vector<int> found;

    // cudaMallocManaged(&cuda_pos, (num_segregrating_Sites * sizeof(int)));
    // cudaMemcpy(cuda_pos, pos, num_segregrating_Sites * sizeof(int), cudaMemcpyHostToDevice);

    cudaMallocManaged(&cuda_pos_Allele, (num_segregrating_Sites * sizeof(int)));
    cudaMallocManaged(&cuda_index_Allele, (num_segregrating_Sites * sizeof(int)));

    cudaMemcpy(cuda_pos_Allele, pos_Allele, num_segregrating_Sites * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(cuda_index_Allele, index_Allele, num_segregrating_Sites * sizeof(int), cudaMemcpyHostToDevice);

    /**
     * @param sequence_Size is used to determine the sequence length. Required for the sequence reconstruction step.
     * @param HAP_ID is used to keep track of the number of unique haplotypes per query region. Ensure they each have an unique ID.
     **/
    int sequence_Size = end_Pos - start_Pos + 1;
    int HAP_ID = 0;

    cout << "STEP 2 OF 2: Detecting unique haplotypes and synthesizing their sequences" << endl;
    /**
     * We go through the collected Haplotypes to identify the unique ones.
     * If a unique haplotype is found that was not previously present its occurrence in the sample space is accounted for.
     * Found haplotypes are tracked by keeping track of their location in the found vector.
     * Since the vector stores numerical locations they can be sorted and fast search algorithms such as binary search can be employed.
     **/
    for (size_t query = 0; query < Haplotypes_All.size(); query++)
    {
        /**
         * If the query haplotype is NOT present in the found vector it is considered as a new unique haplotype.
         * This will trigger the haplotype processing algorithm.
         **/
        if (binary_search(found.begin(), found.end(), query) == false)
        {
            /**
             * HAP_ID is incremented by 1.
             **/
            HAP_ID = HAP_ID + 1;

            /**
             * @param query_core_Count is used to keep track of the number of occurrences of the haplotype in sample space.
             * It is incremented eah time an occurrence of the query haplotype is found.
             **/
            int query_core_Count = 0;
            query_core_Count++;

            /**
             * The newly discovered UNIQUE haplotype is recorded.
             **/
            found.push_back(query);

            /**
             * @param query_Hap is used to capture the query haplotype.
             **/
            string query_Hap = Haplotypes_All[query];
            // cout << query_Hap << endl;

            for (size_t subject = query + 1; subject < Haplotypes_All.size(); subject++)
            {
                /**
                 * Comparison of the query haplotype in the haplotype search space.
                 * We skip over previously accounted for haplotypes to prevent redundancy and improve speed.
                 **/
                if (binary_search(found.begin(), found.end(), subject) == false)
                {
                    /**
                     * @param subject_Hap is used to get the subject haplotype to be compared to the query haplotype (query_Hap).
                     * If they match this haplotype is recorded and accounted for.
                     **/
                    string subject_Hap = Haplotypes_All[subject];
                    if (query_Hap.compare(subject_Hap) == 0)
                    {
                        query_core_Count++;
                        /**
                         * The newly discovered haplotype is recorded.
                         **/
                        found.push_back(subject);
                        /**
                         * The modified vector is now sorted to enable binary search.
                         **/
                        sort(found.begin(), found.end());
                    }
                }
            }

            // Process this Haplotype
            /**
             * Once all occurrences of the unique query haplotype are recorded it will be subjected to processing
             * to reconstruct the complete sequence.
             **/

            /**
             * @param write is used to store the haplotype information to be written to the *.hsum output file.
             **/
            string write = gene_Name + "\t" +
                           chr + ":" + to_string(start_Pos) + ":" + to_string(end_Pos) + "\t" +
                           to_string(HAP_ID) + "\t";

            /**
             * @param mutation_positions is used to record the positions where mutations/polymorphisms have occurred.
             * * TIP: IF mutation_positions remain at "NA" then this means that the haplotype is equal to the reference.
             **/
            string mutation_positions = "NA";

            for (int stride = 0; stride < num_segregrating_Sites; stride++)
            {
                /**
                 * 1 is equal to a mutation. The alternate allele was present at this pos.
                 **/
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
            /**
             * Reconstruction of the FASTA sequence for the query haplotype.
             * @param sequence_ID is used to give an ID to the sequence. It is a combination of the region's gene name, coordinates and haplotype;s unique ID.
             **/
            string sequence_ID = ">" + gene_Name + "_" + chr + ":" + to_string(start_Pos) + ":" + to_string(end_Pos) + "_" + to_string(HAP_ID) + "\n";

            /**
             * @param sequence is used to store the reconstructed FASTA sequence.
             * @param cuda_sequence is used by the CPU. Is a COPY of sequence.
             **/
            char *sequence, *cuda_sequence;
            sequence = (char *)malloc((sequence_Size + 1) * sizeof(char));
            cudaMallocManaged(&cuda_sequence, (sequence_Size + 1) * sizeof(char));

            /**
             * @param haplotye is used to store the query haplotype VCF string sequence.
             * @param cuda_haplotype is used by the CPU. Is a COPY of haplotype.
             **/
            char *haplotye, *cuda_haplotype;
            haplotye = (char *)malloc((num_segregrating_Sites + 1) * sizeof(char));
            cudaMallocManaged(&cuda_haplotype, (num_segregrating_Sites + 1) * sizeof(char));

            strcpy(haplotye, query_Hap.c_str());
            cudaMemcpy(cuda_haplotype, haplotye, (num_segregrating_Sites + 1) * sizeof(char), cudaMemcpyHostToDevice);

            // cuda_sequence_Generation(int sequence_Size, int num_of_Segs, int start, char *ref, char *haplotype, int *pos_Allele, char *index_Allele, char *REF, char *ALT, char *sequence_Full)
            /**
             * CALL THE GPU.
             * * GPU WILL CONDUCT SEQUENCE RECONSTRUCTION.
             **/
            cuda_sequence_Generation<<<tot_Blocks, tot_ThreadsperBlock>>>(sequence_Size, num_segregrating_Sites, start_Pos, cuda_reference, cuda_haplotype, cuda_pos_Allele, cuda_index_Allele, cuda_REF_char, cuda_ALT_char, cuda_sequence);

            if (err != cudaSuccess)
            {
                printf("CUDA Error: %s\n", cudaGetErrorString(err));

                // Possibly: exit(-1) if program cannot continue....
            }
            cudaDeviceSynchronize();

            cudaMemcpy(sequence, cuda_sequence, (sequence_Size + 1) * sizeof(char), cudaMemcpyDeviceToHost);

            /**
             * @param sequence_string is used to convert the char sequence array to a string to be written to the FASTA output file.
             **/
            string sequence_string = sequence;
            sequence_string = sequence_string.substr(0, sequence_Size);

            string sequence_Full = sequence_ID + sequence_string;
            write_Sequences.push_back(sequence_Full);

            cudaFree(cuda_haplotype);
            cudaFree(cuda_sequence);
            free(haplotye);
            free(sequence);
        }

        /**
         * If ALL haplotypes under study are processed then the loop is broken.
         **/
        if (found.size() == Haplotypes_All.size())
        {
            break;
        }
        else
        {
            /**
             * The modified vector is now sorted to enable binary search.
             **/
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