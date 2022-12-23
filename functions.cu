#include "functions.cuh"

functions::functions()
{
}

__global__ void cuda_process_Seg_tajima(char *sites, int *index, int tot_Segregrating_sites, int *VALID_or_NOT, int *MA_count)
{
    /**
     * @param tid is used to get the unique thread ID. In this instance thread ID is used to keep track of the Seg/SNP sites.
     **/
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < tot_Segregrating_sites)
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

        while (column < 9)
        {
            if (sites[i] == '\t')
            {
                column++;
            }
            i++;
        }

        // printf("Column 1: %c\n", sites[i]);

        /**
         * @param ALT_count is used capture number of instances the ALTERNATE allele is present.
         * @param REF_count is used capture number of instances the REFERENCE allele is present.
         **/
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

        /**
         * If either allele has a count of zero this means that the site is not a SEG site for that population.
         **/
        if (ALT_count == 0 || REF_count == 0)
        {
            VALID_or_NOT[tid] = 0;
            MA_count[tid] = 0;
        }
        else
        {
            VALID_or_NOT[tid] = 1;

            if (ALT_count < REF_count)
            {
                MA_count[tid] = ALT_count;
            }
            else
            {
                MA_count[tid] = REF_count;
            }
        }

        tid += blockDim.x * gridDim.x;
    }
}

void functions::process_Seg_sites_tajima(vector<string> &total_Segregrating_sites, int N, int &real_segregrating_Sites, float &tot_pairwise_Differences, int tot_Blocks, int tot_ThreadsperBlock)
{
    /**
     * Administrative function responsible for collection of data required for Tajima's D.
     * 1. Conversion of SNP strings into char pointers for GPU accessability.
     * 2. Call GPU for extracting MAs (Minor allele) and MAF's (Minor Allele Frequencies).
     **/

    cout << "System is processing and filtering segregating site(s)" << endl;
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
    int site_Index[num_segregrating_Sites + 1];
    site_Index[0] = 0;

    /**
     * Conversion of vector SNP information into a 1D array by concatenating the vector data into a single string.
     **/
    for (size_t i = 0; i < num_segregrating_Sites; i++)
    {
        Seg_sites.append(total_Segregrating_sites[i]);
        site_Index[i + 1] = site_Index[i] + total_Segregrating_sites[i].size();
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

    // cuda_process_Seg(char *sites, int *index, int tot_Segregrating_sites, int *VALID_or_NOT, int *MA_count)
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
     * Certain collected SNP sites though present in the VCF, specially if the index was created without a filter,
     * will have MAF of zero. Indicating that region is conserved for that super-pop and therfore not a seg site.
     * These have to be filtered out. This is done using the variables below:
     * @param cuda_VALID_or_NOT is used to determine if a site is seg site which is VALID or NOT.
     * @param VALID_or_NOT is used by the CPU. Is a COPY of cuda_VALID_or_NOT.
     **/
    int *cuda_VALID_or_NOT, *VALID_or_NOT;
    cudaMallocManaged(&cuda_VALID_or_NOT, num_segregrating_Sites * sizeof(int));
    VALID_or_NOT = (int *)malloc(num_segregrating_Sites * sizeof(int));

    /**
     * @param cuda_MA_Count is used to record the MA's count.
     * @param MA_Count is used by the CPU. Is a COPY of cuda_MA_Count.
     **/
    int *cuda_MA_Count, *MA_Count;
    cudaMallocManaged(&cuda_MA_Count, num_segregrating_Sites * sizeof(int));
    MA_Count = (int *)malloc(num_segregrating_Sites * sizeof(int));

    /**
     * CALL THE GPU.
     * * GPU WILL PROCESS THE COLLECTED SEG SITES
     **/
    cuda_process_Seg_tajima<<<tot_Blocks, tot_ThreadsperBlock>>>(cuda_full_Char, cuda_site_Index, num_segregrating_Sites, cuda_VALID_or_NOT, cuda_MA_Count);
    cudaDeviceSynchronize();

    cudaMemcpy(VALID_or_NOT, cuda_VALID_or_NOT, num_segregrating_Sites * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(MA_Count, cuda_MA_Count, num_segregrating_Sites * sizeof(int), cudaMemcpyDeviceToHost);

    free(full_Char);
    cudaFree(cuda_site_Index);
    cudaFree(cuda_MA_Count);
    cudaFree(cuda_VALID_or_NOT);

    real_segregrating_Sites = 0;
    tot_pairwise_Differences = 0;

    /**
     * Processing of each site is done to validate them and extract the data required for the statistical test.
     * Tajima's D required the number of actual segregating sites and the total number of pairwise differences.
     **/
    for (size_t i = 0; i < num_segregrating_Sites; i++)
    {
        if (VALID_or_NOT[i] == 1)
        {
            real_segregrating_Sites = real_segregrating_Sites + 1;
            float MAF = (float)MA_Count[i] / (float)N;
            tot_pairwise_Differences = tot_pairwise_Differences + (MAF * (1 - MAF) * pow(N, 2));
        }
    }

    free(MA_Count);
    free(VALID_or_NOT);

    cout << "System has completed processing the segregrating site(s)" << endl;
}

__global__ void cuda_process_Seg_fu_li(char *sites, int *index, int tot_Segregrating_sites, int *VALID_or_NOT, int *MA_count, int *ne, int *ns)
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

        if (ALT_count == 0 || REF_count == 0)
        {
            VALID_or_NOT[tid] = 0;
            MA_count[tid] = 0;
            ne[tid] = 0;
            ns[tid] = 0;
        }
        else
        {
            VALID_or_NOT[tid] = 1;

            char MA = 'N';

            if (ALT_count < REF_count)
            {
                MA_count[tid] = ALT_count;
                MA = ALT;
            }
            else
            {
                MA = REF;
                MA_count[tid] = REF_count;
            }

            if (MA_count[tid] == 1)
            {
                ns[tid] = 1;

                if (MA == AA)
                {
                    ne[tid] = 0;
                }
                else
                {
                    ne[tid] = 1;
                }
            }
            else
            {
                ne[tid] = 0;
                ns[tid] = 0;
            }
        }

        tid += blockDim.x * gridDim.x;
    }
}

void functions::process_Seg_sites_fu_li(vector<string> &total_Segregrating_sites, float N, int &real_segregrating_Sites, float &tot_pairwise_Differences, int &singletons_ne, int &singletons_ns, int tot_Blocks, int tot_ThreadsperBlock)
{
    cout << "System is processing and filtering segregrating site(s)" << endl;
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
    cudaMemcpy(cuda_full_Char, full_Char, (Seg_sites.size() + 1) * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(cuda_site_Index, site_Index, (num_segregrating_Sites + 1) * sizeof(int), cudaMemcpyHostToDevice);

    int *cuda_VALID_or_NOT, *VALID_or_NOT;
    cudaMallocManaged(&cuda_VALID_or_NOT, num_segregrating_Sites * sizeof(int));
    VALID_or_NOT = (int *)malloc(num_segregrating_Sites * sizeof(int));

    int *cuda_MA_Count, *MA_Count;
    cudaMallocManaged(&cuda_MA_Count, num_segregrating_Sites * sizeof(int));
    MA_Count = (int *)malloc(num_segregrating_Sites * sizeof(int));

    int *ne_CUDA, *ne, *ns_CUDA, *ns;
    cudaMallocManaged(&ne_CUDA, num_segregrating_Sites * sizeof(int));
    cudaMallocManaged(&ns_CUDA, num_segregrating_Sites * sizeof(int));
    ne = (int *)malloc(num_segregrating_Sites * sizeof(int));
    ns = (int *)malloc(num_segregrating_Sites * sizeof(int));

    // cuda_process_Seg_fu_li(char *sites, int *index, int tot_Segregrating_sites, int *VALID_or_NOT, int *MA_count, int *ne, int *ns)
    cuda_process_Seg_fu_li<<<tot_Blocks, tot_ThreadsperBlock>>>(cuda_full_Char, cuda_site_Index, num_segregrating_Sites, cuda_VALID_or_NOT, cuda_MA_Count, ne_CUDA, ns_CUDA);
    cudaDeviceSynchronize();

    cudaMemcpy(VALID_or_NOT, cuda_VALID_or_NOT, num_segregrating_Sites * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(MA_Count, cuda_MA_Count, num_segregrating_Sites * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(ne, ne_CUDA, num_segregrating_Sites * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(ns, ns_CUDA, num_segregrating_Sites * sizeof(int), cudaMemcpyDeviceToHost);

    free(full_Char);
    cudaFree(cuda_full_Char);
    cudaFree(cuda_site_Index);
    cudaFree(cuda_MA_Count);
    cudaFree(cuda_VALID_or_NOT);
    cudaFree(ns_CUDA);
    cudaFree(ne_CUDA);

    real_segregrating_Sites = 0;
    tot_pairwise_Differences = 0;
    singletons_ne = 0;
    singletons_ns = 0;

    for (size_t i = 0; i < num_segregrating_Sites; i++)
    {
        if (VALID_or_NOT[i] == 1)
        {
            real_segregrating_Sites = real_segregrating_Sites + 1;
            float MAF = (float)MA_Count[i] / N;
            tot_pairwise_Differences = tot_pairwise_Differences + (MAF * (1 - MAF) * pow(N, 2));
            singletons_ne = singletons_ne + ne[i];
            singletons_ns = singletons_ns + ns[i];
        }
    }

    free(MA_Count);
    free(VALID_or_NOT);
    free(ns);
    free(ne);

    cout << "System has completed processing the segregrating site(s)" << endl;
}

vector<string> functions::get_Countries(string &input_Folder)
{
    vector<string> folders;
    for (auto &check : std::filesystem::recursive_directory_iterator(input_Folder))
    {
        if (check.is_directory())
        {
            folders.push_back(check.path().string());
        }
    }
    return folders;
}

vector<pair<string, string>> functions::index_Folder(string &country)
{
    cout << "Initiating indexing folder\t: " << country << endl;
    string country_Only = country.substr(country.find_last_of("/") + 1, country.length());

    vector<string> index_pass_1;
    vector<pair<string, string>> file_coordinate;

    for (const auto &entry : filesystem::directory_iterator(country))
    {

        string coordinates = entry.path().string();
        int trim_start = coordinates.find(country_Only + "_") + country_Only.length() + 1;
        int trim_end = coordinates.find_last_of(".") - trim_start;
        string trim = coordinates.substr(trim_start, trim_end);
        index_pass_1.push_back(trim);
        file_coordinate.push_back(make_pair(trim, coordinates));
    }

    vector<pair<int, int>> start_stop;
    vector<int> starts;

    for (string file : index_pass_1)
    {
        vector<string> file_start_end;
        split(file_start_end, file, '_');
        starts.push_back(stoi(file_start_end[0]));
        start_stop.push_back(make_pair(stoi(file_start_end[0]), stoi(file_start_end[1])));
    }

    sort(starts.begin(), starts.end());

    vector<pair<string, string>> sorted_Index;

    for (int nums : starts)
    {
        for (auto index_check : start_stop)
        {
            if (index_check.first == nums)
            {
                string sort_Line = to_string(index_check.first) + "_" + to_string(index_check.second);
                for (auto coordinates : file_coordinate)
                {
                    if (coordinates.first == sort_Line)
                    {
                        sorted_Index.push_back(make_pair(sort_Line, coordinates.second));
                        break;
                    }
                }
                break;
            }
        }
    }

    return sorted_Index;
}

void functions::split(vector<string> &line_Data, string line, char delim)
{
    line_Data.clear();

    // string coordinates_String = line.substr(gene_Combo.find(delim) + 1);
    while (line.find(delim) != string::npos)
    {
        line_Data.push_back(line.substr(0, line.find(delim)));
        line = line.substr(line.find(delim) + 1);
    }

    if (!line.empty())
    {
        if (line[line.length() - 1] == '\r' || line[line.length() - 1] == '\n' || line[line.length() - 1] == '\t')
        {
            line.erase(line.length() - 1);
        }
    }

    // line.erase(remove(line.begin(), line.end(), '\r'), line.cend());

    line_Data.push_back(line);
}

void functions::split_space(vector<string> &line_Data, string line, string delim)
{

    vector<string>().swap(line_Data);
    char *convert;
    string capture(line);
    convert = &capture[0];

    char deliminator[delim.length() + 1];
    strcpy(deliminator, delim.c_str());

    char *split_data;
    split_data = strtok(convert, deliminator);

    while (split_data != NULL)
    {
        string char2string;
        char2string.append(split_data);
        line_Data.push_back(char2string);
        split_data = strtok(NULL, deliminator);
    }
}

int functions::getN_Split(string file)
{
    fstream file_nCount;
    file_nCount.open(file);

    string header;
    getline(file_nCount, header);
    file_nCount.close();

    vector<string> header_Columns;
    split(header_Columns, header, '\t');

    int N = header_Columns.size() - 9;
    // cout << N << endl;
    return N;
}

__global__ void cuda_array_ADD(const float *a, float *out, int arraySize)
{
    int idx = threadIdx.x;
    float sum = 0;
    for (int i = idx; i < arraySize; i += 1024)
        sum += a[i];
    __shared__ float r[1024];
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

float functions::add(int N, float *array)
{

    float *x_Cuda, *y_Cuda;
    float *y_partial = (float *)malloc(sizeof(float));

    cudaMalloc((void **)&x_Cuda, N * sizeof(float));
    cudaMemcpy(x_Cuda, array, N * sizeof(float), cudaMemcpyHostToDevice);

    cudaMalloc((void **)&y_Cuda, sizeof(float));

    cuda_array_ADD<<<1, 1024>>>(x_Cuda, y_Cuda, N);

    cudaDeviceSynchronize();
    cudaMemcpy(y_partial, y_Cuda, sizeof(float), cudaMemcpyDeviceToHost);

    cudaFree(x_Cuda);
    cudaFree(y_Cuda);

    return y_partial[0];

    free(y_partial);
}

long int functions::fact_half(int count)
{
    long int tot = 1;
    for (int i = count; i > count - 2; i--)
    {
        // cout << tot;
        tot = tot * i;
    }
    return tot;
}

long int functions::combos_N(int count)
{
    long int combinations;

    combinations = fact_half(count) / 2;

    return combinations;
}

void functions::createFile(string path)
{
    fstream file;
    file.open(path, ios::out);
    file.close();
}

void functions::createFile(string path, string text)
{
    fstream file;
    file.open(path, ios::out);
    file << text;
    file << "\n";
    file.close();
}

string functions::roundoff(float value, unsigned char prec)
{
    float pow_10 = pow(10.0f, (float)prec);
    return to_string(round(value * pow_10) / pow_10).substr(0, prec + 2);
}

void functions::backward_Search(int pos, vector<pair<string, string>> folder_Index, int start_Co, int end_Co, promise<vector<int>> &backward_Found)
{
    vector<int> backward_get;
    pos = pos - 1;
    vector<string> line_Data_get;
    while (pos >= 0)
    {
        split(line_Data_get, folder_Index[pos].first, '_');
        int low_Value_atpos = stoi(line_Data_get[0]);
        int high_Value_atpos = stoi(line_Data_get[1]);

        if (start_Co > high_Value_atpos)
        {
            break;
        }

        if ((start_Co >= low_Value_atpos && start_Co <= high_Value_atpos) || (start_Co <= low_Value_atpos && end_Co <= high_Value_atpos) || (start_Co <= low_Value_atpos && end_Co >= high_Value_atpos))
        {
            backward_get.push_back(pos);
        }

        pos = pos - 1;
    }
    backward_Found.set_value(backward_get);
    // return backward_get;
}

void functions::forward_Search(int pos, vector<pair<string, string>> folder_Index, int start_Co, int end_Co, promise<vector<int>> &forward_Found)
{
    vector<int> forward_get;
    pos = pos + 1;
    vector<string> line_Data_get;
    while (pos < folder_Index.size())
    {
        split(line_Data_get, folder_Index[pos].first, '_');
        int low_Value_atpos = stoi(line_Data_get[0]);
        int high_Value_atpos = stoi(line_Data_get[1]);

        if (end_Co < low_Value_atpos)
        {
            break;
        }

        if ((start_Co >= low_Value_atpos && start_Co <= high_Value_atpos) || (start_Co <= low_Value_atpos && end_Co <= high_Value_atpos) || (start_Co <= low_Value_atpos && end_Co >= high_Value_atpos))
        {
            forward_get.push_back(pos);
        }

        pos = pos + 1;
    }
    forward_Found.set_value(forward_get);
    // return forward_get;
}

int functions::single_segment_retrieval(int SNP_position, vector<pair<string, string>> &folder_Index)
{
    int segment_Position = -1;

    vector<string> line_Data;
    split(line_Data, folder_Index[0].first, '_');
    int low_Value = stoi(line_Data[0]);
    split(line_Data, folder_Index[folder_Index.size() - 1].first, '_');
    int high_Value = stoi(line_Data[1]);

    int start = 0;
    int end = folder_Index.size() - 1;

    while (start <= end && SNP_position <= high_Value)
    {
        vector<string> line_Data_get;

        int pos = start + ((((double)(end - start) / (high_Value - low_Value)) * (SNP_position - low_Value)));

        split(line_Data_get, folder_Index[pos].first, '_');
        int low_Value_atpos = stoi(line_Data_get[0]);
        int high_Value_atpos = stoi(line_Data_get[1]);

        if ((SNP_position >= low_Value_atpos) && (SNP_position <= high_Value_atpos))
        {
            segment_Position = pos;
            break;
        }
        else if (SNP_position > low_Value_atpos)
        {
            // cout << pos << "\t";
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

        split(line_Data_get, folder_Index[start].first, '_');
        low_Value = stoi(line_Data_get[0]);

        split(line_Data_get, folder_Index[end].first, '_');
        high_Value = stoi(line_Data_get[1]);
    }

    return segment_Position;
}

vector<string> functions::compound_interpolationSearch_ordered(vector<pair<string, string>> &folder_Index, int &start_Co, int &end_Co)
{
    vector<string> file_List;

    vector<string> line_Data;
    split(line_Data, folder_Index[0].first, '_');
    int low_Value = stoi(line_Data[0]);
    split(line_Data, folder_Index[folder_Index.size() - 1].first, '_');
    int high_Value = stoi(line_Data[1]);

    int start = 0;
    int end = folder_Index.size() - 1;

    while (start <= end && start_Co <= high_Value)
    {
        vector<string> line_Data_get;
        int pos = start + ((((double)(end - start) / (high_Value - low_Value)) * (start_Co - low_Value)));

        split(line_Data_get, folder_Index[pos].first, '_');
        int low_Value_atpos = stoi(line_Data_get[0]);
        int high_Value_atpos = stoi(line_Data_get[1]);

        // cout << low_Value_atpos << "\t" << high_Value_atpos << endl;

        if ((start_Co >= low_Value_atpos && start_Co <= high_Value_atpos) || (start_Co <= low_Value_atpos && end_Co <= high_Value_atpos) || (start_Co <= low_Value_atpos && end_Co >= high_Value_atpos))
        {

            promise<vector<int>> backward;
            promise<vector<int>> forward;

            future<vector<int>> fut_Backward = backward.get_future();
            future<vector<int>> fut_Forward = forward.get_future();

            thread backward_thread{&functions::backward_Search, this, pos, folder_Index, start_Co, end_Co, ref(backward)};
            thread forward_thread{&functions::forward_Search, this, pos, folder_Index, start_Co, end_Co, ref(forward)};

            vector<int> backward_get = fut_Backward.get();
            vector<int> forward_get = fut_Forward.get();

            backward_thread.join();
            forward_thread.join();

            // auto backward_thread = async(&functions::backward_Search, this, pos, folder_Index, start_Co, end_Co);
            // auto forward_thread = async(&functions::forward_Search, this, pos, folder_Index, start_Co, end_Co);

            // backward_get = backward_thread.get();
            // forward_get = forward_thread.get();

            for (size_t i = backward_get.size() - 1; i <= 0; i--)
            {
                file_List.push_back(folder_Index[backward_get[i]].second);
            }

            // for (auto positions : backward_get)
            // {
            //     file_List.push_back(folder_Index[positions].second);
            // }

            file_List.push_back(folder_Index[pos].second);

            for (auto positions : forward_get)
            {
                file_List.push_back(folder_Index[positions].second);
            }

            break;
        }
        else if (start_Co > low_Value_atpos)
        {
            // cout << pos << "\t";
            int new_pos = pos;

            do
            {
                new_pos = new_pos + 1;
            } while (new_pos <= start);

            start = new_pos;

            // cout << start << endl;
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

        split(line_Data_get, folder_Index[start].first, '_');
        low_Value = stoi(line_Data_get[0]);

        split(line_Data_get, folder_Index[end].first, '_');
        high_Value = stoi(line_Data_get[1]);
    }

    return file_List;
}

vector<string> functions::compound_interpolationSearch(vector<pair<string, string>> &folder_Index, int &start_Co, int &end_Co)
{
    vector<string> file_List;

    vector<string> line_Data;
    split(line_Data, folder_Index[0].first, '_');
    int low_Value = stoi(line_Data[0]);
    split(line_Data, folder_Index[folder_Index.size() - 1].first, '_');
    int high_Value = stoi(line_Data[1]);

    int start = 0;
    int end = folder_Index.size() - 1;

    while (start <= end && start_Co <= high_Value)
    {
        vector<string> line_Data_get;
        int pos = start + ((((double)(end - start) / (high_Value - low_Value)) * (start_Co - low_Value)));

        split(line_Data_get, folder_Index[pos].first, '_');
        int low_Value_atpos = stoi(line_Data_get[0]);
        int high_Value_atpos = stoi(line_Data_get[1]);

        // cout << low_Value_atpos << "\t" << high_Value_atpos << endl;

        if ((start_Co >= low_Value_atpos && start_Co <= high_Value_atpos) || (start_Co <= low_Value_atpos && end_Co <= high_Value_atpos) || (start_Co <= low_Value_atpos && end_Co >= high_Value_atpos))
        {

            promise<vector<int>> backward;
            promise<vector<int>> forward;

            future<vector<int>> fut_Backward = backward.get_future();
            future<vector<int>> fut_Forward = forward.get_future();

            thread backward_thread{&functions::backward_Search, this, pos, folder_Index, start_Co, end_Co, ref(backward)};
            thread forward_thread{&functions::forward_Search, this, pos, folder_Index, start_Co, end_Co, ref(forward)};

            vector<int> backward_get = fut_Backward.get();
            vector<int> forward_get = fut_Forward.get();

            backward_thread.join();
            forward_thread.join();

            // auto backward_thread = async(&functions::backward_Search, this, pos, folder_Index, start_Co, end_Co);
            // auto forward_thread = async(&functions::forward_Search, this, pos, folder_Index, start_Co, end_Co);

            // backward_get = backward_thread.get();
            // forward_get = forward_thread.get();

            for (auto positions : backward_get)
            {
                file_List.push_back(folder_Index[positions].second);
            }

            file_List.push_back(folder_Index[pos].second);

            for (auto positions : forward_get)
            {
                file_List.push_back(folder_Index[positions].second);
            }

            break;
        }
        else if (start_Co > low_Value_atpos)
        {
            // cout << pos << "\t";
            int new_pos = pos;

            do
            {
                new_pos = new_pos + 1;
            } while (new_pos <= start);

            start = new_pos;

            // cout << start << endl;
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

        split(line_Data_get, folder_Index[start].first, '_');
        low_Value = stoi(line_Data_get[0]);

        split(line_Data_get, folder_Index[end].first, '_');
        high_Value = stoi(line_Data_get[1]);
    }

    return file_List;
}

// void functions::split_getPos(vector<string> &line_Data, string line, string delim)
// {
//     vector<string>().swap(line_Data);
//     char *convert;
//     string capture(line);
//     convert = &capture[0];

//     char deliminator[delim.length() + 1];
//     strcpy(deliminator, delim.c_str());

//     char *split_data;
//     split_data = strtok(convert, deliminator);
//     int count = 0;

//     while (split_data != NULL)
//     {
//         string char2string;
//         char2string.append(split_data);
//         line_Data.push_back(char2string);
//         if (count == 7)
//         {
//             break;
//         }
//         split_data = strtok(NULL, deliminator);
//         count++;
//     }
// }

int functions::get_Valid(string line)
{
    int VALID = -1;
    int POS = -1;
    int count = 0;
    char delim = '\t';

    while (line.find(delim) != string::npos)
    {
        if (count == 1)
        {
            // line_Data.push_back(line.substr(0, line.find(delim)));
            POS = stoi(line.substr(0, line.find(delim)));
        }

        if (count == 9)
        {
            line = line.substr(line.find(delim) + 1);
            break;
        }

        line = line.substr(line.find(delim) + 1);
        count++;
    }

    int found_0 = 0;
    int found_1 = 0;
    // cout << POS << endl;

    for (int i = 0; i < line.size(); i++)
    {
        char val = line.at(i);

        if (val == '0')
        {
            if (found_0 == 0)
            {
                found_0 = 1;
            }
        }
        else if (val == '1')
        {
            if (found_1 == 0)
            {
                found_1 = 1;
            }
        }

        // cout << val << endl;
        // cout << val << "\t" << found_0 << "\t" << found_1 << endl;

        if ((found_0 == 1) && (found_1 == 1))
        {
            VALID = POS;
            break;
        }
    }

    return VALID;
}

void functions::split_getPos_ONLY(vector<string> &line_Data, string line, char delim)
{
    line_Data.clear();
    int count = 0;

    // string coordinates_String = line.substr(gene_Combo.find(delim) + 1);
    while (line.find(delim) != string::npos)
    {
        line_Data.push_back(line.substr(0, line.find(delim)));

        if (count == 1)
        {
            break;
        }

        line = line.substr(line.find(delim) + 1);
        count++;
    }

    // vector<string>().swap(line_Data);
    // char *convert;
    // string capture(line);
    // convert = &capture[0];

    // char deliminator[delim.length() + 1];
    // strcpy(deliminator, delim.c_str());

    // char *split_data;
    // split_data = strtok(convert, deliminator);
    // int count = 0;

    // while (split_data != NULL)
    // {
    //     string char2string;
    //     char2string.append(split_data);
    //     line_Data.push_back(char2string);
    //     if (count == 1)
    //     {
    //         break;
    //     }
    //     split_data = strtok(NULL, deliminator);
    //     count++;
    // }
}

__global__ void pairwise_Cuda_(int N, int *SNP, int *differences)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < N)
    {
        int tot = 0;
        for (int i = tid + 1; i < N; i++)
        {
            if (SNP[tid] != SNP[i])
            {
                tot = tot + 1;
            }
        }
        differences[tid] = tot;
        tid += blockDim.x * gridDim.x;
    }
}

int functions::calc_Pairwise(string &line, int N, int tot_Blocks, int tot_ThreadsperBlock)
{
    int pairwise_Differences = 0;

    int *line_temp = new int[N];
    split_Convert(line_temp, line, "\t");

    int *cuda_line_Data;
    cudaMallocManaged(&cuda_line_Data, N * sizeof(int));

    int *differences, *cuda_Differences;
    cudaMallocManaged(&cuda_Differences, N * sizeof(int));
    differences = (int *)malloc(N * sizeof(int));

    cudaMemcpy(cuda_line_Data, line_temp, (N * sizeof(int)), cudaMemcpyHostToDevice);

    pairwise_Cuda_<<<tot_Blocks, tot_ThreadsperBlock>>>(N, cuda_line_Data, cuda_Differences);
    cudaDeviceSynchronize();

    cudaMemcpy(differences, cuda_Differences, N * sizeof(int), cudaMemcpyDeviceToHost);

    cudaFree(cuda_line_Data);
    cudaFree(cuda_Differences);

    for (int i = 0; i < N; i++)
    {
        pairwise_Differences = pairwise_Differences + differences[i];
    }

    free(differences);
    free(line_temp);

    return pairwise_Differences;
}

void functions::split_Convert(int *line_temp, string line, string delim)
{
    char *convert;
    string capture(line);
    convert = &capture[0];

    char deliminator[delim.length() + 1];
    strcpy(deliminator, delim.c_str());

    char *split_data;
    split_data = strtok(convert, deliminator);
    int count = 0;

    while (split_data != NULL)
    {
        string char2string;
        char2string.append(split_data);
        if (count >= 9)
        {
            if (char2string.substr(0, 3) == "0|0")
            {
                line_temp[count - 9] = 0;
            }
            else if (char2string.substr(0, 3) == "0|1")
            {
                line_temp[count - 9] = 1;
            }
            else if (char2string.substr(0, 3) == "1|0")
            {
                line_temp[count - 9] = 2;
            }
            else if (char2string.substr(0, 3) == "1|1")
            {
                line_temp[count - 9] = 3;
            }
        }
        count++;
        split_data = strtok(NULL, deliminator);
    }
}

void functions::split_to_MA(vector<string> &line_Data, string line, char delim)
{
    line_Data.clear();
    int count = 0;

    // string coordinates_String = line.substr(gene_Combo.find(delim) + 1);
    while (line.find(delim) != string::npos)
    {
        line_Data.push_back(line.substr(0, line.find(delim)));

        if (count == 10)
        {
            break;
        }

        line = line.substr(line.find(delim) + 1);
        count++;
    }
    // vector<string>().swap(line_Data);
    // char *convert;
    // string capture(line);
    // convert = &capture[0];
    // // cout<<convert;

    // char deliminator[delim.length() + 1];
    // strcpy(deliminator, delim.c_str());

    // char *split_data;
    // split_data = strtok(convert, deliminator);
    // int count = 0;

    // while (split_data != NULL)
    // {
    //     // cout<<split_data<<endl;
    //     string char2string;
    //     char2string.append(split_data);
    //     // cout << char2string << endl;
    //     line_Data.push_back(char2string);
    //     if (count == 10)
    //     {
    //         break;
    //     }
    //     split_data = strtok(NULL, deliminator);
    //     count++;
    // }
}