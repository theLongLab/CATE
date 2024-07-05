#include "functions_library.cuh"

functions_library::functions_library()
{
}

functions_library::functions_library(int tot_Blocks, int tot_ThreadsperBlock, int gpu_Limit, int CPU_cores)
{
    this->tot_Blocks = tot_Blocks;
    this->tot_ThreadsperBlock = tot_ThreadsperBlock;

    this->gpu_Limit = gpu_Limit;
    this->CPU_cores = CPU_cores;
}

functions_library::functions_library(int *tot_Blocks_array, int *tot_ThreadsperBlock_array, int *CUDA_device_IDs, int num_Cuda_devices, int gpu_Limit, int CPU_cores)
{
    this->tot_Blocks_array = tot_Blocks_array;
    this->tot_ThreadsperBlock_array = tot_ThreadsperBlock_array;
    this->CUDA_device_IDs = CUDA_device_IDs;
    this->num_Cuda_devices = num_Cuda_devices;
    this->gpu_Limit = gpu_Limit;
    this->CPU_cores = CPU_cores;

    tot_Blocks = tot_Blocks_array[0];
    tot_ThreadsperBlock = tot_Blocks_array[0];
}

string functions_library::to_Upper_Case(const string &text)
{
    string result = text;

    for (char &c : result)
    {
        c = toupper(c);
    }

    return result;
}

void functions_library::print_Cuda_device(int cuda_ID, int &tot_Blocks, int &tot_ThreadsperBlock)
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
    tot_Blocks = prop.maxBlocksPerMultiProcessor;
    tot_ThreadsperBlock = prop.maxThreadsPerBlock;
    cout << "GPU thread(s) per block\t: " << tot_ThreadsperBlock << endl;
}

void functions_library::print_Cuda_devices(vector<string> cuda_IDs, int *CUDA_device_IDs, int num_Cuda_devices, int *tot_Blocks, int *tot_ThreadsperBlock)
{
    int nDevices;
    cudaGetDeviceCount(&nDevices);

    if (nDevices >= num_Cuda_devices)
    {
        cout << "Properties of selected " << num_Cuda_devices << " CUDA GPU(s):" << endl;
        for (int device = 0; device < num_Cuda_devices; device++)
        {
            cudaDeviceProp prop;
            CUDA_device_IDs[device] = stoi(cuda_IDs[device]);
            cudaGetDeviceProperties(&prop, CUDA_device_IDs[device]);
            cout << "\nGPU number\t: " << CUDA_device_IDs[device] << endl;
            cout << "GPU name\t: " << prop.name << endl;
            size_t l_free = 0;
            size_t l_Total = 0;
            cudaError_t error_id = cudaMemGetInfo(&l_free, &l_Total);
            cout << "GPU memory (GB)\t: " << l_Total / (1000 * 1000 * 1000) << endl;
            cout << "GPU number of multiprocessor(s)\t: " << prop.multiProcessorCount << endl;
            tot_Blocks[device] = prop.maxBlocksPerMultiProcessor;
            cout << "GPU block(s) per multiprocessor\t: " << tot_Blocks[device] << endl;
            tot_ThreadsperBlock[device] = prop.maxThreadsPerBlock;
            cout << "GPU thread(s) per block\t: " << tot_ThreadsperBlock[device] << endl;
        }
    }
    else
    {
        cout << "ERROR: THERE MORE CUDA DEVICES THAN PRESENT HAVE BEEN SELECTED\n";
        cout << "USER HAS SELECTED " << num_Cuda_devices << " BUT THERE IS/ ARE ONLY " << nDevices << " PRESENT IN THE SYSTEM\n";
        exit(-1);
    }
}

float functions_library::beta_Distribution(float &alpha, float &beta, mt19937 &gen)
{
    gamma_distribution<float> gammaAlpha(alpha, 1.0);
    gamma_distribution<float> gammaBeta(beta, 1.0);
    float x = gammaAlpha(gen);
    float y = gammaBeta(gen);

    // Calculate beta distributed random variable using gamma variables
    return x / (x + y);
}

void functions_library::config_Folder(string location, string type_Folder)
{
    filesystem::path folderPath = location;

    if (filesystem::exists(folderPath) && filesystem::is_directory(folderPath))
    {
        cout << type_Folder << " folder: " << location;
    }
    else
    {
        cout << "Creating " << type_Folder << " folder: " << location;
        filesystem::create_directory(location);
    }
    cout << "\n";
}

string functions_library::read_Reference(string file_location, string &header, int &genome_Size)
{
    string reference_Genome = "";

    header = "";

    cout << "Processing reference file: " << file_location << endl;

    fstream reference_File;
    reference_File.open(file_location, ios::in);

    if (reference_File.is_open())
    {
        string line;

        // read header line;
        getline(reference_File, header);
        cout << "Reference genome header: " << header.substr(1) << endl;

        while (getline(reference_File, line))
        {
            reference_Genome.append(line);
        }

        reference_File.close();
    }

    transform(reference_Genome.begin(), reference_Genome.end(), reference_Genome.begin(), ::toupper);
    genome_Size = reference_Genome.length();

    cout << "Reference genome size (bp): " << genome_Size << "\n\n";

    return reference_Genome;
}

vector<string> functions_library::get_Files(string folder, string extension)
{
    vector<string> files;

    if (extension.at(0) != '.')
    {
        extension = "." + extension;
    }

    for (const auto &entry : filesystem::directory_iterator(folder))
    {
        if (entry.is_regular_file())
        {
            // cout << entry.path() << endl;

            string file = entry.path().string();

            if (file.substr(file.find_last_of('.'), file.length()) == extension)
            {
                cout << file << endl;
                files.push_back(file);
            }
        }
    }

    return files;
}

string functions_library::read_Reference(string file_location, int &genome_Size)
{
    string reference_Genome = "";

    string header;

    cout << "Processing reference file: " << file_location << endl;

    fstream reference_File;
    reference_File.open(file_location, ios::in);

    if (reference_File.is_open())
    {
        string line;

        // read header line;
        getline(reference_File, header);
        cout << "Reference genome header: " << header.substr(1) << endl;

        while (getline(reference_File, line))
        {
            reference_Genome.append(line);
        }

        reference_File.close();
    }

    transform(reference_Genome.begin(), reference_Genome.end(), reference_Genome.begin(), ::toupper);
    genome_Size = reference_Genome.length();

    cout << "Reference genome size (bp): " << genome_Size << "\n\n";

    return reference_Genome;
}

void functions_library::split(vector<string> &line_Data, string line, char delim)
{
    line_Data.clear();

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

    line_Data.push_back(line);
}

int **functions_library::create_INT_2D_arrays(int rows, int columns)
{
    int **array_2D = (int **)malloc(rows * sizeof(int *));

    for (int i = 0; i < rows; i++)
    {
        array_2D[i] = (int *)malloc((columns + 1) * sizeof(int));
    }

    if (rows > 0)
    {
        array_2D[0][0] = 0;
    }

    return array_2D;
    free(array_2D);
}

int **functions_library::create_INT_2D_arrays_for_GPU(int rows, int columns)
{
    int **array_2D = (int **)malloc(rows * sizeof(int *));

    for (int i = 0; i < rows; i++)
    {
        array_2D[i] = (int *)malloc((columns + 1) * sizeof(int));
    }

    return array_2D;
    free(array_2D);
}

float **functions_library::create_FLOAT_2D_arrays_for_GPU(int rows, int columns)
{
    float **array_2D = (float **)malloc(rows * sizeof(float *));

    for (int i = 0; i < rows; i++)
    {
        array_2D[i] = (float *)malloc((columns + 1) * sizeof(float));
    }

    return array_2D;
}

float **functions_library::create_FLOAT_2D_arrays(int rows, int columns)
{
    float **array_2D = (float **)malloc(rows * sizeof(float *));

    for (int i = 0; i < rows; i++)
    {
        array_2D[i] = (float *)malloc((columns + 1) * sizeof(float));
    }

    if (rows > 0)
    {
        array_2D[0][0] = 0;
    }

    return array_2D;
}

void functions_library::array_Copy_Float(float **target_2D, float *D1_array, int row, int num_Values)
{
    for (size_t i = 0; i < num_Values; i++)
    {
        target_2D[row][i] = D1_array[i];
    }
}

__global__ void CUDA_normal_distribution(curandState *state, int num_Values, float *values, float mean, float st_deviation, int start_Index)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < num_Values)
    {
        curand_init(clock64(), tid, 0, &state[tid]);

        float u1 = curand_uniform(&state[tid]);
        float u2 = curand_uniform(&state[tid]);
        float z1 = sqrtf(-2.0f * logf(u1)) * cosf(2.0f * M_PI * u2);

        // printf("%d\n", start_Index);

        values[tid + start_Index] = mean + st_deviation * z1;

        tid += blockDim.x * gridDim.x;
    }
}

__global__ void CUDA_poisson_distribution(curandState *state, int num_Values, int *values, float mean, int start_Index)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < num_Values)
    {
        curand_init(clock64(), tid, 0, &state[tid]);

        values[tid + start_Index] = curand_poisson(&state[tid], mean);

        tid += blockDim.x * gridDim.x;
    }
}

// __device__ float rand_uniform(curandState *state)
// {
//     return curand_uniform(state);
// }

// __device__ float rand_normal(curandState *state)
// {
//     float u1 = rand_uniform(state);
//     float u2 = rand_uniform(state);
//     return sqrtf(-2.0f * logf(u1)) * cosf(2.0f * M_PI * u2);
// }

// __device__ float rand_gamma(curandState *state, float alpha, float beta)
// {
//     if (alpha <= 0.0f || beta <= 0.0f)
//     {
//         return 0.0f;
//     }
//     if (alpha < 1.0f)
//     {
//         float u = rand_uniform(state);
//         return rand_gamma(state, alpha + 1.0f, beta) * powf(u, 1.0f / alpha);
//     }
//     float d = alpha - 1.0f / 3.0f;
//     float c = 1.0f / sqrtf(9.0f * d);
//     while (true)
//     {
//         float x = rand_normal(state);
//         float v = 1.0f + c * x;
//         if (v <= 0.0f)
//         {
//             continue;
//         }
//         float u = rand_uniform(state);
//         float xx = x * x;
//         if (u < 1.0f - 0.0331f * xx * xx)
//         {
//             return beta * d * v * v * v;
//         }
//         if (logf(u) < 0.5f * xx + d * (1.0f - v + logf(v)))
//         {
//             return beta * d * v * v * v;
//         }
//     }
// }

// __global__ void CUDA_gamma_distribution(curandState *state, int num_Values, float *values, float shape, float scale, int start_Index)
// {
//     int tid = threadIdx.x + blockIdx.x * blockDim.x;

//     while (tid < num_Values)
//     {
//         curand_init(clock64(), tid, 0, &state[tid]);

//         values[tid + start_Index] = rand_gamma(&state[tid], shape, scale);

//         tid += blockDim.x * gridDim.x;
//     }
// }

// __global__ void CUDA_gamma_distribution_PROGENY(curandState *state, int num_Values, int **Progeny_values, float shape, float scale, int start_Index, int hotspot_Number, float **CUDA_current_gen_Parent_data)
// {
//     int tid = threadIdx.x + blockIdx.x * blockDim.x;

//     while (tid < num_Values)
//     {
//         curand_init(clock64(), tid, 0, &state[tid]);
//         curandState local_state = state[tid];

//         int parent_ID = tid + start_Index;

//         // float progeny_base = rand_gamma(&state[tid], shape, scale);

//         float progeny_Base_fitness = CUDA_current_gen_Parent_data[parent_ID][0];
//         // int fitness_Index = (hotspot * 3) + 2;

//         for (int hotspot = 0; hotspot < hotspot_Number; hotspot++)
//         {
//             progeny_Base_fitness = progeny_Base_fitness * CUDA_current_gen_Parent_data[parent_ID][(hotspot * 3) + 2];
//         }

//         Progeny_values[0][parent_ID] = (int)(rand_gamma(&state[tid], shape, scale) * progeny_Base_fitness);

//         for (int hotspot = 0; hotspot < hotspot_Number; hotspot++)
//         {
//             int ratio_Index = (hotspot * 3) + 3;
//             int count = 0;

//             if (Progeny_values[0][parent_ID] == 0)
//             {
//                 int selectivity_Index = (hotspot * 3) + 1;
//                 CUDA_current_gen_Parent_data[parent_ID][selectivity_Index] = 0;
//             }

//             for (int i = 0; i < Progeny_values[0][parent_ID]; i++)
//             {
//                 if (curand_uniform(&local_state) < CUDA_current_gen_Parent_data[parent_ID][ratio_Index])
//                 {
//                     count++;
//                 }
//             }

//             Progeny_values[hotspot + 1][parent_ID] = count;
//         }

//         state[tid] = local_state;

//         tid += blockDim.x * gridDim.x;
//     }
// }
__global__ void CUDA_binomial_distribution(curandState *state, int num_Values, int *values, float prob, int trials, int start_Index)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < num_Values)
    {
        curand_init(clock64(), tid, 0, &state[tid]);
        curandState local_state = state[tid];

        int count = 0;

        for (int i = 0; i < trials; i++)
        {
            if (curand_uniform(&local_state) < prob)
            {
                count++;
            }
        }

        values[tid + start_Index] = count;

        state[tid] = local_state;

        tid += blockDim.x * gridDim.x;
    }
}

__global__ void CUDA_negative_binomial_distribution(curandState *state, int num_Values, int *values, int r, float p, int start_Index)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < num_Values)
    {
        curand_init(clock64(), tid, 0, &state[tid]);
        curandState local_state = state[tid];

        int k = 0;
        float q = 1 - p;

        while (k < r)
        {
            float u = curand_uniform(&local_state);
            float v = curand_uniform(&local_state);
            float x = ceil(log(u) / log(q));
            float y = 1 - pow(q, x);
            if (v <= y)
            {
                values[tid + start_Index] = k + x;
                break;
            }
            k += x;
        }

        state[tid] = local_state;

        tid += blockDim.x * gridDim.x;
    }
}

__global__ void CUDA_negative_binomial_distribution_PROGENY(curandState *state, int num_Values, int *Progeny_values, int r, float p, int start_Index, float *cuda__Parent_finess)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < num_Values)
    {
        curand_init(clock64(), tid, 0, &state[tid]);
        curandState local_state = state[tid];

        int parent_ID = tid + start_Index;

        int k = 0;
        float q = 1 - p;

        while (k < r)
        {
            float u = curand_uniform(&local_state);
            float v = curand_uniform(&local_state);
            float x = ceil(log(u) / log(q));
            float y = 1 - pow(q, x);
            if (v <= y)
            {
                Progeny_values[parent_ID] = (int)((k + x) * cuda__Parent_finess[parent_ID]);
                break;
            }
            k += x;
        }

        state[tid] = local_state;

        tid += blockDim.x * gridDim.x;
    }
}

int *functions_library::negative_binomial_distribution_CUDA(int num_of_values, float mean, float dispersion_Parameter)
{
    cout << "Negative binomial distribution: Mean (" << mean << ") Dispersion (" << dispersion_Parameter << ")" << endl;

    // dispersion_paramter = r
    float prob = dispersion_Parameter / (mean + dispersion_Parameter);
    int r = round((mean * prob) / (1 - prob));

    cout << "R: " << r << endl;
    cout << "Probability: " << prob << endl;

    int *values, *cuda_Values;
    cudaMallocManaged(&cuda_Values, num_of_values * sizeof(int));
    values = (int *)malloc(num_of_values * sizeof(int));

    int full_Rounds = num_of_values / this->gpu_Limit;
    int partial_Rounds = num_of_values % this->gpu_Limit;

    vector<pair<int, int>> start_stops;

    for (int full = 0; full < full_Rounds; full++)
    {
        int start = full * this->gpu_Limit;
        int stop = start + this->gpu_Limit;
        start_stops.push_back(make_pair(start, stop));
    }

    if (partial_Rounds != 0)
    {
        int start = num_of_values - partial_Rounds;
        start_stops.push_back(make_pair(start, num_of_values));
    }

    for (size_t i = 0; i < start_stops.size(); i++)
    {
        // cout << "Round " << i + 1 << " of " << start_stops.size() << endl;
        int num_of_values_current = start_stops[i].second - start_stops[i].first;

        // cout << start_stops[i].first << "\t" << start_stops[i].second << endl;

        curandState *state;
        cudaMalloc((void **)&state, num_of_values_current * sizeof(curandState));
        // CUDA_negative_binomial_distribution(curandState *state, int num_Values, int *values, float r, float p, int start_Index)
        CUDA_negative_binomial_distribution<<<tot_Blocks, tot_ThreadsperBlock>>>(state, num_of_values_current, cuda_Values, r, prob, start_stops[i].first);

        cudaError_t err = cudaGetLastError();

        if (err != cudaSuccess)
        {
            printf("CUDA Error: %s\n", cudaGetErrorString(err));

            // Possibly: exit(-1) if program cannot continue....
        }
        cudaDeviceSynchronize();

        cudaFree(state);
        // break;
    }

    cudaMemcpy(values, cuda_Values, num_of_values * sizeof(int), cudaMemcpyDeviceToHost);

    cudaFree(cuda_Values);
    cout << "Completed generation via " << start_stops.size() << " GPU rounds" << endl;

    return values;
}

int *functions_library::binomial_distribution_CUDA(int num_of_values, float prob, int progeny_Number)
{
    cout << "Binomial distribution: Probability of sucess (" << prob << ") Trials (" << progeny_Number << ")" << endl;

    int *values, *cuda_Values;
    cudaMallocManaged(&cuda_Values, num_of_values * sizeof(int));
    values = (int *)malloc(num_of_values * sizeof(int));

    int full_Rounds = num_of_values / this->gpu_Limit;
    int partial_Rounds = num_of_values % this->gpu_Limit;

    vector<pair<int, int>> start_stops;

    for (int full = 0; full < full_Rounds; full++)
    {
        int start = full * this->gpu_Limit;
        int stop = start + this->gpu_Limit;
        start_stops.push_back(make_pair(start, stop));
    }

    if (partial_Rounds != 0)
    {
        int start = num_of_values - partial_Rounds;
        start_stops.push_back(make_pair(start, num_of_values));
    }

    for (size_t i = 0; i < start_stops.size(); i++)
    {
        // cout << "Round " << i + 1 << " of " << start_stops.size() << endl;
        int num_of_values_current = start_stops[i].second - start_stops[i].first;

        // cout << start_stops[i].first << "\t" << start_stops[i].second << endl;

        curandState *state;
        cudaMalloc((void **)&state, num_of_values_current * sizeof(curandState));

        CUDA_binomial_distribution<<<tot_Blocks, tot_ThreadsperBlock>>>(state, num_of_values_current, cuda_Values, prob, progeny_Number, start_stops[i].first);

        cudaError_t err = cudaGetLastError();

        if (err != cudaSuccess)
        {
            printf("CUDA Error: %s\n", cudaGetErrorString(err));

            // Possibly: exit(-1) if program cannot continue....
        }
        cudaDeviceSynchronize();

        cudaFree(state);
        // break;
    }

    cudaMemcpy(values, cuda_Values, num_of_values * sizeof(int), cudaMemcpyDeviceToHost);

    cudaFree(cuda_Values);
    cout << "Completed generation via " << start_stops.size() << " GPU rounds" << endl;
    return values;
}

// float *functions_library::gamma_distribution_CUDA(int num_of_values, float shape, float scale)
// {
//     // Marsaglia and Tsang method
//     cout << "Gamma distribution: Shape (" << shape << ") Scale (" << scale << ")" << endl;
//     float *values, *cuda_Values;
//     cudaMallocManaged(&cuda_Values, num_of_values * sizeof(float));
//     values = (float *)malloc(num_of_values * sizeof(float));

//     int full_Rounds = num_of_values / this->gpu_Limit;
//     int partial_Rounds = num_of_values % this->gpu_Limit;

//     vector<pair<int, int>> start_stops;

//     for (int full = 0; full < full_Rounds; full++)
//     {
//         int start = full * this->gpu_Limit;
//         int stop = start + this->gpu_Limit;
//         start_stops.push_back(make_pair(start, stop));
//     }

//     if (partial_Rounds != 0)
//     {
//         int start = num_of_values - partial_Rounds;
//         start_stops.push_back(make_pair(start, num_of_values));
//     }

//     for (size_t i = 0; i < start_stops.size(); i++)
//     {

//         // cout << "Round " << i + 1 << " of " << start_stops.size() << endl;
//         int num_of_values_current = start_stops[i].second - start_stops[i].first;

//         curandState *state;
//         cudaMalloc((void **)&state, num_of_values_current * sizeof(curandState));

//         // cout << start_stops[i].first << "\t" << start_stops[i].second << endl;

//         //CUDA_gamma_distribution<<<tot_Blocks, tot_ThreadsperBlock>>>(state, num_of_values_current, cuda_Values, shape, scale, start_stops[i].first);

//         cudaError_t err = cudaGetLastError();

//         if (err != cudaSuccess)
//         {
//             printf("CUDA Error: %s\n", cudaGetErrorString(err));

//             // Possibly: exit(-1) if program cannot continue....
//         }
//         cudaDeviceSynchronize();
//         // break;

//         cudaFree(state);
//     }

//     cout << "Completed generation via " << start_stops.size() << " GPU rounds" << endl;

//     cudaMemcpy(values, cuda_Values, num_of_values * sizeof(float), cudaMemcpyDeviceToHost);

//     cudaFree(cuda_Values);

//     return values;
// }

int *functions_library::poisson_distribution_CUDA(int num_of_values, float mean)
{
    cout << "Poisson distribution: Lamda (" << mean << ")" << endl;
    int *values, *cuda_Values;
    cudaMallocManaged(&cuda_Values, num_of_values * sizeof(int));
    values = (int *)malloc(num_of_values * sizeof(int));

    int full_Rounds = num_of_values / this->gpu_Limit;
    int partial_Rounds = num_of_values % this->gpu_Limit;

    vector<pair<int, int>> start_stops;

    for (int full = 0; full < full_Rounds; full++)
    {
        int start = full * this->gpu_Limit;
        int stop = start + this->gpu_Limit;
        start_stops.push_back(make_pair(start, stop));
    }

    if (partial_Rounds != 0)
    {
        int start = num_of_values - partial_Rounds;
        start_stops.push_back(make_pair(start, num_of_values));
    }

    for (size_t i = 0; i < start_stops.size(); i++)
    {

        // cout << "Round " << i + 1 << " of " << start_stops.size() << endl;
        int num_of_values_current = start_stops[i].second - start_stops[i].first;

        curandState *state;
        cudaMalloc((void **)&state, num_of_values_current * sizeof(curandState));

        // cout << start_stops[i].first << "\t" << start_stops[i].second << endl;

        CUDA_poisson_distribution<<<tot_Blocks, tot_ThreadsperBlock>>>(state, num_of_values_current, cuda_Values, mean, start_stops[i].first);

        cudaError_t err = cudaGetLastError();

        if (err != cudaSuccess)
        {
            printf("CUDA Error: %s\n", cudaGetErrorString(err));

            // Possibly: exit(-1) if program cannot continue....
        }
        cudaDeviceSynchronize();
        // break;

        cudaFree(state);
    }

    cudaMemcpy(values, cuda_Values, num_of_values * sizeof(int), cudaMemcpyDeviceToHost);

    cudaFree(cuda_Values);
    cout << "Completed generation via " << start_stops.size() << " GPU rounds" << endl;
    return values;
}

int *functions_library::copy_1D_to_CUDA_INT(int *host_Array, int num_Values)
{
    int *cuda_Device_array;
    cudaMallocManaged(&cuda_Device_array, num_Values * sizeof(int));

    cudaMemcpy(cuda_Device_array, host_Array, num_Values * sizeof(int), cudaMemcpyHostToDevice);

    return cuda_Device_array;
}

float *functions_library::copy_1D_to_CUDA_FLOAT(float *host_Array, int num_Values)
{
    float *cuda_Device_array;
    cudaMallocManaged(&cuda_Device_array, num_Values * sizeof(float));

    cudaMemcpy(cuda_Device_array, host_Array, num_Values * sizeof(float), cudaMemcpyHostToDevice);

    return cuda_Device_array;
}

float *functions_library::normal_distribution_CUDA(int num_of_values, float mean, float st_deviation)
{
    cout << "Normal distribution: Mean (" << mean << ") Standard deviation (" << st_deviation << ")" << endl;
    // Uses Box-Muller transformation
    float *values, *cuda_Values;
    cudaMallocManaged(&cuda_Values, num_of_values * sizeof(float));
    values = (float *)malloc(num_of_values * sizeof(float));

    int full_Rounds = num_of_values / this->gpu_Limit;
    int partial_Rounds = num_of_values % this->gpu_Limit;

    vector<pair<int, int>> start_stops;

    for (int full = 0; full < full_Rounds; full++)
    {
        int start = full * this->gpu_Limit;
        int stop = start + this->gpu_Limit;
        start_stops.push_back(make_pair(start, stop));
    }

    if (partial_Rounds != 0)
    {
        int start = num_of_values - partial_Rounds;
        start_stops.push_back(make_pair(start, num_of_values));
    }

    for (size_t i = 0; i < start_stops.size(); i++)
    {

        // cout << "Round " << i + 1 << " of " << start_stops.size() << endl;
        int num_of_values_current = start_stops[i].second - start_stops[i].first;

        curandState *state;
        cudaMalloc((void **)&state, num_of_values_current * sizeof(curandState));

        // cout << start_stops[i].first << "\t" << start_stops[i].second << endl;

        CUDA_normal_distribution<<<tot_Blocks, tot_ThreadsperBlock>>>(state, num_of_values_current, cuda_Values, mean, st_deviation, start_stops[i].first);

        cudaError_t err = cudaGetLastError();

        if (err != cudaSuccess)
        {
            printf("CUDA Error: %s\n", cudaGetErrorString(err));

            // Possibly: exit(-1) if program cannot continue....
        }
        cudaDeviceSynchronize();

        cudaFree(state);
        // break;
    }

    cudaMemcpy(values, cuda_Values, num_of_values * sizeof(float), cudaMemcpyDeviceToHost);

    cudaFree(cuda_Values);
    cout << "Completed generation via " << start_stops.size() << " GPU rounds" << endl;
    return values;
}

int **functions_library::progeny_distribution_CUDA(string &distribution_Type, int &num_of_parents, float &shape, float &scale, float &mean, float &dispersion_Parameter, float *cuda__Parent_finess, float **CUDA_current_gen_Parent_data, int recombination_hotspots)
{
    cout << "Determining progeny numbers" << endl;

    // int rows = this->recombination_hotspots;
    if (recombination_hotspots == -1)
    {
        recombination_hotspots = 0;
    }

    cout << "Generating Progeny numbers " << endl;

    int **cuda_Progeny_numbers;

    cudaMallocManaged(&cuda_Progeny_numbers, (num_of_parents + 1) * (recombination_hotspots + 1) * sizeof(int));

    int **tmp = (int **)malloc((recombination_hotspots + 1) * sizeof(tmp[0]));

    for (int i = 0; i < (recombination_hotspots + 1); i++)
    {
        // cout << i << endl;
        cudaMalloc((void **)&tmp[i], (num_of_parents + 1) * sizeof(tmp[0][0]));
    }
    // cout << "run" << endl;
    cudaMemcpy(cuda_Progeny_numbers, tmp, (recombination_hotspots + 1) * sizeof(int *), cudaMemcpyHostToDevice);
    free(tmp);

    int full_Rounds = num_of_parents / this->gpu_Limit;
    int partial_Rounds = num_of_parents % this->gpu_Limit;

    float prob;
    int r;

    if (distribution_Type == "Negative binomial")
    {
        cout << "Negative binomial distribution: Mean (" << mean << ") Dispersion (" << dispersion_Parameter << ")" << endl;
        prob = dispersion_Parameter / (mean + dispersion_Parameter);
        r = round((mean * prob) / (1 - prob));

        cout << "R: " << r << endl;
        cout << "Probability: " << prob << endl;
    }
    else
    {
        cout << "Gamma distribution: Shape (" << shape << ") Scale (" << scale << ")" << endl;
    }

    vector<pair<int, int>> start_stops;

    for (int full = 0; full < full_Rounds; full++)
    {
        int start = full * this->gpu_Limit;
        int stop = start + this->gpu_Limit;
        start_stops.push_back(make_pair(start, stop));
    }

    if (partial_Rounds != 0)
    {
        int start = num_of_parents - partial_Rounds;
        start_stops.push_back(make_pair(start, num_of_parents));
    }

    for (size_t i = 0; i < start_stops.size(); i++)
    {

        // cout << "Round " << i + 1 << " of " << start_stops.size() << endl;
        int num_of_values_current = start_stops[i].second - start_stops[i].first;

        curandState *state;
        cudaMalloc((void **)&state, num_of_values_current * sizeof(curandState));

        if (distribution_Type == "Gamma")
        {
            // CUDA_gamma_distribution_PROGENY(curandState *state, int num_Values, int **Progeny_values, float shape, float scale, int start_Index, float *cuda__Parent_finess, int hotspot_Number, float **CUDA_current_gen_Parent_data)
            // CUDA_gamma_distribution_PROGENY<<<tot_Blocks, tot_ThreadsperBlock>>>(state, num_of_values_current, cuda_Progeny_numbers, shape, scale, start_stops[i].first, cuda__Parent_finess, recombination_hotspots, CUDA_current_gen_Parent_data);
        }
        else if (distribution_Type == "Negative binomial")
        {
            // CUDA_negative_binomial_distribution_PROGENY<<<tot_Blocks, tot_ThreadsperBlock>>>(state, num_of_values_current, cuda_Progeny_numbers, r, prob, start_stops[i].first, cuda__Parent_finess);
            // CUDA_negative_binomial_distribution_PROGENY
        }

        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess)
        {
            printf("CUDA Error: %s\n", cudaGetErrorString(err));

            // Possibly: exit(-1) if program cannot continue....
        }
        cudaDeviceSynchronize();
        // break;

        cudaFree(state);
    }

    cout << "Completed generation via " << start_stops.size() << " GPU rounds" << endl;

    return cuda_Progeny_numbers;
}

int functions_library::get_base_Index(string base)
{
    if (base == "A" || base == "a")
    {
        return 0;
    }
    else if (base == "T" || base == "t")
    {
        return 1;
    }
    else if (base == "G" || base == "g")
    {
        return 2;
    }
    else if (base == "C" || base == "c")
    {
        return 3;
    }
    else
    {
        cout << "BASES SHOULD BE EITHER A, T , G OR C\n";
        exit(-1);
    }
}

int functions_library::binary_Search(vector<int> &values, int value)
{
    int top = 0;
    int bottom = values.size() - 1;

    while (top <= bottom)
    {
        int mid = top + (bottom - top) / 2;

        if (values[mid] == value)
        {
            return mid;
        }
        else if (values[mid] < value)
        {
            top = mid + 1;
        }
        else
        {
            bottom = mid - 1;
        }
    }

    return -1;
}

void functions_library::get_base_mutation(string query, char &base, int &mutation)
{
    base = toupper(query.at(0));
    char mutation_Query = query.at(1);

    if (mutation_Query == 'A' || mutation_Query == 'a')
    {
        mutation = 0;
    }
    else if (mutation_Query == 'T' || mutation_Query == 't')
    {
        mutation = 1;
    }
    else if (mutation_Query == 'G' || mutation_Query == 'g')
    {
        mutation = 2;
    }
    else if (mutation_Query == 'C' || mutation_Query == 'c')
    {
        mutation = 3;
    }
    else
    {
        cout << "ERROR IN MUTATION QUERY " << mutation_Query << endl;
        exit(-1);
    }
}

__global__ void cuda_Fill_2D_array(int total, int columns, int fill_Value, int **array_2D)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < total)
    {
        int row = tid / columns;
        int column = tid % columns;

        array_2D[row][column] = fill_Value;

        tid += blockDim.x * gridDim.x;
    }
}

__global__ void cuda_Fill_2D_array_Float(int total, int columns, float fill_Value, float **array_2D)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < total)
    {
        int row = tid / columns;
        int column = tid % columns;

        array_2D[row][column] = fill_Value;

        tid += blockDim.x * gridDim.x;
    }
}

__global__ void cuda_Fill_2D_array_INT_Device(int total, int columns, int fill_Value, int **array_2D, int start_Index)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < total)
    {
        int cell_Index = tid + start_Index;
        int row = cell_Index / columns;
        int column = cell_Index % columns;

        array_2D[row][column] = fill_Value;

        tid += blockDim.x * gridDim.x;
    }
}

__global__ void progeny_Array_CUDA(int total_Progeny_current, int parent_ID, int num_Recombination_hotspots, int **cuda_Progeny_numbers, int **cuda_progeny_Array, int fill_Value, int start_Index, float **CUDA_fitness_distribution, int count_Parents)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < total_Progeny_current)
    {
        // int parent_Progeny_Index = progeny_Index + parent_Start;
        // printf("ID: %d start: %d\n", tid, start_Index);
        // if (parent_ID != 1)
        //{
        int progeny_Index = tid + start_Index;
        cuda_progeny_Array[progeny_Index][0] = parent_ID;
        for (int hotspot = 0; hotspot < num_Recombination_hotspots; hotspot++)
        {
            int number_of_recombinant_Progeny = cuda_Progeny_numbers[hotspot + 1][parent_ID];

            if (progeny_Index < number_of_recombinant_Progeny)
            {
                curandState state;
                curand_init(clock64(), tid, 0, &state);
                float rand_num = curand_uniform(&state);
                float cumulative_prob = 0.0f;

                int parent = -1;

                for (int check = 0; check < count_Parents; check++)
                {
                    cumulative_prob += CUDA_fitness_distribution[hotspot][check];
                    if (rand_num < cumulative_prob)
                    {
                        parent = check;
                        break;
                    }
                }

                cuda_progeny_Array[progeny_Index][hotspot + 1] = parent;
            }
            else
            {
                cuda_progeny_Array[progeny_Index][hotspot + 1] = fill_Value;
            }
        }
        //}
        tid += blockDim.x * gridDim.x;
    }
}

__global__ void progeny_Array_CUDA_CELLS(int total_Progeny_current, int parent_Start, int start_Index, int parent_ID, int num_Recombination_hotspots, int **cuda_Progeny_numbers, int **cuda_progeny_Array, int fill_Value, float **CUDA_fitness_distribution, int *CUDA_per_Cell_parents_Stride, int **CUDA_parent_Indexes_parents_and_their_Cells)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < total_Progeny_current)
    {
        // int parent_Progeny_Index = progeny_Index + parent_Start;
        // printf("ID: %d start: %d\n", tid, start_Index);
        // if (parent_ID != 1)
        //{
        int progeny_Index = tid + start_Index;
        int progeny_Array_Index = parent_Start + progeny_Index;

        cuda_progeny_Array[progeny_Array_Index][0] = parent_ID;

        int cell_ID = CUDA_parent_Indexes_parents_and_their_Cells[1][parent_ID];
        int start_Cell = CUDA_per_Cell_parents_Stride[cell_ID];
        int stop_Cell = CUDA_per_Cell_parents_Stride[cell_ID + 1];

        for (int hotspot = 0; hotspot < num_Recombination_hotspots; hotspot++)
        {
            int number_of_recombinant_Progeny = cuda_Progeny_numbers[hotspot + 1][parent_ID];

            if (progeny_Index < number_of_recombinant_Progeny)
            {
                curandState state;
                curand_init(clock64(), tid, 0, &state);
                float rand_num = curand_uniform(&state);
                float cumulative_prob = 0.0f;

                int parent = -1;
                //  int index_Parent = 0;

                for (int check = start_Cell; check < stop_Cell; check++)
                {
                    cumulative_prob += CUDA_fitness_distribution[check][hotspot];
                    if (rand_num < cumulative_prob)
                    {
                        parent = check;
                        break;
                    }
                    // index_Parent++;
                }

                cuda_progeny_Array[progeny_Array_Index][hotspot + 1] = parent;
            }
            else
            {
                cuda_progeny_Array[progeny_Array_Index][hotspot + 1] = fill_Value;
            }
        }
        //}
        tid += blockDim.x * gridDim.x;
    }
}

__global__ void CUDA_summation_Selectivity(int num_Hotspots, int count_Parents, float **CUDA_current_gen_Parent_data, float *CUDA_hotspot_selectivity_Summations)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < num_Hotspots)
    {
        int selectivity_Index = (tid * 3) + 1;
        float total_Selectivity = 0;

        for (int parent = 0; parent < count_Parents; parent++)
        {
            total_Selectivity = total_Selectivity + CUDA_current_gen_Parent_data[parent][selectivity_Index];
        }

        CUDA_hotspot_selectivity_Summations[tid] = total_Selectivity;

        tid += blockDim.x * gridDim.x;
    }
}

__global__ void CUDA_summation_Selectivity_CELLS(int cells, int start_Index, int num_Hotspots, float **CUDA_current_gen_Parent_data, float **CUDA_hotspot_selectivity_Summations, int *CUDA_per_Cell_parents_Stride)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < cells)
    {
        int cell_ID = tid + start_Index;

        int start = CUDA_per_Cell_parents_Stride[cell_ID];
        int stop = CUDA_per_Cell_parents_Stride[cell_ID + 1];

        for (int hotspot = 0; hotspot < num_Hotspots; hotspot++)
        {
            int selectivity_Index = (hotspot * 3) + 1;
            float total_Selectivity = 0;

            for (int parent = start; parent < stop; parent++)
            {
                total_Selectivity = total_Selectivity + CUDA_current_gen_Parent_data[parent][selectivity_Index];
            }
            CUDA_hotspot_selectivity_Summations[cell_ID][hotspot] = total_Selectivity;
        }

        // for (int parent = 0; parent < count_Parents; parent++)
        // {
        //     total_Selectivity = total_Selectivity + CUDA_current_gen_Parent_data[parent][selectivity_Index];
        // }

        // CUDA_hotspot_selectivity_Summations[tid] = total_Selectivity;

        tid += blockDim.x * gridDim.x;
    }
}

__global__ void CUDA_Selectivity_Distribution_CELLS(int num_Hotspots, int count_Parents, float **CUDA_current_gen_Parent_data, float **CUDA_hotspot_selectivity_Summations, float **CUDA_fitness_distribution, int start_Index, int **CUDA_parent_Indexes_parents_and_their_Cells)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < count_Parents)
    {
        int parent_Index = tid + start_Index;
        int cell_ID = CUDA_parent_Indexes_parents_and_their_Cells[1][parent_Index];

        for (int hotspot = 0; hotspot < num_Hotspots; hotspot++)
        {
            int selectivity_Index = (hotspot * 3) + 1;
            CUDA_fitness_distribution[parent_Index][hotspot] = CUDA_current_gen_Parent_data[parent_Index][selectivity_Index] / CUDA_hotspot_selectivity_Summations[cell_ID][hotspot];
        }

        tid += blockDim.x * gridDim.x;
    }
}

__global__ void CUDA_Selectivity_Distribution(int num_Hotspots, int count_Parents, float **CUDA_current_gen_Parent_data, float *CUDA_hotspot_selectivity_Summations, float **CUDA_fitness_distribution, int start_Index)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < count_Parents)
    {
        int parent_Index = tid + start_Index;

        for (int hotspot = 0; hotspot < num_Hotspots; hotspot++)
        {
            int selectivity_Index = (hotspot * 3) + 1;
            CUDA_fitness_distribution[hotspot][parent_Index] = CUDA_current_gen_Parent_data[parent_Index][selectivity_Index] / CUDA_hotspot_selectivity_Summations[hotspot];
        }

        tid += blockDim.x * gridDim.x;
    }
}

__global__ void CUDA_assign_parents_Recombination(curandState *states, int total_Progeny, int **progeny_recom_Index_Cuda, int num_Hotspots, float **cell_hotspot_selectivity_distribution, int *parent_and_their_cell_CUDA, int **cell_and_their_viruses_CUDA, int *per_Cell_max_viruses_CUDA, float **CUDA_current_gen_Parent_data, int start_Index)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < total_Progeny)
    {
        int progeny_Index = tid + start_Index;
        int parent = progeny_recom_Index_Cuda[progeny_Index][0];
        int cell_of_Parent = parent_and_their_cell_CUDA[parent];
        int putative_parents_in_Cell = per_Cell_max_viruses_CUDA[cell_of_Parent];

        for (int recom_Hotspot = 0; recom_Hotspot < num_Hotspots; recom_Hotspot++)
        {
            if (progeny_recom_Index_Cuda[progeny_Index][recom_Hotspot + 1] != -1)
            {
                int hotspot_start = recom_Hotspot * putative_parents_in_Cell;
                int hotspot_end = hotspot_start + putative_parents_in_Cell;

                curand_init(clock64(), tid, 0, &states[tid + recom_Hotspot]);
                float rand_num = curand_uniform(&states[tid + recom_Hotspot]);
                float cumulative_prob = 0.0f;

                int virus_Parent_Index = -1;

                int vp_Index_Count = 0;

                for (int i = hotspot_start; i < hotspot_end; i++)
                {
                    cumulative_prob += cell_hotspot_selectivity_distribution[cell_of_Parent][i];
                    if (rand_num < cumulative_prob)
                    {
                        virus_Parent_Index = vp_Index_Count;
                        break;
                    }
                    vp_Index_Count++;
                }

                progeny_recom_Index_Cuda[progeny_Index][recom_Hotspot + 1] = cell_and_their_viruses_CUDA[cell_of_Parent][virus_Parent_Index];
            }
        }

        tid += blockDim.x * gridDim.x;
    }
}

void functions_library::progeny_Recombination_parents_array(int **progeny_recom_Index_Cuda, int total_Progeny, int num_Hotspots,
                                                            int *parent_and_their_cell_CUDA, int **cell_and_their_viruses_CUDA, int *per_Cell_max_viruses_CUDA,
                                                            float **CUDA_current_gen_Parent_data, int max_Count, int num_Unique_cells)
{
    cout << "Configuring recombinant progeny's parents" << endl;
    float **cell_hotspot_selectivity_distribution_CUDA;

    cudaMallocManaged(&cell_hotspot_selectivity_distribution_CUDA, ((max_Count * num_Hotspots) + 1) * num_Unique_cells * sizeof(int));

    float **tmp = (float **)malloc(num_Unique_cells * sizeof(tmp[0]));

    for (int i = 0; i < num_Unique_cells; i++)
    {
        // cout << i << endl;
        cudaMalloc((void **)&tmp[i], ((max_Count * num_Hotspots) + 1) * sizeof(tmp[0][0]));
    }
    // cout << "run" << endl;
    cudaMemcpy(cell_hotspot_selectivity_distribution_CUDA, tmp, num_Unique_cells * sizeof(float *), cudaMemcpyHostToDevice);
    free(tmp);

    cout << "Creating per cell per recombinant hotspot distributions" << endl;

    int total_Cells = num_Unique_cells;

    int full_Rounds = total_Cells / this->gpu_Limit;
    int partial_Rounds = total_Cells % this->gpu_Limit;

    vector<pair<int, int>> start_stops;

    for (int full = 0; full < full_Rounds; full++)
    {
        int start = full * this->gpu_Limit;
        int stop = start + this->gpu_Limit;
        start_stops.push_back(make_pair(start, stop));
    }

    if (partial_Rounds != 0)
    {
        int start = total_Cells - partial_Rounds;
        start_stops.push_back(make_pair(start, total_Cells));
    }

    for (size_t i = 0; i < start_stops.size(); i++)
    {
        int num_of_values_current = start_stops[i].second - start_stops[i].first;

        // CUDA_distribution_Cells_Selectivity(int num_Cells, int num_Hotspots, int **cell_and_their_viruses_CUDA, int *per_Cell_max_viruses_CUDA, float **CUDA_current_gen_Parent_data, float **cell_hotspot_selectivity_distribution, int start_Index)
        // CUDA_distribution_Cells_Selectivity<<<tot_Blocks, tot_ThreadsperBlock>>>(num_of_values_current, num_Hotspots, cell_and_their_viruses_CUDA, per_Cell_max_viruses_CUDA, CUDA_current_gen_Parent_data, cell_hotspot_selectivity_distribution_CUDA, start_stops[i].first);
        cudaError_t err = cudaGetLastError();

        if (err != cudaSuccess)
        {
            printf("CUDA Error: %s\n", cudaGetErrorString(err));

            // Possibly: exit(-1) if program cannot continue....
        }
        cudaDeviceSynchronize();
    }

    cout << "Distributions created after " << start_stops.size() << " round(s)." << endl;

    // float **test = create_FLOAT_2D_arrays(num_Unique_cells, ((max_Count * num_Hotspots) + 1));
    // // cout << "ok" << endl;
    // for (size_t i = 0; i < num_Unique_cells; i++)
    // {
    //     // cout << i << endl;
    //     cudaMemcpy(test[i], cell_hotspot_selectivity_distribution_CUDA[i], ((max_Count * num_Hotspots) + 1) * sizeof(cell_hotspot_selectivity_distribution_CUDA[0][0]), cudaMemcpyDeviceToHost);
    // }

    // cout << "ok" << endl;

    // for (int row = 0; row < num_Unique_cells; row++)
    // {
    //     for (size_t column = 0; column < 6; column++)
    //     {
    //         cout << test[row][column] << " ";
    //     }
    //     cout << endl;
    // }

    start_stops.clear();

    total_Cells = total_Progeny;

    full_Rounds = total_Cells / this->gpu_Limit;
    partial_Rounds = total_Cells % this->gpu_Limit;

    for (int full = 0; full < full_Rounds; full++)
    {
        int start = full * this->gpu_Limit;
        int stop = start + this->gpu_Limit;
        start_stops.push_back(make_pair(start, stop));
    }

    if (partial_Rounds != 0)
    {
        int start = total_Cells - partial_Rounds;
        start_stops.push_back(make_pair(start, total_Cells));
    }

    for (size_t i = 0; i < start_stops.size(); i++)
    {
        int num_of_values_current = start_stops[i].second - start_stops[i].first;

        curandState *state;
        cudaMalloc((void **)&state, num_of_values_current * num_Hotspots * sizeof(curandState));

        // CUDA_assign_parents_Recombination(curandState *states, int total_Progeny, int **progeny_recom_Index_Cuda, int num_Hotspots, float **cell_hotspot_selectivity_distribution, int *parent_and_their_cell_CUDA, int **cell_and_their_viruses_CUDA, int *per_Cell_max_viruses_CUDA, float **CUDA_current_gen_Parent_data, int start_Index)
        // CUDA_assign_parents_Recombination<<<tot_Blocks, tot_ThreadsperBlock>>>(state, num_of_values_current, progeny_recom_Index_Cuda, num_Hotspots, cell_hotspot_selectivity_distribution_CUDA, parent_and_their_cell_CUDA, cell_and_their_viruses_CUDA, per_Cell_max_viruses_CUDA, CUDA_current_gen_Parent_data, start_stops[i].first);
        cudaError_t err = cudaGetLastError();

        if (err != cudaSuccess)
        {
            printf("CUDA Error: %s\n", cudaGetErrorString(err));

            // Possibly: exit(-1) if program cannot continue....
        }
        cudaDeviceSynchronize();

        cudaFree(state);
    }

    cout << "Recombinant parents assigned after " << start_stops.size() << " round(s)." << endl;
}

__global__ void CUDA_progeny_shuffle(int num_Hotspots, int **progeny_recom_Index_Cuda, int num_Progeny)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < num_Hotspots)
    {

        // curand_init(clock64(), tid, 0, &states[tid]);

        for (int i = num_Progeny - 1; i >= 0; i--)
        {
            curandState state;
            curand_init(clock64(), tid, 0, &state);
            int j = curand(&state) % (i + 1);

            // swap elements at indices i and j
            int temp = progeny_recom_Index_Cuda[i][tid + 1];
            progeny_recom_Index_Cuda[i][tid + 1] = progeny_recom_Index_Cuda[j][tid + 1];
            progeny_recom_Index_Cuda[j][tid + 1] = temp;
        }

        tid += blockDim.x * gridDim.x;
    }
}

__global__ void CUDA_progeny_shuffle_CELLs(int parents, int start_Index, int num_Hotspots, int **progeny_recom_Index_Cuda, int *CUDA_progeny_Stride_Index)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < parents)
    {
        int parent_Index = tid + start_Index;

        int start = CUDA_progeny_Stride_Index[parent_Index];
        int stop = CUDA_progeny_Stride_Index[parent_Index + 1];

        for (int hotspot = 0; hotspot < num_Hotspots; hotspot++)
        {
            curandState state;
            curand_init(clock64(), tid, 0, &state);
            for (int i = start; i < stop - 1; i++)
            {

                int j = curand(&state) % (stop - i) + i;

                int temp = progeny_recom_Index_Cuda[i][hotspot + 1];
                progeny_recom_Index_Cuda[i][hotspot + 1] = progeny_recom_Index_Cuda[j][hotspot + 1];
                progeny_recom_Index_Cuda[j][hotspot + 1] = temp;
            }
        }

        tid += blockDim.x * gridDim.x;
    }
}

void functions_library::progeny_Shuffle(int **progeny_recom_Index_Cuda, int num_Hotspots, int parents_in_current_generation, int *stride_Progeny_Index_CUDA)
{
    cout << "Final configuration of recombinant parents" << endl;

    int total_Cells = parents_in_current_generation * num_Hotspots;

    int full_Rounds = total_Cells / this->gpu_Limit;
    int partial_Rounds = total_Cells % this->gpu_Limit;

    vector<pair<int, int>> start_stops;

    for (int full = 0; full < full_Rounds; full++)
    {
        int start = full * this->gpu_Limit;
        int stop = start + this->gpu_Limit;
        start_stops.push_back(make_pair(start, stop));
    }

    if (partial_Rounds != 0)
    {
        int start = total_Cells - partial_Rounds;
        start_stops.push_back(make_pair(start, total_Cells));
    }

    for (size_t i = 0; i < start_stops.size(); i++)
    {
        int num_of_values_current = start_stops[i].second - start_stops[i].first;

        curandState *state;
        cudaMalloc((void **)&state, num_of_values_current * num_Hotspots * sizeof(curandState));

        // GPU shuffle function
        // CUDA_progeny_shuffle(curandState *states, int iterations, int **progeny_recom_Index_Cuda, int num_Hotspots, int *stride_Progeny_Index_CUDA, int start_Index)
        // CUDA_progeny_shuffle<<<tot_Blocks, tot_ThreadsperBlock>>>(state, num_of_values_current, progeny_recom_Index_Cuda, num_Hotspots, stride_Progeny_Index_CUDA, start_stops[i].first);
        cudaError_t err = cudaGetLastError();

        if (err != cudaSuccess)
        {
            printf("CUDA Error: %s\n", cudaGetErrorString(err));

            // Possibly: exit(-1) if program cannot continue....
        }
        cudaDeviceSynchronize();

        cudaFree(state);
    }

    cout << "Final configuration of recombinant parents completed after " << start_stops.size() << " round(s)." << endl;
}

int **functions_library::create_Progeny_Array(int parents_in_current_generation, int *stride_Progeny_Index, int total_Progeny, int num_Hotspots, int **cuda_Progeny_numbers, int fill_Value)
{
    int **cuda_progeny_Array;
    if (num_Hotspots == -1)
    {
        num_Hotspots = 0;
    }
    cudaMallocManaged(&cuda_progeny_Array, (num_Hotspots + 2) * total_Progeny * sizeof(int));

    int **tmp = (int **)malloc(total_Progeny * sizeof(tmp[0]));

    for (int i = 0; i < total_Progeny; i++)
    {
        // cout << i << endl;
        cudaMalloc((void **)&tmp[i], (num_Hotspots + 2) * sizeof(tmp[0][0]));
    }
    // cout << "run" << endl;
    cudaMemcpy(cuda_progeny_Array, tmp, total_Progeny * sizeof(int *), cudaMemcpyHostToDevice);
    free(tmp);

    int **Progeny_numbers = create_INT_2D_arrays((num_Hotspots + 1), parents_in_current_generation + 1);

    for (size_t i = 0; i < (num_Hotspots + 1); i++)
    {
        // cout << i << endl;
        cudaMemcpy(Progeny_numbers[i], cuda_Progeny_numbers[i], (parents_in_current_generation + 1) * sizeof(cuda_Progeny_numbers[0][0]), cudaMemcpyDeviceToHost);
    }

    // cout << endl;

    int parent_Start = 0;
    stride_Progeny_Index[0] = parent_Start;

    for (int parent = 0; parent < parents_in_current_generation; parent++)
    {
        int total_Cells = Progeny_numbers[0][parent];

        int full_Rounds = total_Cells / this->gpu_Limit;
        int partial_Rounds = total_Cells % this->gpu_Limit;

        vector<pair<int, int>> start_stops;

        for (int full = 0; full < full_Rounds; full++)
        {
            int start = full * this->gpu_Limit;
            int stop = start + this->gpu_Limit;
            start_stops.push_back(make_pair(start, stop));
        }

        if (partial_Rounds != 0)
        {
            int start = total_Cells - partial_Rounds;
            start_stops.push_back(make_pair(start, total_Cells));
        }

        for (size_t i = 0; i < start_stops.size(); i++)
        {
            int num_of_values_current = start_stops[i].second - start_stops[i].first;
            // progeny_Array_CUDA(int total_Progeny, int parent_ID, int parent_Start, int num_Recombination_hotspots, int **cuda_Progeny_numbers, int **cuda_progeny_Array, int fill_Value, int start_Index)
            // progeny_Array_CUDA<<<tot_Blocks, tot_ThreadsperBlock>>>(num_of_values_current, parent, parent_Start, num_Hotspots, cuda_Progeny_numbers, cuda_progeny_Array, fill_Value, start_stops[i].first);
            cudaError_t err = cudaGetLastError();

            if (err != cudaSuccess)
            {
                printf("CUDA Error: %s\n", cudaGetErrorString(err));

                // Possibly: exit(-1) if program cannot continue....
            }
            cudaDeviceSynchronize();
        }

        cout << "Created progeny index for parent " << parent + 1 << " after " << start_stops.size() << " round(s)." << endl;
        parent_Start = parent_Start + total_Cells;
        stride_Progeny_Index[parent + 1] = parent_Start;
    }

    // cout << "Total test: " << parent_Start << endl;

    free(Progeny_numbers);

    return cuda_progeny_Array;
}

__global__ void CUDA_progeny_Profiles_fill(int total_Progeny, float **CUDA_current_gen_Progeny_data, int num_Hotspots, int **progeny_recom_Index_Cuda, float **CUDA_current_gen_Parent_data, float *CUDA_parent_Proof_reading_probability, float *CUDA_progeny_Proof_reading_probability, int start_Index, int proof_reading_Activate_parent)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < total_Progeny)
    {
        int progeny_Index = tid + start_Index;
        int parent = progeny_recom_Index_Cuda[progeny_Index][0];

        CUDA_current_gen_Progeny_data[tid][0] = CUDA_current_gen_Parent_data[parent][0];
        if (proof_reading_Activate_parent != 0)
        {
            CUDA_progeny_Proof_reading_probability[tid] = CUDA_parent_Proof_reading_probability[parent];
        }
        for (int hotspot = 0; hotspot < num_Hotspots; hotspot++)
        {
            int hotspot_Parent = progeny_recom_Index_Cuda[progeny_Index][hotspot + 1];
            if (hotspot_Parent == -1)
            {
                hotspot_Parent = parent;
            }

            int hotspot_start = (hotspot * 3) + 1;
            int hotspot_stop = hotspot_start + 3;
            for (int i = hotspot_start; i < hotspot_stop; i++)
            {
                CUDA_current_gen_Progeny_data[tid][i] = CUDA_current_gen_Parent_data[hotspot_Parent][i];
            }
        }

        tid += blockDim.x * gridDim.x;
    }
}

__global__ void CUDA_progeny_Profiles_fill_CELLs(int total_Progeny, float **CUDA_current_gen_Progeny_data, int num_Hotspots, int **progeny_recom_Index_Cuda, float **CUDA_current_gen_Parent_data, float *CUDA_parent_Proof_reading_probability, float *CUDA_progeny_Proof_reading_probability, int start_Index, int proof_reading_Activate_parent, float **CUDA_parent_survivability_Probabilities, float **CUDA_progeny_survivability_Probabilities)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < total_Progeny)
    {
        int progeny_Index = tid + start_Index;
        int parent = progeny_recom_Index_Cuda[progeny_Index][0];

        CUDA_current_gen_Progeny_data[progeny_Index][0] = CUDA_current_gen_Parent_data[parent][0];
        if (proof_reading_Activate_parent != 0)
        {
            CUDA_progeny_Proof_reading_probability[progeny_Index] = CUDA_parent_Proof_reading_probability[parent];
        }

        CUDA_progeny_survivability_Probabilities[progeny_Index][0] = CUDA_parent_survivability_Probabilities[parent][0];

        for (int hotspot = 0; hotspot < num_Hotspots; hotspot++)
        {
            int hotspot_Parent = progeny_recom_Index_Cuda[progeny_Index][hotspot + 1];
            if (hotspot_Parent == -1)
            {
                hotspot_Parent = parent;
            }

            CUDA_progeny_survivability_Probabilities[progeny_Index][hotspot + 1] = CUDA_parent_survivability_Probabilities[hotspot_Parent][hotspot + 1];

            int hotspot_start = (hotspot * 3) + 1;
            int hotspot_stop = hotspot_start + 3;
            for (int i = hotspot_start; i < hotspot_stop; i++)
            {
                CUDA_current_gen_Progeny_data[progeny_Index][i] = CUDA_current_gen_Parent_data[hotspot_Parent][i];
            }
        }

        tid += blockDim.x * gridDim.x;
    }
}

float **functions_library::create_current_Progeny_data(int total_Progeny, int num_Hotspots, int **progeny_recom_Index_Cuda, float **CUDA_current_gen_Parent_data, float *CUDA_parent_Proof_reading_probability, float *CUDA_progeny_Proof_reading_probability)
{
    cout << "Configuring progeny profiles" << endl;
    float **CUDA_current_gen_Progeny_data;

    cudaMallocManaged(&CUDA_progeny_Proof_reading_probability, total_Progeny * sizeof(float));

    // cudaMallocManaged(&CUDA_progeny_Proof_reading_probability, total_Progeny * sizeof(float));

    cudaMallocManaged(&CUDA_current_gen_Progeny_data, ((1 + (3 * num_Hotspots)) + 1) * total_Progeny * sizeof(float));
    float **tmp = (float **)malloc(total_Progeny * sizeof(tmp[0]));

    for (int i = 0; i < total_Progeny; i++)
    {
        // cout << i << endl;
        cudaMalloc((void **)&tmp[i], ((1 + (3 * num_Hotspots)) + 1) * sizeof(tmp[0][0]));
    }
    // cout << "run" << endl;
    cudaMemcpy(CUDA_current_gen_Progeny_data, tmp, total_Progeny * sizeof(float *), cudaMemcpyHostToDevice);
    free(tmp);

    int total_Cells = total_Progeny;

    int full_Rounds = total_Cells / this->gpu_Limit;
    int partial_Rounds = total_Cells % this->gpu_Limit;

    vector<pair<int, int>> start_stops;

    for (int full = 0; full < full_Rounds; full++)
    {
        int start = full * this->gpu_Limit;
        int stop = start + this->gpu_Limit;
        start_stops.push_back(make_pair(start, stop));
    }

    if (partial_Rounds != 0)
    {
        int start = total_Cells - partial_Rounds;
        start_stops.push_back(make_pair(start, total_Cells));
    }

    // cout << start_stops.size();

    // cout << total_Progeny << endl;

    for (size_t i = 0; i < start_stops.size(); i++)
    {
        int num_of_values_current = start_stops[i].second - start_stops[i].first;
        // cout << num_of_values_current << "\t" << start_stops[i].second << "\t" << start_stops[i].first << endl;
        // CUDA_progeny_Profiles_fill(int total_Progeny, float **CUDA_current_gen_Progeny_data, int num_Hotspots, int **progeny_recom_Index_Cuda, float **CUDA_current_gen_Parent_data, int start_Index, float *CUDA_parent_Proof_reading_probability, float *CUDA_progeny_Proof_reading_probability)
        // CUDA_progeny_Profiles_fill<<<tot_Blocks, tot_ThreadsperBlock>>>(num_of_values_current, CUDA_current_gen_Progeny_data, num_Hotspots, progeny_recom_Index_Cuda, CUDA_current_gen_Parent_data, start_stops[i].first, CUDA_parent_Proof_reading_probability, CUDA_progeny_Proof_reading_probability);

        cudaError_t err = cudaGetLastError();

        if (err != cudaSuccess)
        {
            printf("CUDA Error: %s\n", cudaGetErrorString(err));
            // cout << i << endl;
            exit(-1);
            // Possibly: exit(-1) if program cannot continue....
        }
        // else
        // {
        //     cout << "Sucess " << i << endl;
        // }
        cudaDeviceSynchronize();
    }

    cout << "Configured progeny profiles after " << start_stops.size() << " round(s)." << endl;

    return CUDA_current_gen_Progeny_data;
}

__global__ void CUDA_progeny_sequence_generation(int total_Progeny, int start_Index, int **cuda_parent_sequences, int **cuda_progeny_Sequences, int **progeny_recom_Index_Cuda, int num_Hotspots, int genome_SIZE, int **CUDA_recombination_hotspots_start_stop)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < total_Progeny)
    {
        int progeny_Index = tid + start_Index;
        int parent = progeny_recom_Index_Cuda[progeny_Index][0];

        // if (parent_Sequence_Index == -1)
        // {
        //     printf("%d %d %d %d \n", parent, CUDA_parent_Indexes[0], num_Unique_Parents, parent_Sequence_Index);
        // }

        for (int base = 0; base < genome_SIZE; base++)
        {
            cuda_progeny_Sequences[tid][base] = cuda_parent_sequences[parent][base];
        }

        for (int hotspot = 0; hotspot < num_Hotspots; hotspot++)
        {
            int hotspot_Parent = progeny_recom_Index_Cuda[progeny_Index][hotspot + 1];

            if (hotspot_Parent != -1 && hotspot_Parent != parent)
            {

                int start_Hotspot = CUDA_recombination_hotspots_start_stop[hotspot][0] - 1;
                int stop_Hotspot = CUDA_recombination_hotspots_start_stop[hotspot][1];

                for (int base = start_Hotspot; base < stop_Hotspot; base++)
                {
                    cuda_progeny_Sequences[tid][base] = cuda_parent_sequences[hotspot_Parent][base];
                }

                // if (parent_Sequence_Index == -1)
                // {
                //     printf("2: %d\n", parent_Sequence_Index);
                // }
            }
        }

        tid += blockDim.x * gridDim.x;
    }
}

__global__ void CUDA_progeny_sequence_generation_CELLS(int total_Progeny, int start_Index, int **cuda_parent_sequences, int **cuda_progeny_Sequences, int **progeny_recom_Index_Cuda, int num_Hotspots, int genome_SIZE, int **CUDA_recombination_hotspots_start_stop)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < total_Progeny)
    {
        int progeny_Index = tid + start_Index;
        int parent = progeny_recom_Index_Cuda[progeny_Index][0];

        // if (parent_Sequence_Index == -1)
        // {
        //     printf("%d %d %d %d \n", parent, CUDA_parent_Indexes[0], num_Unique_Parents, parent_Sequence_Index);
        // }

        for (int base = 0; base < genome_SIZE; base++)
        {
            cuda_progeny_Sequences[progeny_Index][base] = cuda_parent_sequences[parent][base];
        }

        for (int hotspot = 0; hotspot < num_Hotspots; hotspot++)
        {
            int hotspot_Parent = progeny_recom_Index_Cuda[progeny_Index][hotspot + 1];

            if (hotspot_Parent != -1 && hotspot_Parent != parent)
            {

                int start_Hotspot = CUDA_recombination_hotspots_start_stop[hotspot][0] - 1;
                int stop_Hotspot = CUDA_recombination_hotspots_start_stop[hotspot][1];

                for (int base = start_Hotspot; base < stop_Hotspot; base++)
                {
                    cuda_progeny_Sequences[progeny_Index][base] = cuda_parent_sequences[hotspot_Parent][base];
                }

                // if (parent_Sequence_Index == -1)
                // {
                //     printf("2: %d\n", parent_Sequence_Index);
                // }
            }
        }

        tid += blockDim.x * gridDim.x;
    }
}

int **functions_library::create_progeny_Sequences(int **cuda_parent_Sequences, int **progeny_recom_Index_Cuda, int num_Hotspots, int total_Progeny, int genome_SIZE, int **CUDA_recombination_hotspots_start_stop)
{
    cout << "Generating progeny sequences" << endl;

    int **cuda_progeny_Sequences;
    cudaMallocManaged(&cuda_progeny_Sequences, (genome_SIZE + 1) * total_Progeny * sizeof(int));
    int **tmp = (int **)malloc(total_Progeny * sizeof(tmp[0]));
    for (int i = 0; i < total_Progeny; i++)
    {
        cudaMalloc((void **)&tmp[i], (genome_SIZE + 1) * sizeof(tmp[0][0]));
    }
    cudaMemcpy(cuda_progeny_Sequences, tmp, total_Progeny * sizeof(int *), cudaMemcpyHostToDevice);
    free(tmp);

    int total_Cells = total_Progeny;

    int full_Rounds = total_Cells / this->gpu_Limit;
    int partial_Rounds = total_Cells % this->gpu_Limit;

    vector<pair<int, int>> start_stops;

    for (int full = 0; full < full_Rounds; full++)
    {
        int start = full * this->gpu_Limit;
        int stop = start + this->gpu_Limit;
        start_stops.push_back(make_pair(start, stop));
    }

    if (partial_Rounds != 0)
    {
        int start = total_Cells - partial_Rounds;
        start_stops.push_back(make_pair(start, total_Cells));
    }
    for (size_t i = 0; i < start_stops.size(); i++)
    {
        int num_of_values_current = start_stops[i].second - start_stops[i].first;

        // CUDA_progeny_sequence_generation(int total_Progeny, int start_Index, int **cuda_parent_sequences, int **cuda_progeny_Sequences, int **progeny_recom_Index_Cuda, int num_Hotspots, int genome_SIZE, int **CUDA_recombination_hotspots_start_stop)
        // CUDA_progeny_sequence_generation<<<tot_Blocks, tot_ThreadsperBlock>>>(num_of_values_current, start_stops[i].first, cuda_parent_Sequences, cuda_progeny_Sequences, progeny_recom_Index_Cuda, num_Hotspots, genome_SIZE, CUDA_recombination_hotspots_start_stop);
        cudaError_t err = cudaGetLastError();

        if (err != cudaSuccess)
        {
            printf("CUDA Error: %s\n", cudaGetErrorString(err));

            // Possibly: exit(-1) if program cannot continue....
        }
        cudaDeviceSynchronize();
    }

    cout << "Generated progeny sequences after " << start_stops.size() << " round(s)." << endl;

    return cuda_progeny_Sequences;
}

int **functions_library::Fill_2D_array_CUDA(int rows, int columns, int fill_Value, int **cuda_Progeny_numbers)
{
    int **cuda_Array_2D;

    cudaMallocManaged(&cuda_Array_2D, rows * sizeof(int *));
    for (int row = 0; row < rows; row++)
    {
        cudaMalloc((void **)&(cuda_Array_2D[row]), (1 + columns) * sizeof(int));
    }

    // cudaMallocManaged(&cuda_Array_2D, (columns + 1) * rows * sizeof(int));

    // int **tmp = (int **)malloc(rows * sizeof(tmp[0]));

    // for (int i = 0; i < rows; i++)
    // {
    //     // cout << i << endl;
    //     cudaMalloc((void **)&tmp[i], (columns + 1) * sizeof(tmp[0][0]));
    // }
    // // cout << "run" << endl;
    // cudaMemcpy(cuda_Array_2D, tmp, rows * sizeof(int *), cudaMemcpyHostToDevice);
    // free(tmp);

    int total_Cells = rows * columns;

    int full_Rounds = total_Cells / this->gpu_Limit;
    int partial_Rounds = total_Cells % this->gpu_Limit;

    vector<pair<int, int>> start_stops;

    for (int full = 0; full < full_Rounds; full++)
    {
        int start = full * this->gpu_Limit;
        int stop = start + this->gpu_Limit;
        start_stops.push_back(make_pair(start, stop));
    }

    if (partial_Rounds != 0)
    {
        int start = total_Cells - partial_Rounds;
        start_stops.push_back(make_pair(start, total_Cells));
    }

    // cout << start_stops.size();

    for (size_t i = 0; i < start_stops.size(); i++)
    {
        int num_of_values_current = start_stops[i].second - start_stops[i].first;

        // cout << start_stops[i].first << "\t" << start_stops[i].second << endl;

        // cuda_Fill_2D_array_INT_Device(int total, int columns, int fill_Value, int **array_2D, int start_Index)
        cuda_Fill_2D_array_INT_Device<<<tot_Blocks, tot_ThreadsperBlock>>>(num_of_values_current, columns, fill_Value, cuda_Array_2D, start_stops[i].first);
        cudaError_t err = cudaGetLastError();

        if (err != cudaSuccess)
        {
            printf("CUDA Error: %s\n", cudaGetErrorString(err));

            // Possibly: exit(-1) if program cannot continue....
        }
        cudaDeviceSynchronize();
    }

    cout << "Created recombination index after " << start_stops.size() << " round(s)." << endl;

    return cuda_Array_2D;
}

int **functions_library::create_Fill_2D_array(int rows, int columns, int fill_Value)
{
    cudaSetDevice(CUDA_device_IDs[0]);
    int **array_2D;
    // create_INT_2D_arrays(rows, (columns + 1));
    array_2D = (int **)malloc(rows * sizeof(int *));
    for (int row = 0; row < rows; row++)
    {
        array_2D[row] = (int *)malloc((columns + 1) * sizeof(int));
    }

    int **cuda_Array_2D;

    cudaMallocManaged(&cuda_Array_2D, rows * sizeof(int *));
    for (int row = 0; row < rows; row++)
    {
        cudaMalloc((void **)&(cuda_Array_2D[row]), (1 + columns) * sizeof(int));
    }

    // cudaMallocManaged(&cuda_Array_2D, (columns + 1) * rows * sizeof(int));

    // int **tmp = (int **)malloc(rows * sizeof(tmp[0]));

    // for (int i = 0; i < rows; i++)
    // {
    //     // cout << i << endl;
    //     cudaMalloc((void **)&tmp[i], (columns + 1) * sizeof(tmp[0][0]));
    // }
    // // cout << "run" << endl;
    // cudaMemcpy(cuda_Array_2D, tmp, rows * sizeof(int *), cudaMemcpyHostToDevice);
    // free(tmp);

    cuda_Fill_2D_array<<<tot_Blocks, tot_ThreadsperBlock>>>((rows * columns), columns, fill_Value, cuda_Array_2D);

    cudaError_t err = cudaGetLastError();

    if (err != cudaSuccess)
    {
        printf("CUDA Error: %s\n", cudaGetErrorString(err));
        exit(-1);
        // Possibly: exit(-1) if program cannot continue....
    }
    cudaDeviceSynchronize();
    //  cout << "run" << endl;

    // cudaMemcpy(array_2D, cuda_Array_2D, rows * sizeof(int), cudaMemcpyDeviceToHost);

    for (size_t i = 0; i < rows; i++)
    {
        // cout << i << endl;
        cudaMemcpy(array_2D[i], cuda_Array_2D[i], (columns + 1) * sizeof(int), cudaMemcpyDeviceToHost);
    }
    // cout << "run" << endl;
    for (int row = 0; row < rows; row++)
    {
        cudaFree(cuda_Array_2D[row]);
    }

    // Free the array of pointers
    cudaFree(cuda_Array_2D);

    return array_2D;
}

float **functions_library::float_2D_Array_load_to_CUDA(float **host_Array, int rows, int columns)
{
    float **cuda_Array_2D;
    // cudaMallocManaged(&cuda_Array_2D, (columns + 1) * rows * sizeof(float));
    // float **tmp = (float **)malloc(rows * sizeof(tmp[0]));

    // for (int i = 0; i < rows; i++)
    // {
    //     // cout << i << endl;
    //     cudaMalloc((void **)&tmp[i], (columns + 1) * sizeof(tmp[0][0]));
    // }
    // // cout << "run" << endl;
    // cudaMemcpy(cuda_Array_2D, tmp, rows * sizeof(float *), cudaMemcpyHostToDevice);

    // for (size_t i = 0; i < rows; i++)
    // {
    //     cudaMemcpy(tmp[i], host_Array[i], (columns + 1) * sizeof(cuda_Array_2D[0][0]), cudaMemcpyHostToDevice);
    // }

    // free(tmp);

    cudaMallocManaged(&cuda_Array_2D, rows * sizeof(float *));
    for (int row = 0; row < rows; row++)
    {
        cudaMalloc((void **)&(cuda_Array_2D[row]), columns * sizeof(float));
    }

    for (int row = 0; row < rows; row++)
    {
        cudaMemcpy(cuda_Array_2D[row], host_Array[row], columns * sizeof(float), cudaMemcpyHostToDevice);
    }

    return cuda_Array_2D;
}

int **functions_library::int_2D_Array_load_to_CUDA(int **host_Array, int rows, int columns)
{
    int **cuda_Array_2D;
    // cudaMallocManaged(&cuda_Array_2D, (columns + 1) * rows * sizeof(int));
    // int **tmp = (int **)malloc(rows * sizeof(tmp[0]));

    // for (int i = 0; i < rows; i++)
    // {
    //     // cout << i << endl;
    //     cudaMalloc((void **)&tmp[i], (columns + 1) * sizeof(tmp[0][0]));
    // }
    // // cout << "run" << endl;
    // cudaMemcpy(cuda_Array_2D, tmp, rows * sizeof(int *), cudaMemcpyHostToDevice);

    // for (size_t i = 0; i < rows; i++)
    // {
    //     cudaMemcpy(tmp[i], host_Array[i], (columns + 1) * sizeof(cuda_Array_2D[0][0]), cudaMemcpyHostToDevice);
    // }

    // free(tmp);

    cudaMallocManaged(&cuda_Array_2D, rows * sizeof(int *));
    for (int row = 0; row < rows; row++)
    {
        cudaMalloc((void **)&(cuda_Array_2D[row]), columns * sizeof(int));
    }

    for (int row = 0; row < rows; row++)
    {
        cudaMemcpy(cuda_Array_2D[row], host_Array[row], columns * sizeof(int), cudaMemcpyHostToDevice);
    }

    return cuda_Array_2D;
}

float **functions_library::create_Fill_2D_array_FLOAT(int rows, int columns, float fill_Value)
{
    cudaSetDevice(CUDA_device_IDs[0]);
    float **array_2D;
    // = create_FLOAT_2D_arrays(rows, (columns + 1));
    array_2D = (float **)malloc(rows * sizeof(float *));
    for (int row = 0; row < rows; row++)
    {
        array_2D[row] = (float *)malloc((columns + 1) * sizeof(float));
    }

    float **cuda_Array_2D;

    cudaMallocManaged(&cuda_Array_2D, rows * sizeof(float *));
    for (int row = 0; row < rows; row++)
    {
        cudaMalloc((void **)&(cuda_Array_2D[row]), (1 + columns) * sizeof(float));
    }

    // cudaMallocManaged(&cuda_Array_2D, (columns + 1) * rows * sizeof(float));

    // float **tmp = (float **)malloc(rows * sizeof(tmp[0]));

    // for (int i = 0; i < rows; i++)
    // {
    //     // cout << i << endl;
    //     cudaMalloc((void **)&tmp[i], (columns + 1) * sizeof(tmp[0][0]));
    // }
    // // cout << "run" << endl;
    // cudaMemcpy(cuda_Array_2D, tmp, rows * sizeof(float *), cudaMemcpyHostToDevice);
    // free(tmp);

    cuda_Fill_2D_array_Float<<<tot_Blocks, tot_ThreadsperBlock>>>((rows * columns), columns, fill_Value, cuda_Array_2D);

    cudaError_t err = cudaGetLastError();

    if (err != cudaSuccess)
    {
        printf("CUDA Error: %s\n", cudaGetErrorString(err));
        exit(-1);
        // Possibly: exit(-1) if program cannot continue....
    }
    cudaDeviceSynchronize();
    //  cout << "run" << endl;

    // cudaMemcpy(array_2D, cuda_Array_2D, rows * sizeof(int), cudaMemcpyDeviceToHost);

    for (size_t i = 0; i < rows; i++)
    {
        // cout << i << endl;
        cudaMemcpy(array_2D[i], cuda_Array_2D[i], (columns + 1) * sizeof(float), cudaMemcpyDeviceToHost);
    }
    // cout << "run" << endl;
    // cudaFree(cuda_Array_2D);

    for (int row = 0; row < rows; row++)
    {
        cudaFree(cuda_Array_2D[row]);
    }

    // Free the array of pointers
    cudaFree(cuda_Array_2D);

    return array_2D;
}

__global__ void CUDA_sum(int **input, int *output, int size)
{
    __shared__ float shared[1024];
    int tid = threadIdx.x;
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    shared[tid] = (i < size) ? input[0][i] : 0;
    __syncthreads();

    for (int s = blockDim.x / 2; s > 0; s >>= 1)
    {
        if (tid < s)
        {
            shared[tid] += shared[tid + s];
        }
        __syncthreads();
    }

    if (tid == 0)
    {
        output[blockIdx.x] = shared[0];
    }
}

int functions_library::sum_CUDA(int **cuda_Array_input, int num_Elements)
{
    int *host_output = (int *)malloc(num_Elements * sizeof(int));
    int *CUDA_output;

    cudaMalloc((void **)&CUDA_output, num_Elements * sizeof(float));

    int threads_per_block = 1024;
    int blocks_per_grid = ceil(float(num_Elements) / threads_per_block);

    CUDA_sum<<<blocks_per_grid, threads_per_block>>>(cuda_Array_input, CUDA_output, num_Elements);
    cudaError_t err = cudaGetLastError();

    if (err != cudaSuccess)
    {
        printf("CUDA Error: %s\n", cudaGetErrorString(err));

        // Possibly: exit(-1) if program cannot continue....
    }
    cudaDeviceSynchronize();

    cudaMemcpy(host_output, CUDA_output, blocks_per_grid * sizeof(int), cudaMemcpyDeviceToHost);

    // Final reduction on CPU

    int sum_Total = 0;

    for (int i = 0; i < blocks_per_grid; i++)
    {
        sum_Total += host_output[i];
    }

    cudaFree(CUDA_output);
    free(host_output);

    return sum_Total;
}

__global__ void CUDA_mutate_Progeny_CELLS(int total_Progeny, int num_mutation_Hotspots, int generation,
                                          float **CUDA_mutation_rates_Hotspot_generation, int **CUDA_mutation_Regions_start_stop,
                                          int **cuda_progeny_Sequences, float *CUDA_progeny_Proof_reading_probability, int proof_reading_Activate,
                                          float **CUDA_A_0_mutation, float **CUDA_T_1_mutation, float **CUDA_G_2_mutation, float **CUDA_C_3_mutation,
                                          int **CUDA_sequence_Mutation_tracker, float **CUDA_current_gen_Progeny_data,
                                          float **CUDA_A_0_fitness, float **CUDA_T_1_fitness, float **CUDA_G_2_fitness, float **CUDA_C_3_fitness,
                                          float **CUDA_A_0_probability_Proof_reading, float **CUDA_T_1_probability_Proof_reading, float **CUDA_G_2_probability_Proof_reading, float **CUDA_C_3_probability_Proof_reading,
                                          float **CUDA_A_0_Recombination, float **CUDA_T_1_Recombination, float **CUDA_G_2_Recombination, float **CUDA_C_3_Recombination,
                                          int *CUDA_stride_Array, int start_Index,
                                          float **CUDA_progeny_survivability_Probabilities,
                                          float **CUDA_A_0_survivability, float **CUDA_T_1_survivability, float **CUDA_G_2_survivability, float **CUDA_C_3_survivability)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < total_Progeny)
    {
        int progeny_Index = tid + start_Index;

        for (int mutation_Hotspot = 0; mutation_Hotspot < num_mutation_Hotspots; mutation_Hotspot++)
        {
            float mean = CUDA_mutation_rates_Hotspot_generation[mutation_Hotspot][generation];
            int start = CUDA_mutation_Regions_start_stop[mutation_Hotspot][0] - 1;
            int stop = CUDA_mutation_Regions_start_stop[mutation_Hotspot][1] - 1;

            curandState state;
            curand_init(clock64(), tid, 0, &state);
            int num_of_Mutations = curand_poisson(&state, mean);

            if (proof_reading_Activate == 1)
            {
                float proof_Reading_probability = CUDA_progeny_Proof_reading_probability[progeny_Index];

                int corrections = 0;

                for (int i = 0; i < num_of_Mutations; i++)
                {
                    if (curand_uniform(&state) < proof_Reading_probability)
                    {
                        corrections++;
                    }
                }

                num_of_Mutations = num_of_Mutations - corrections;

                // if (num_of_Mutations < 0)
                // {
                //     num_of_Mutations = 0;
                // }
            }

            for (int mutation_Number = 0; mutation_Number < num_of_Mutations; mutation_Number++)
            {
                curand_init(clock64(), tid, 0, &state);
                int position = (int)(curand_uniform(&state) * (stop - start + 1)) + start;

                // now do the mutations
                int original_BASE = cuda_progeny_Sequences[progeny_Index][position];
                int new_BASE = -1;

                float rand_num = curand_uniform(&state);
                float cumulative_prob = 0.0f;

                int probable_Base = 0;

                for (int base = 0; base < 4; base++)
                {

                    cumulative_prob += (original_BASE == 0)   ? CUDA_A_0_mutation[mutation_Hotspot][base]
                                       : (original_BASE == 1) ? CUDA_T_1_mutation[mutation_Hotspot][base]
                                       : (original_BASE == 2) ? CUDA_G_2_mutation[mutation_Hotspot][base]
                                       : (original_BASE == 3) ? CUDA_C_3_mutation[mutation_Hotspot][base]
                                                              : 0.0f;

                    // if (original_BASE == 0)
                    // {
                    //     cumulative_prob += CUDA_A_0_mutation[mutation_Hotspot][base];
                    // }
                    // else if (original_BASE == 1)
                    // {
                    //     cumulative_prob += CUDA_T_1_mutation[mutation_Hotspot][base];
                    // }
                    // else if (original_BASE == 2)
                    // {
                    //     cumulative_prob += CUDA_G_2_mutation[mutation_Hotspot][base];
                    // }
                    // else if (original_BASE == 3)
                    // {
                    //     cumulative_prob += CUDA_C_3_mutation[mutation_Hotspot][base];
                    // }

                    if (rand_num < cumulative_prob)
                    {
                        new_BASE = probable_Base;
                        break;
                    }

                    probable_Base++;
                }

                if (original_BASE != new_BASE)
                {
                    // if mutation has occured
                    // change base
                    // see what changes have to be done to sequence profile
                    // dont forget to check proof reading
                    // do fitness and proof reading first

                    // CHANGE to new base

                    cuda_progeny_Sequences[progeny_Index][position] = new_BASE;

                    for (int mutation_Tracker = 0; mutation_Tracker < 7; mutation_Tracker++)
                    {
                        if (CUDA_sequence_Mutation_tracker[mutation_Tracker][position] != -1)
                        {
                            if (mutation_Tracker == 0)
                            {
                                // overal sequence fitness;
                                int fitness_Point = CUDA_sequence_Mutation_tracker[mutation_Tracker][position];
                                float fitness_Change;

                                fitness_Change = (original_BASE == 0) ? CUDA_A_0_fitness[fitness_Point][new_BASE] : (original_BASE == 1) ? CUDA_T_1_fitness[fitness_Point][new_BASE]
                                                                                                                : (original_BASE == 2)   ? CUDA_G_2_fitness[fitness_Point][new_BASE]
                                                                                                                : (original_BASE == 3)   ? CUDA_C_3_fitness[fitness_Point][new_BASE]
                                                                                                                                         : 1.0f;

                                CUDA_current_gen_Progeny_data[progeny_Index][0] = CUDA_current_gen_Progeny_data[progeny_Index][0] * fitness_Change;
                            }
                            else if (mutation_Tracker == 4)
                            {
                                int proof_Point = CUDA_sequence_Mutation_tracker[mutation_Tracker][position];
                                float proof_Change;

                                proof_Change = (original_BASE == 0) ? CUDA_A_0_probability_Proof_reading[proof_Point][new_BASE] : (original_BASE == 1) ? CUDA_T_1_probability_Proof_reading[proof_Point][new_BASE]
                                                                                                                              : (original_BASE == 2)   ? CUDA_G_2_probability_Proof_reading[proof_Point][new_BASE]
                                                                                                                              : (original_BASE == 3)   ? CUDA_C_3_probability_Proof_reading[proof_Point][new_BASE]
                                                                                                                                                       : 0.00f;

                                CUDA_progeny_Proof_reading_probability[progeny_Index] = CUDA_progeny_Proof_reading_probability[progeny_Index] + proof_Change;

                                if (CUDA_progeny_Proof_reading_probability[progeny_Index] < 0)
                                {
                                    CUDA_progeny_Proof_reading_probability[progeny_Index] = 0;
                                }
                                else if (CUDA_progeny_Proof_reading_probability[progeny_Index] > 1)
                                {
                                    CUDA_progeny_Proof_reading_probability[progeny_Index] = 1;
                                }
                            }
                            else if (mutation_Tracker == 5)
                            {
                                // overall survivability
                                int proof_Point = CUDA_sequence_Mutation_tracker[mutation_Tracker][position];
                                float proof_Change;

                                proof_Change = (original_BASE == 0) ? CUDA_A_0_survivability[proof_Point][new_BASE] : (original_BASE == 1) ? CUDA_T_1_survivability[proof_Point][new_BASE]
                                                                                                                  : (original_BASE == 2)   ? CUDA_G_2_survivability[proof_Point][new_BASE]
                                                                                                                  : (original_BASE == 3)   ? CUDA_C_3_survivability[proof_Point][new_BASE]
                                                                                                                                           : 0.00f;
                                CUDA_progeny_survivability_Probabilities[progeny_Index][0] = CUDA_progeny_survivability_Probabilities[progeny_Index][0] + proof_Change;

                                // if (CUDA_progeny_survivability_Probabilities[progeny_Index][0] < 0)
                                // {
                                //     CUDA_progeny_survivability_Probabilities[progeny_Index][0] = 0;
                                // }
                                // else if (CUDA_progeny_survivability_Probabilities[progeny_Index][0] > 1)
                                // {
                                //     CUDA_progeny_survivability_Probabilities[progeny_Index][0] = 1;
                                // }
                            }
                            else
                            {
                                int shift = (mutation_Tracker == 1) ? 3 : (mutation_Tracker == 2) ? 1
                                                                                                  : 2;

                                int num_of_Recombination_events = CUDA_sequence_Mutation_tracker[mutation_Tracker][position];
                                int start;
                                int stop;

                                if (mutation_Tracker == 6)
                                {
                                    start = CUDA_stride_Array[3];
                                    stop = CUDA_stride_Array[4];
                                }
                                else
                                {

                                    start = CUDA_stride_Array[mutation_Tracker - 1];
                                    stop = CUDA_stride_Array[mutation_Tracker];
                                }

                                int pro_track = 0;

                                for (int rows = start; rows < stop; rows++)
                                {
                                    int hotspot = -1;
                                    float change;

                                    if (original_BASE == 0 && CUDA_A_0_Recombination[rows][1] == (position + 1))
                                    {
                                        hotspot = (int)CUDA_A_0_Recombination[rows][0];
                                        change = CUDA_A_0_Recombination[rows][new_BASE + 2];

                                        // int index_Change = (hotspot * 3) + shift;

                                        // CUDA_current_gen_Progeny_data[progeny_Index][index_Change] = CUDA_current_gen_Progeny_data[progeny_Index][index_Change] + probability_Change;

                                        // pro_track++;
                                    }
                                    else if (original_BASE == 1 && CUDA_T_1_Recombination[rows][1] == (position + 1))
                                    {
                                        hotspot = (int)CUDA_T_1_Recombination[rows][0];
                                        change = CUDA_T_1_Recombination[rows][new_BASE + 2];

                                        // int index_Change = (hotspot * 3) + shift;

                                        // CUDA_current_gen_Progeny_data[progeny_Index][index_Change] = CUDA_current_gen_Progeny_data[progeny_Index][index_Change] + probability_Change;

                                        // pro_track++;
                                    }
                                    else if (original_BASE == 2 && CUDA_G_2_Recombination[rows][1] == (position + 1))
                                    {
                                        hotspot = (int)CUDA_G_2_Recombination[rows][0];
                                        change = CUDA_G_2_Recombination[rows][new_BASE + 2];

                                        // int index_Change = (hotspot * 3) + shift;

                                        // CUDA_current_gen_Progeny_data[progeny_Index][index_Change] = CUDA_current_gen_Progeny_data[progeny_Index][index_Change] + probability_Change;

                                        // pro_track++;
                                    }
                                    else if (original_BASE == 2 && CUDA_C_3_Recombination[rows][1] == (position + 1))
                                    {
                                        hotspot = (int)CUDA_C_3_Recombination[rows][0];
                                        change = CUDA_C_3_Recombination[rows][new_BASE + 2];

                                        // int index_Change = (hotspot * 3) + shift;

                                        // CUDA_current_gen_Progeny_data[progeny_Index][index_Change] = CUDA_current_gen_Progeny_data[progeny_Index][index_Change] + probability_Change;

                                        // pro_track++;
                                    }

                                    if (hotspot != -1)
                                    {
                                        int index_Change = (hotspot * 3) + shift;
                                        if (mutation_Tracker == 1)
                                        {
                                            CUDA_current_gen_Progeny_data[progeny_Index][index_Change] = CUDA_current_gen_Progeny_data[progeny_Index][index_Change] + change;
                                            if (CUDA_current_gen_Progeny_data[progeny_Index][index_Change] < 0)
                                            {
                                                CUDA_current_gen_Progeny_data[progeny_Index][index_Change] = 0;
                                            }
                                            else if (CUDA_current_gen_Progeny_data[progeny_Index][index_Change] > 1)
                                            {
                                                CUDA_current_gen_Progeny_data[progeny_Index][index_Change] = 1;
                                            }
                                        }
                                        else if (mutation_Tracker == 6)
                                        {
                                            CUDA_progeny_survivability_Probabilities[progeny_Index][hotspot + 1] = CUDA_progeny_survivability_Probabilities[progeny_Index][hotspot + 1] + change;
                                            // if (CUDA_progeny_survivability_Probabilities[progeny_Index][hotspot + 1] < 0)
                                            // {
                                            //     CUDA_progeny_survivability_Probabilities[progeny_Index][hotspot + 1] = 0;
                                            // }
                                            // else if (CUDA_progeny_survivability_Probabilities[progeny_Index][hotspot + 1] > 1)
                                            // {
                                            //     CUDA_progeny_survivability_Probabilities[progeny_Index][hotspot + 1] = 1;
                                            // }
                                        }
                                        else
                                        {
                                            CUDA_current_gen_Progeny_data[progeny_Index][index_Change] = CUDA_current_gen_Progeny_data[progeny_Index][index_Change] * change;
                                        }
                                        pro_track++;
                                    }

                                    if (pro_track == num_of_Recombination_events)
                                    {
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        tid += blockDim.x * gridDim.x;
    }
}

__global__ void CUDA_mutate_Progeny(int total_Progeny, int num_mutation_Hotspots, int generation,
                                    float **CUDA_mutation_rates_Hotspot_generation, int **CUDA_mutation_Regions_start_stop,
                                    int **cuda_progeny_Sequences, float *CUDA_progeny_Proof_reading_probability, int proof_reading_Activate,
                                    float **CUDA_A_0_mutation, float **CUDA_T_1_mutation, float **CUDA_G_2_mutation, float **CUDA_C_3_mutation,
                                    int **CUDA_sequence_Mutation_tracker, float **CUDA_current_gen_Progeny_data,
                                    float **CUDA_A_0_fitness, float **CUDA_T_1_fitness, float **CUDA_G_2_fitness, float **CUDA_C_3_fitness,
                                    float **CUDA_A_0_probability_Proof_reading, float **CUDA_T_1_probability_Proof_reading, float **CUDA_G_2_probability_Proof_reading, float **CUDA_C_3_probability_Proof_reading,
                                    float **CUDA_A_0_Recombination, float **CUDA_T_1_Recombination, float **CUDA_G_2_Recombination, float **CUDA_C_3_Recombination,
                                    int *CUDA_stride_Array)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < total_Progeny)
    {
        // int progeny_Index = tid + start_Index;

        for (int mutation_Hotspot = 0; mutation_Hotspot < num_mutation_Hotspots; mutation_Hotspot++)
        {
            float mean = CUDA_mutation_rates_Hotspot_generation[mutation_Hotspot][generation];
            int start = CUDA_mutation_Regions_start_stop[mutation_Hotspot][0] - 1;
            int stop = CUDA_mutation_Regions_start_stop[mutation_Hotspot][1] - 1;

            curandState state;
            curand_init(clock64(), tid, 0, &state);
            int num_of_Mutations = curand_poisson(&state, mean);

            if (proof_reading_Activate == 1)
            {
                float proof_Reading_probability = CUDA_progeny_Proof_reading_probability[tid];

                int corrections = 0;

                for (int i = 0; i < num_of_Mutations; i++)
                {
                    if (curand_uniform(&state) < proof_Reading_probability)
                    {
                        corrections++;
                    }
                }

                num_of_Mutations = num_of_Mutations - corrections;

                // if (num_of_Mutations < 0)
                // {
                //     num_of_Mutations = 0;
                // }
            }

            for (int mutation_Number = 0; mutation_Number < num_of_Mutations; mutation_Number++)
            {
                curand_init(clock64(), tid, 0, &state);
                int position = (int)(curand_uniform(&state) * (stop - start + 1)) + start;

                // now do the mutations
                int original_BASE = cuda_progeny_Sequences[tid][position];
                int new_BASE = -1;

                float rand_num = curand_uniform(&state);
                float cumulative_prob = 0.0f;

                int probable_Base = 0;

                for (int base = 0; base < 4; base++)
                {

                    cumulative_prob += (original_BASE == 0)   ? CUDA_A_0_mutation[mutation_Hotspot][base]
                                       : (original_BASE == 1) ? CUDA_T_1_mutation[mutation_Hotspot][base]
                                       : (original_BASE == 2) ? CUDA_G_2_mutation[mutation_Hotspot][base]
                                       : (original_BASE == 3) ? CUDA_C_3_mutation[mutation_Hotspot][base]
                                                              : 0.0f;

                    // if (original_BASE == 0)
                    // {
                    //     cumulative_prob += CUDA_A_0_mutation[mutation_Hotspot][base];
                    // }
                    // else if (original_BASE == 1)
                    // {
                    //     cumulative_prob += CUDA_T_1_mutation[mutation_Hotspot][base];
                    // }
                    // else if (original_BASE == 2)
                    // {
                    //     cumulative_prob += CUDA_G_2_mutation[mutation_Hotspot][base];
                    // }
                    // else if (original_BASE == 3)
                    // {
                    //     cumulative_prob += CUDA_C_3_mutation[mutation_Hotspot][base];
                    // }

                    if (rand_num < cumulative_prob)
                    {
                        new_BASE = probable_Base;
                        break;
                    }

                    probable_Base++;
                }

                if (original_BASE != new_BASE)
                {
                    // if mutation has occured
                    // change base
                    // see what changes have to be done to sequence profile
                    // dont forget to check proof reading
                    // do fitness and proof reading first

                    // CHANGE to new base

                    cuda_progeny_Sequences[tid][position] = new_BASE;

                    for (int mutation_Tracker = 0; mutation_Tracker < 5; mutation_Tracker++)
                    {
                        if (CUDA_sequence_Mutation_tracker[mutation_Tracker][position] != -1)
                        {
                            if (mutation_Tracker == 0)
                            {
                                // overal sequence fitness;
                                int fitness_Point = CUDA_sequence_Mutation_tracker[mutation_Tracker][position];
                                float fitness_Change;

                                fitness_Change = (original_BASE == 0) ? CUDA_A_0_fitness[fitness_Point][new_BASE] : (original_BASE == 1) ? CUDA_T_1_fitness[fitness_Point][new_BASE]
                                                                                                                : (original_BASE == 2)   ? CUDA_G_2_fitness[fitness_Point][new_BASE]
                                                                                                                : (original_BASE == 3)   ? CUDA_C_3_fitness[fitness_Point][new_BASE]
                                                                                                                                         : 1.0f;

                                CUDA_current_gen_Progeny_data[tid][0] = CUDA_current_gen_Progeny_data[tid][0] * fitness_Change;
                            }
                            else if (mutation_Tracker == 4)
                            {
                                int proof_Point = CUDA_sequence_Mutation_tracker[mutation_Tracker][position];
                                float proof_Change;

                                proof_Change = (original_BASE == 0) ? CUDA_A_0_probability_Proof_reading[proof_Point][new_BASE] : (original_BASE == 1) ? CUDA_T_1_probability_Proof_reading[proof_Point][new_BASE]
                                                                                                                              : (original_BASE == 2)   ? CUDA_G_2_probability_Proof_reading[proof_Point][new_BASE]
                                                                                                                              : (original_BASE == 3)   ? CUDA_C_3_probability_Proof_reading[proof_Point][new_BASE]
                                                                                                                                                       : 0.00f;

                                CUDA_progeny_Proof_reading_probability[tid] = CUDA_progeny_Proof_reading_probability[tid] + proof_Change;
                                if (CUDA_progeny_Proof_reading_probability[tid] < 0)
                                {
                                    CUDA_progeny_Proof_reading_probability[tid] = 0;
                                }
                                else if (CUDA_progeny_Proof_reading_probability[tid] > 1)
                                {
                                    CUDA_progeny_Proof_reading_probability[tid] = 1;
                                }
                            }
                            else
                            {
                                int shift = (mutation_Tracker == 1) ? 3 : (mutation_Tracker == 2) ? 1
                                                                                                  : 2;

                                int num_of_Recombination_events = CUDA_sequence_Mutation_tracker[mutation_Tracker][position];
                                int start = CUDA_stride_Array[mutation_Tracker - 1];
                                int stop = CUDA_stride_Array[mutation_Tracker];

                                int pro_track = 0;

                                for (int rows = start; rows < stop; rows++)
                                {
                                    int hotspot = -1;
                                    float change;

                                    if (original_BASE == 0 && CUDA_A_0_Recombination[rows][1] == (position + 1))
                                    {
                                        hotspot = (int)CUDA_A_0_Recombination[rows][0];
                                        change = CUDA_A_0_Recombination[rows][new_BASE + 2];

                                        // int index_Change = (hotspot * 3) + shift;

                                        // CUDA_current_gen_Progeny_data[progeny_Index][index_Change] = CUDA_current_gen_Progeny_data[progeny_Index][index_Change] + probability_Change;

                                        // pro_track++;
                                    }
                                    else if (original_BASE == 1 && CUDA_T_1_Recombination[rows][1] == (position + 1))
                                    {
                                        hotspot = (int)CUDA_T_1_Recombination[rows][0];
                                        change = CUDA_T_1_Recombination[rows][new_BASE + 2];

                                        // int index_Change = (hotspot * 3) + shift;

                                        // CUDA_current_gen_Progeny_data[progeny_Index][index_Change] = CUDA_current_gen_Progeny_data[progeny_Index][index_Change] + probability_Change;

                                        // pro_track++;
                                    }
                                    else if (original_BASE == 2 && CUDA_G_2_Recombination[rows][1] == (position + 1))
                                    {
                                        hotspot = (int)CUDA_G_2_Recombination[rows][0];
                                        change = CUDA_G_2_Recombination[rows][new_BASE + 2];

                                        // int index_Change = (hotspot * 3) + shift;

                                        // CUDA_current_gen_Progeny_data[progeny_Index][index_Change] = CUDA_current_gen_Progeny_data[progeny_Index][index_Change] + probability_Change;

                                        // pro_track++;
                                    }
                                    else if (original_BASE == 2 && CUDA_C_3_Recombination[rows][1] == (position + 1))
                                    {
                                        hotspot = (int)CUDA_C_3_Recombination[rows][0];
                                        change = CUDA_C_3_Recombination[rows][new_BASE + 2];

                                        // int index_Change = (hotspot * 3) + shift;

                                        // CUDA_current_gen_Progeny_data[progeny_Index][index_Change] = CUDA_current_gen_Progeny_data[progeny_Index][index_Change] + probability_Change;

                                        // pro_track++;
                                    }

                                    if (hotspot != -1)
                                    {
                                        int index_Change = (hotspot * 3) + shift;
                                        if (mutation_Tracker == 1)
                                        {
                                            CUDA_current_gen_Progeny_data[tid][index_Change] = CUDA_current_gen_Progeny_data[tid][index_Change] + change;
                                            if (CUDA_current_gen_Progeny_data[tid][index_Change] < 0)
                                            {
                                                CUDA_current_gen_Progeny_data[tid][index_Change] = 0;
                                            }
                                            else if (CUDA_current_gen_Progeny_data[tid][index_Change] > 1)
                                            {
                                                CUDA_current_gen_Progeny_data[tid][index_Change] = 1;
                                            }
                                        }
                                        else
                                        {
                                            CUDA_current_gen_Progeny_data[tid][index_Change] = CUDA_current_gen_Progeny_data[tid][index_Change] * change;
                                        }
                                        pro_track++;
                                    }

                                    if (pro_track == num_of_Recombination_events)
                                    {
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        tid += blockDim.x * gridDim.x;
    }
}

void functions_library::mutate_Sequences(int **cuda_progeny_Sequences, float **CUDA_current_gen_Progeny_data, float *CUDA_progeny_Proof_reading_probability,
                                         int current_Generation, int num_Mutation_hotspots,
                                         float **CUDA_mutation_rates_Hotspot_generation, int **CUDA_mutation_Regions_start_stop, int **CUDA_sequence_Mutation_tracker,
                                         int total_Progeny, int proof_Reading_Activate)
{
    cout << "Mutating progeny" << endl;

    int total_Cells = total_Progeny;
    int full_Rounds = total_Cells / this->gpu_Limit;
    int partial_Rounds = total_Cells % this->gpu_Limit;

    vector<pair<int, int>> start_stops;

    for (int full = 0; full < full_Rounds; full++)
    {
        int start = full * this->gpu_Limit;
        int stop = start + this->gpu_Limit;
        start_stops.push_back(make_pair(start, stop));
    }

    if (partial_Rounds != 0)
    {
        int start = total_Cells - partial_Rounds;
        start_stops.push_back(make_pair(start, total_Cells));
    }

    for (size_t i = 0; i < start_stops.size(); i++)
    {
        int num_of_values_current = start_stops[i].second - start_stops[i].first;

        // CUDA_mutate_Progeny(int total_Progeny, int start_Index, int num_mutation_Hotspots, int generation,
        //                             float **CUDA_mutation_rates_Hotspot_generation, int **CUDA_mutation_Regions_start_stop,
        //                             int **cuda_progeny_Sequences, float *CUDA_progeny_Proof_reading_probability, int proof_reading_Activate,
        //                             float **CUDA_A_0_mutation, float **CUDA_T_1_mutation, float **CUDA_G_2_mutation, float **CUDA_C_3_mutation,
        //                             int **CUDA_sequence_Mutation_tracker, float **CUDA_current_gen_Progeny_data,
        //                             float **CUDA_A_0_fitness, float **CUDA_T_1_fitness, float **CUDA_G_2_fitness, float **CUDA_C_3_fitness,
        //                             float **CUDA_A_0_probability_Proof_reading, float **CUDA_T_1_probability_Proof_reading, float **CUDA_G_2_probability_Proof_reading, float **CUDA_C_3_probability_Proof_reading,
        //                             float **CUDA_A_0_Recombination, float **CUDA_T_1_Recombination, float **CUDA_G_2_Recombination, float **CUDA_C_3_Recombination,
        //                             int *CUDA_stride_Array)

        // CUDA_mutate_Progeny<<<tot_Blocks, tot_ThreadsperBlock>>>(num_of_values_current, start_stops[i].first, num_Mutation_hotspots, current_Generation,
        //                                                          CUDA_mutation_rates_Hotspot_generation, CUDA_mutation_Regions_start_stop,
        //                                                          cuda_progeny_Sequences, CUDA_progeny_Proof_reading_probability, proof_Reading_Activate,
        //                                                          CUDA_A_0_mutation, CUDA_T_1_mutation, CUDA_G_2_mutation, CUDA_C_3_mutation,
        //                                                          CUDA_sequence_Mutation_tracker, CUDA_current_gen_Progeny_data,
        //                                                          CUDA_A_0_fitness, CUDA_T_1_fitness, CUDA_G_2_fitness, CUDA_C_3_fitness,
        //                                                          CUDA_A_0_probability_Proof_reading, CUDA_T_1_probability_Proof_reading, CUDA_G_2_probability_Proof_reading, CUDA_C_3_probability_Proof_reading,
        //                                                          CUDA_A_0_Recombination, CUDA_T_1_Recombination, CUDA_G_2_Recombination, CUDA_C_3_Recombination,
        //                                                          CUDA_stride_Array);

        cudaError_t err = cudaGetLastError();

        if (err != cudaSuccess)
        {
            printf("CUDA Error: %s\n", cudaGetErrorString(err));

            // Possibly: exit(-1) if program cannot continue....
        }
        cudaDeviceSynchronize();
    }

    cout << "Mutated progeny after " << start_stops.size() << " round(s)." << endl;
}

void functions_library::config_File_crash(string location, string headers)
{
    if (filesystem::exists(location))
    {
        cout << "ERROR FILE: " << location << " ALREADY EXISTS." << endl;
        exit(-1);
    }
    else
    {
        create_File(location, headers);
    }
    cout << endl;
}

void functions_library::create_File(string location, string headers)
{
    cout << "File being created: " << location << endl;
    fstream file;
    file.open(location, ios::out);

    if (file.is_open())
    {
        file << headers << "\n";
        file.close();
    }
}

void functions_library::config_File_delete_create(string location, string headers)
{
    if (filesystem::exists(location))
    {
        filesystem::remove(location);
        cout << "Existing file was deleted: " << location << endl;
    }

    create_File(location, headers);

    cout << endl;
}

void functions_library::create_File(string location)
{
    cout << "File being created: " << location << endl;
    fstream file;
    file.open(location, ios::out);

    if (file.is_open())
    {
        file.close();
    }
}

void functions_library::config_File_delete_create(string location)
{
    if (filesystem::exists(location))
    {
        filesystem::remove(location);
        cout << "Existing file was deleted: " << location << endl;
    }

    create_File(location);

    cout << endl;
}

__global__ void unique_Values(int total, int start_Index, int columns_Total, int **array_2D, int *unique_values_Array)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < total)
    {
        int row_Index = tid + start_Index;
        int start = tid * columns_Total;

        for (int column = 0; column < columns_Total; column++)
        {
            unique_values_Array[start] = array_2D[row_Index][column];
            start++;
        }

        // int cell_Index = tid + (start_Index * columns_Total);
        // int row = cell_Index / columns_Total;
        // int column = cell_Index % columns_Total;

        // int val = array_2D[row][column];

        // unique_values_Array[tid] = val;

        tid += blockDim.x * gridDim.x;
    }
}

__global__ void cuda_fill_parent_Sequences(int genome_Size, int **master_Sequences, char *cuda_reference, int ref_Num)
{

    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < genome_Size)
    {
        // A = 0
        // T = 1
        // G = 2
        // C = 3

        char base = cuda_reference[tid];

        master_Sequences[ref_Num][tid] = (base == '0') ? 0 : (base == '1') ? 1
                                                         : (base == '2')   ? 2
                                                                           : 3;

        // if (base == '0')
        // {
        //     master_Sequence[ref_Num][tid] = 0;
        // }
        // else if (base == '1')
        // {
        //     master_Sequence[ref_Num][tid] = 1;
        // }
        // else if (base == '2')
        // {
        //     master_Sequence[ref_Num][tid] = 2;
        // }
        // else if (base == '3')
        // {
        //     master_Sequence[ref_Num][tid] = 3;
        // }

        tid += blockDim.x * gridDim.x;
    }
}

void functions_library::find_Unique_values(int **progeny_recom_Index_Cuda, int total_Elements, int num_Recom_hotspots,
                                           int start_Index,
                                           string multi_READ, vector<string> &parent_IDs, string &parent_Sequences_Store,
                                           int &genome_Size,
                                           string &progeny_Sequences_Store,
                                           int **CUDA_recombination_hotspots_start_stop,
                                           int mutation_Activate,
                                           float **CUDA_current_gen_Progeny_data, float *CUDA_progeny_Proof_reading_probability, int current_Generation, float **CUDA_mutation_rates_Hotspot_generation, int **CUDA_mutation_Regions_start_stop, int **CUDA_sequence_Mutation_tracker, int proof_Reading_Activate,
                                           vector<string> &progeny_IDs)
{
    cout << "Finding unique parent list" << endl;
    // int total = 5 * 20;
    // int start_Index = 0;
    int columns_Total = num_Recom_hotspots + 1;

    // int **test_Data = create_INT_2D_arrays(5, 20);

    // for (int row = 0; row < 5; row++)
    // {
    //     for (int col = 0; col < 20; col++)
    //     {
    //         if (col == 10)
    //         {
    //             test_Data[row][col] = -1;
    //         }
    //         else
    //         {
    //             test_Data[row][col] = col;
    //         }
    //     }
    // }

    // int **cuda_Test_Array = int_2D_Array_load_to_CUDA(Array_2D, 5, 20);

    // int num_Unique = 0;
    // int *CUDA_num_Unique;
    // cudaMalloc(&CUDA_num_Unique, sizeof(int));
    int *CUDA_unique_values_Array, *unique_values_Array;

    // cudaMemcpy(CUDA_num_Unique, &num_Unique, sizeof(int), cudaMemcpyHostToDevice);

    unique_values_Array = (int *)malloc((total_Elements * columns_Total) * sizeof(int));
    cudaMalloc(&CUDA_unique_values_Array, (total_Elements * columns_Total) * sizeof(int));

    // cout << (total_Elements * columns_Total) << endl;

    unique_Values<<<tot_Blocks, tot_ThreadsperBlock>>>(total_Elements, start_Index, columns_Total, progeny_recom_Index_Cuda, CUDA_unique_values_Array);
    cudaError_t err = cudaGetLastError();

    if (err != cudaSuccess)
    {
        printf("CUDA Error: %s\n", cudaGetErrorString(err));

        // Possibly: exit(-1) if program cannot continue....
    }
    cudaDeviceSynchronize();

    cudaMemcpy(unique_values_Array, CUDA_unique_values_Array, sizeof(int) * (total_Elements * columns_Total), cudaMemcpyDeviceToHost);
    // cudaMemcpy(&num_Unique, CUDA_num_Unique, sizeof(int), cudaMemcpyDeviceToHost);

    cudaFree(CUDA_unique_values_Array);
    // cudaFree(CUDA_num_Unique);

    // for (int i = 0; i < num_Unique; i++)
    // {
    //     printf("%d ", unique_values_Array[i]);
    // }

    int num_per_Core = (total_Elements * columns_Total) / this->CPU_cores;
    int remainder = (total_Elements * columns_Total) % this->CPU_cores;

    vector<thread> threads_vec;

    for (int core_ID = 0; core_ID < this->CPU_cores; core_ID++)
    {
        int start_Cell = core_ID * num_per_Core;
        int stop_Cell = start_Cell + num_per_Core;

        threads_vec.push_back(thread{&functions_library::unique_Collect_Threads, this, unique_values_Array, start_Cell, stop_Cell});
    }

    if (remainder != 0)
    {
        int start_Cell = (total_Elements * columns_Total) - remainder;
        int stop_Cell = (total_Elements * columns_Total);

        threads_vec.push_back(thread{&functions_library::unique_Collect_Threads, this, unique_values_Array, start_Cell, stop_Cell});
    }

    for (thread &t : threads_vec)
    {
        if (t.joinable())
        {
            t.join();
        }
    }

    threads_vec.clear();

    free(unique_values_Array);

    // for (auto it = parent_Indexes_Unique.begin(); it != parent_Indexes_Unique.end(); ++it)
    // {
    //     std::cout << *it << " ";
    // }

    int num_Unique_Parents = parent_Indexes_Unique.size();

    cout << num_Unique_Parents << " unique parents being loaded" << endl;

    // exit(-1);

    int *parent_Indexes = (int *)malloc(num_Unique_Parents * sizeof(int));

    if (multi_READ == "YES")
    {
        cout << "Initiating multi read" << endl;
        int index = 0;
        for (auto it = parent_Indexes_Unique.begin(); it != parent_Indexes_Unique.end(); ++it)
        {
            parent_Indexes[index] = *it;
            sequences.push_back("");
            index++;
        }

        //  int num_per_Core = num_Unique_Parents / this->CPU_cores;
        // int remainder = num_Unique_Parents % this->CPU_cores;

        for (int core_ID = 0; core_ID < this->CPU_cores; core_ID++)
        {
            // int start_Cell = core_ID * num_per_Core;
            // int stop_Cell = start_Cell + num_per_Core;

            //  threads_vec.push_back(thread{&functions_library::read_nFASTA_multi_Read, this, parent_Sequences_Store, parent_IDs, genome_Size, start_Cell, stop_Cell, parent_Indexes});
        }
        if (remainder != 0)
        {
            // int start_Cell = num_Unique_Parents - remainder;
            // int stop_Cell = num_Unique_Parents;

            // threads_vec.push_back(thread{&functions_library::read_nFASTA_multi_Read, this, parent_Sequences_Store, parent_IDs, genome_Size, start_Cell, stop_Cell, parent_Indexes});
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
    else
    {
        cout << "Initiating single read" << endl;
        int index = 0;
        for (auto it = parent_Indexes_Unique.begin(); it != parent_Indexes_Unique.end(); ++it)
        {
            parent_Indexes[index] = *it;
            string name_nFASTA = parent_IDs[parent_Indexes[index]];
            string location = parent_Sequences_Store + "/" + name_nFASTA + ".nfasta";

            sequences.push_back(read_nFASTA(location, genome_Size));

            index++;
        }
    }

    parent_Indexes_Unique.clear();

    // cout << endl;

    // for (size_t i = 0; i < num_Unique_Parents; i++)
    // {
    //     cout << "Index " << i << endl;
    //     cout << sequences[i] << endl;
    //     cout << endl;
    // }

    int **CUDA_all_Parent_Sequences;

    cudaMallocManaged(&CUDA_all_Parent_Sequences, (genome_Size + 1) * num_Unique_Parents * sizeof(int));
    int **tmp = (int **)malloc(num_Unique_Parents * sizeof(tmp[0]));
    for (int i = 0; i < num_Unique_Parents; i++)
    {
        cudaMalloc((void **)&tmp[i], (genome_Size + 1) * sizeof(tmp[0][0]));
    }
    cudaMemcpy(CUDA_all_Parent_Sequences, tmp, num_Unique_Parents * sizeof(int *), cudaMemcpyHostToDevice);

    free(tmp);

    cout << "Loading parents to the GPU" << endl;

    for (int sequence_i = 0; sequence_i < num_Unique_Parents; sequence_i++)
    {
        string sequence = sequences[sequence_i];

        char *reference_full, *cuda_reference;
        reference_full = (char *)malloc((sequence.size() + 1) * sizeof(char));
        cudaMallocManaged(&cuda_reference, (sequence.size() + 1) * sizeof(char));
        strcpy(reference_full, sequence.c_str());
        cudaMemcpy(cuda_reference, reference_full, (sequence.size() + 1) * sizeof(char), cudaMemcpyHostToDevice);

        cuda_fill_parent_Sequences<<<tot_Blocks, tot_ThreadsperBlock>>>(genome_Size, CUDA_all_Parent_Sequences, cuda_reference, sequence_i);

        cudaError_t err = cudaGetLastError();

        if (err != cudaSuccess)
        {
            printf("CUDA Error: %s\n", cudaGetErrorString(err));

            // Possibly: exit(-1) if program cannot continue....
        }
        cudaDeviceSynchronize();

        free(reference_full);
    }

    sequences.clear();
    cout << num_Unique_Parents << " parents loaded to the GPU" << endl;

    // int **test = create_INT_2D_arrays(num_Unique_Parents, genome_Size);
    // cout << "ok" << endl;
    // for (size_t i = 0; i < num_Unique_Parents; i++)
    // {
    //     // cout << i << endl;
    //     cudaMemcpy(test[i], CUDA_all_Parent_Sequences[i], (genome_Size + 1) * sizeof(CUDA_all_Parent_Sequences[0][0]), cudaMemcpyDeviceToHost);
    // }

    // cout << "ok_shuffle" << endl;

    // for (int row = 0; row < num_Unique_Parents; row++)
    // {
    //     for (int column = 0; column < genome_Size; column++)
    //     {
    //         cout << test[row][column];
    //     }
    //     cout << endl
    //          << endl;
    // }

    cout << "Create progeny sequences" << endl;

    int *CUDA_parent_Indexes;
    cudaMallocManaged(&CUDA_parent_Indexes, num_Unique_Parents * sizeof(int));
    cudaMemcpy(CUDA_parent_Indexes, parent_Indexes, num_Unique_Parents * sizeof(int), cudaMemcpyHostToDevice);

    // if (current_Generation == 1)
    // {
    //     int **progeny_recom_Index = create_INT_2D_arrays(total_Elements, 4);

    //     // cout << "ok" << endl;
    //     for (size_t i = 0; i < total_Elements; i++)
    //     {
    //         cudaMemcpy(progeny_recom_Index[i], progeny_recom_Index_Cuda[i], (4) * sizeof(progeny_recom_Index_Cuda[0][0]), cudaMemcpyDeviceToHost);
    //     }

    //     for (int row = 0; row < total_Elements; row++)
    //     {
    //         for (int column = 0; column < 4; column++)
    //         {
    //             cout << progeny_recom_Index[row][column] << " ";
    //         }
    //         cout << endl;
    //     }
    //     for (size_t i = 0; i < num_Unique_Parents; i++)
    //     {
    //         cout << parent_Indexes[i] << endl;
    //     }
    //     // exit(-1);
    // }

    free(parent_Indexes);

    int **cuda_progeny_Sequences;
    cudaMallocManaged(&cuda_progeny_Sequences, (genome_Size + 1) * total_Elements * sizeof(int));
    tmp = (int **)malloc(total_Elements * sizeof(tmp[0]));
    for (int i = 0; i < total_Elements; i++)
    {
        cudaMalloc((void **)&tmp[i], (genome_Size + 1) * sizeof(tmp[0][0]));
    }
    cudaMemcpy(cuda_progeny_Sequences, tmp, total_Elements * sizeof(int *), cudaMemcpyHostToDevice);
    free(tmp);

    // CUDA_progeny_sequence_generation<<<tot_Blocks, tot_ThreadsperBlock>>>(total_Elements, start_Index, CUDA_all_Parent_Sequences, cuda_progeny_Sequences, progeny_recom_Index_Cuda, num_Recom_hotspots, genome_Size, CUDA_recombination_hotspots_start_stop, CUDA_parent_Indexes, num_Unique_Parents);
    err = cudaGetLastError();

    if (err != cudaSuccess)
    {
        printf("CUDA Error: %s\n", cudaGetErrorString(err));

        // Possibly: exit(-1) if program cannot continue....
    }
    cudaDeviceSynchronize();

    cout << "Generated progeny sequences" << endl;

    // if (current_Generation == 1)
    // {
    //     exit(-1);
    // }

    // if (current_Generation == 1)
    // {
    //     int **progeny_Sequences = create_INT_2D_arrays(total_Elements, genome_Size);

    //     // cout << "ok" << endl;
    //     for (size_t i = 0; i < total_Elements; i++)
    //     {
    //         cudaMemcpy(progeny_Sequences[i], cuda_progeny_Sequences[i], (genome_Size + 1) * sizeof(cuda_progeny_Sequences[0][0]), cudaMemcpyDeviceToHost);
    //     }

    //     for (int row = 0; row < total_Elements; row++)
    //     {
    //         for (int column = 0; column < genome_Size; column++)
    //         {
    //             cout << progeny_Sequences[row][column];
    //         }
    //         cout << endl
    //              << endl;
    //     }

    //     exit(-1);
    // }

    cudaFree(CUDA_all_Parent_Sequences);
    cudaFree(CUDA_parent_Indexes);

    if (mutation_Activate != 0)
    {
        cout << "Mutating progeny sequences" << endl;
        // CUDA_mutate_Progeny<<<tot_Blocks, tot_ThreadsperBlock>>>(total_Elements, start_Index, mutation_Activate, current_Generation,
        //                                                          CUDA_mutation_rates_Hotspot_generation, CUDA_mutation_Regions_start_stop,
        //                                                          cuda_progeny_Sequences, CUDA_progeny_Proof_reading_probability, proof_Reading_Activate,
        //                                                          CUDA_A_0_mutation, CUDA_T_1_mutation, CUDA_G_2_mutation, CUDA_C_3_mutation,
        //                                                          CUDA_sequence_Mutation_tracker, CUDA_current_gen_Progeny_data,
        //                                                          CUDA_A_0_fitness, CUDA_T_1_fitness, CUDA_G_2_fitness, CUDA_C_3_fitness,
        //                                                          CUDA_A_0_probability_Proof_reading, CUDA_T_1_probability_Proof_reading, CUDA_G_2_probability_Proof_reading, CUDA_C_3_probability_Proof_reading,
        //                                                          CUDA_A_0_Recombination, CUDA_T_1_Recombination, CUDA_G_2_Recombination, CUDA_C_3_Recombination,
        //                                                          CUDA_stride_Array);

        err = cudaGetLastError();

        if (err != cudaSuccess)
        {
            printf("CUDA Error: %s\n", cudaGetErrorString(err));
            exit(-1);

            // Possibly: exit(-1) if program cannot continue....
        }
        cudaDeviceSynchronize();

        cout << "Progeny sequences mutated" << endl;
    }

    cout << "Writing progeny sequences" << endl;

    int **progeny_Sequences = create_INT_2D_arrays(total_Elements, genome_Size);

    // cout << "ok" << endl;
    for (size_t i = 0; i < total_Elements; i++)
    {
        progeny_IDs.push_back(to_string(current_Generation + 1) + "_" + to_string(start_Index + i));
        // cout << i << endl;
        cudaMemcpy(progeny_Sequences[i], cuda_progeny_Sequences[i], (genome_Size + 1) * sizeof(cuda_progeny_Sequences[0][0]), cudaMemcpyDeviceToHost);
    }

    cudaFree(cuda_progeny_Sequences);

    if (multi_READ == "YES")
    {
        cout << "Initiating multi write" << endl;

        // int num_per_Core = total_Elements / this->CPU_cores;
        int remainder = total_Elements % this->CPU_cores;

        for (int core_ID = 0; core_ID < this->CPU_cores; core_ID++)
        {
            // int start_Cell = core_ID * num_per_Core;
            // int stop_Cell = start_Cell + num_per_Core;

            //  threads_vec.push_back(thread{&functions_library::write_nFASTA_multi_WRITE, this, progeny_Sequences_Store, progeny_IDs, genome_Size, start_Cell, stop_Cell, progeny_Sequences});
        }

        if (remainder != 0)
        {
            // int start_Cell = total_Elements - remainder;
            // int stop_Cell = total_Elements;

            // threads_vec.push_back(thread{&functions_library::write_nFASTA_multi_WRITE, this, progeny_Sequences_Store, progeny_IDs, genome_Size, start_Cell, stop_Cell, progeny_Sequences});
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
    else
    {
        cout << "Initiating single write" << endl;
        for (int progeny_Num = 0; progeny_Num < total_Elements; progeny_Num++)
        {
            string file_Location = progeny_Sequences_Store + "/" + progeny_IDs[progeny_Num] + ".nfasta";
            config_File_progeny(file_Location);

            fstream nFASTA;
            nFASTA.open(file_Location, ios::app);

            if (nFASTA.is_open())
            {
                for (int base = 0; base < genome_Size; base++)
                {
                    nFASTA << progeny_Sequences[progeny_Num][base];
                }
                nFASTA.close();
            }
        }
    }

    free(progeny_Sequences);

    cout << "Progeny sequences were written: " << progeny_Sequences_Store << endl;

    // for (int row = 0; row < total_Elements; row++)
    // {
    //     for (int column = 0; column < genome_Size; column++)
    //     {
    //         cout << progeny_Sequences[row][column];
    //     }
    //     cout << endl
    //          << endl;
    // }
}

void functions_library::config_File_progeny(string location)
{
    if (filesystem::exists(location))
    {
        filesystem::remove(location);
        // cout << "Existing file was deleted: " << location << endl;
    }

    fstream file;
    file.open(location, ios::out);

    if (file.is_open())
    {
        file.close();
    }

    // cout << endl;
}

void functions_library::write_nFASTA_multi_WRITE(string progeny_Sequences_Store, int sum_Progeny_in_Generation, int genome_Size, int start, int stop, int **progeny_Sequences, int generation_Current)
{
    for (int i = start; i < stop; i++)
    {
        string progeny_Name = to_string(generation_Current + 1) + "_" + to_string(i + sum_Progeny_in_Generation);
        string file_Location = progeny_Sequences_Store + "/" + progeny_Name + ".nfasta";

        if (filesystem::exists(file_Location))
        {
            filesystem::remove(file_Location);
            cout << "Existing file was deleted: " << file_Location << endl;
        }

        fstream nFASTA;
        nFASTA.open(file_Location, ios::app);

        if (nFASTA.is_open())
        {
            for (int base = 0; base < genome_Size; base++)
            {
                nFASTA << progeny_Sequences[i][base];
            }
            nFASTA.close();
        }
        // cout << "Wrote: " << progeny_Name << endl;
    }
}

void functions_library::read_nFASTA_multi_Read(string parent_Sequences_Store, int genome_Size, int start, int stop, int *parent_Indexes, int generation_Current)
{
    vector<string> sequences_Read;

    for (int i = start; i < stop; i++)
    {
        string name_nFASTA = to_string(generation_Current) + "_" + to_string(parent_Indexes[i]);
        string file_location = parent_Sequences_Store + "/" + name_nFASTA + ".nfasta";

        // cout << "Processing parent file: " << file_location << endl;

        string reference_Genome = "";

        fstream reference_File;
        reference_File.open(file_location, ios::in);

        if (reference_File.is_open())
        {
            string line;

            while (getline(reference_File, line))
            {
                reference_Genome.append(line);
            }

            reference_File.close();
        }

        reference_Genome = reference_Genome.substr(0, genome_Size);

        sequences_Read.push_back(reference_Genome);
    }

    unique_lock<shared_mutex> ul(g_mutex);
    int index = 0;
    for (int i = start; i < stop; i++)
    {
        sequences[i] = sequences_Read[index];
        index++;
    }
}

void functions_library::read_nFASTA_multi_Read_CELLS(string parent_Sequences_Store, int genome_Size, int start, int stop, int **parent_Indexes_parents_and_their_Cells, int generation_Current)
{
    vector<string> sequences_Read;

    for (int i = start; i < stop; i++)
    {
        string name_nFASTA = to_string(generation_Current) + "_" + to_string(parent_Indexes_parents_and_their_Cells[0][i]);
        string file_location = parent_Sequences_Store + "/" + name_nFASTA + ".nfasta";

        // cout << "Processing parent file: " << file_location << endl;

        string reference_Genome = "";

        fstream reference_File;
        reference_File.open(file_location, ios::in);

        if (reference_File.is_open())
        {
            string line;

            while (getline(reference_File, line))
            {
                reference_Genome.append(line);
            }

            reference_File.close();
        }

        reference_Genome = reference_Genome.substr(0, genome_Size);

        sequences_Read.push_back(reference_Genome);
    }

    unique_lock<shared_mutex> ul(g_mutex);
    int index = 0;
    for (int i = start; i < stop; i++)
    {
        sequences[i] = sequences_Read[index];
        index++;
    }
}

string functions_library::read_nFASTA(string file_location, int genome_Size)
{
    string reference_Genome = "";

    string header;

    // cout << "Processing parent file: " << file_location << endl;

    fstream reference_File;
    reference_File.open(file_location, ios::in);

    if (reference_File.is_open())
    {
        string line;

        // read header line;
        // getline(reference_File, header);
        // cout << "Reference genome header: " << header.substr(1) << endl;

        while (getline(reference_File, line))
        {
            reference_Genome.append(line);
        }

        reference_File.close();
    }

    // transform(reference_Genome.begin(), reference_Genome.end(), reference_Genome.begin(), ::toupper);
    // genome_Size = reference_Genome.length();

    // cout << "Reference genome size (bp): " << genome_Size << "\n\n";

    reference_Genome = reference_Genome.substr(0, genome_Size);

    return reference_Genome;
}

void functions_library::unique_Collect_Threads(int *unique_values_Array, int start, int stop)
{
    set<int> parent_Indexes;

    for (int i = start; i < stop; i++)
    {
        if (unique_values_Array[i] != -1)
        {
            parent_Indexes.insert(unique_values_Array[i]);
        }
    }

    unique_lock<shared_mutex> ul(g_mutex);

    for (auto it = parent_Indexes.begin(); it != parent_Indexes.end(); ++it)
    {
        parent_Indexes_Unique.insert(*it);
    }
}

void functions_library::hard_Load_progeny(int generation, string &parent_Sequences_Store, int genome_Size,
                                          vector<string> &parent_IDs,
                                          int total_Progeny,
                                          int **progeny_recom_Index_Cuda, int num_Hotspots, int **CUDA_recombination_hotspots_start_stop,
                                          string multi_READ,
                                          string &progeny_Sequences_Store,
                                          int mutation_Activate,
                                          float **CUDA_current_gen_Progeny_data, float *CUDA_progeny_Proof_reading_probability, float **CUDA_mutation_rates_Hotspot_generation, int **CUDA_mutation_Regions_start_stop, int **CUDA_sequence_Mutation_tracker, int proof_Reading_Activate,
                                          vector<string> &progeny_IDs)
{

    cout << "Creating progeny sequences" << endl;

    if (num_Hotspots == -1)
    {
        num_Hotspots = 0;
    }
    if (mutation_Activate == -1)
    {
        mutation_Activate = 0;
    }

    int total_Cells = total_Progeny;
    int full_Rounds = total_Cells / this->gpu_Limit;
    int partial_Rounds = total_Cells % this->gpu_Limit;

    vector<pair<int, int>> start_stops;

    for (int full = 0; full < full_Rounds; full++)
    {
        int start = full * this->gpu_Limit;
        int stop = start + this->gpu_Limit;
        start_stops.push_back(make_pair(start, stop));
    }

    if (partial_Rounds != 0)
    {
        int start = total_Cells - partial_Rounds;
        start_stops.push_back(make_pair(start, total_Cells));
    }

    for (size_t i = 0; i < start_stops.size(); i++)
    {
        cout << "\nRound " << i + 1 << " of " << start_stops.size() << endl;
        int num_of_values_current = start_stops[i].second - start_stops[i].first;

        // find_Unique_values(int **CUDA_Array_2D, int total_Elements, int num_Recom_hotspots,
        //                                            int start_Index,
        //                                            string multi_READ, vector<string> &parent_IDs, string &parent_Sequences_Store,
        //                                            int &genome_Size)

        find_Unique_values(progeny_recom_Index_Cuda, num_of_values_current, num_Hotspots,
                           start_stops[i].first,
                           multi_READ, parent_IDs, parent_Sequences_Store,
                           genome_Size,
                           progeny_Sequences_Store,
                           CUDA_recombination_hotspots_start_stop,
                           mutation_Activate,
                           CUDA_current_gen_Progeny_data, CUDA_progeny_Proof_reading_probability, generation, CUDA_mutation_rates_Hotspot_generation, CUDA_mutation_Regions_start_stop, CUDA_sequence_Mutation_tracker, proof_Reading_Activate,
                           progeny_IDs);
        cout << endl;
    }

    // cout << "\nCreated progeny sequences after " << start_stops.size() << " rounds." << endl;
}

void functions_library::clear_Array_INT(int **CUDA_2D_array, int rows)
{
    cudaSetDevice(CUDA_device_IDs[0]);
    for (int i = 0; i < rows; i++)
    {
        cudaFree(CUDA_2D_array[i]);
    }
    cudaFree(CUDA_2D_array);
}

void functions_library::clear_Array_FLOAT(float **CUDA_2D_array, int rows)
{
    cudaSetDevice(CUDA_device_IDs[0]);
    for (int i = 0; i < rows; i++)
    {
        cudaFree(CUDA_2D_array[i]);
    }
    cudaFree(CUDA_2D_array);
}

void functions_library::clear_Array_float_CPU(float **array_cpu, int rows)
{
    for (int i = 0; i < rows; i++)
    {
        free(array_cpu[i]);
    }
    free(array_cpu);
}

void functions_library::clear_Array_int_CPU(int **array_cpu, int rows)
{
    for (int i = 0; i < rows; i++)
    {
        free(array_cpu[i]);
    }
    free(array_cpu);
}

void functions_library::read_Profiles_multi_Thread_CELL(int start, int stop, int **parent_Indexes_parents_and_their_Cells,
                                                        float **current_gen_Parent_data, float *parent_Proof_reading_probability,
                                                        int generation_Current, string parent_Profiles_Store,
                                                        float **parent_survivability_Probabilities)
{
    for (int parent = start; parent < stop; parent++)
    {
        string profile_Name = parent_Profiles_Store + "/" + to_string(generation_Current) + "_" + to_string(parent_Indexes_parents_and_their_Cells[0][parent]) + ".profile";
        fstream parent_Profile;
        parent_Profile.open(profile_Name, ios::in);

        string parent_Profile_Line;
        getline(parent_Profile, parent_Profile_Line);
        vector<string> line_Data;
        split(line_Data, parent_Profile_Line, '\t');

        for (int i = 0; i < line_Data.size(); i++)
        {
            current_gen_Parent_data[parent][i] = stof(line_Data[i]);
        }
        parent_Profile.close();

        if (proof_reading_Activate_parent != 0)
        {
            fstream parent_Proof_Reading_prob;
            parent_Proof_Reading_prob.open(profile_Name + "_prob", ios::in);
            string parent_Proof_Prob_Line;
            getline(parent_Proof_Reading_prob, parent_Proof_Prob_Line);
            parent_Proof_reading_probability[parent] = stof(parent_Proof_Prob_Line);
            parent_Proof_Reading_prob.close();
        }

        fstream parent_Survivability_prob;
        parent_Survivability_prob.open(profile_Name + "_surv", ios::in);
        string parent_survivability_Line;
        getline(parent_Survivability_prob, parent_survivability_Line);
        split(line_Data, parent_survivability_Line, '\t');
        for (int i = 0; i < line_Data.size(); i++)
        {
            parent_survivability_Probabilities[parent][i] = stof(line_Data[i]);
        }
        parent_Survivability_prob.close();
    }
}

void functions_library::process_Cells(string &multi_READ, int &generation_Current, int &sum_Progeny_in_Generation,
                                      int &num_of_Cells, int &start, int &stop,
                                      vector<int> &cells_Start_Stop,
                                      string &parent_Profiles_Store, int *parents,
                                      string &parent_Sequences_Store,
                                      string &progeny_File, string &progeny_Recombination_File, string &sequence_Profiles,
                                      string &progeny_Sequences_Store, string progeny_Profile_Store,
                                      int &processed_Cells, vector<int> &surviving_Progeny)
{

    vector<thread> threads_vec;

    if (this->mode != "PARENT")
    {
        cout << "\nPer cell mode" << endl;

        int *per_Cell_parents_Stride = (int *)malloc((num_of_Cells + 1) * sizeof(int));
        per_Cell_parents_Stride[0] = 0;

        // int count_Parents_test = 0;

        int cell_index = 0;
        for (int cell = start; cell < stop; cell++)
        {
            per_Cell_parents_Stride[cell_index + 1] = per_Cell_parents_Stride[cell_index] + (cells_Start_Stop[cell + 1] - cells_Start_Stop[cell]);
            // count_Parents_test = count_Parents_test + (cells_Start_Stop[cell + 1] - cells_Start_Stop[cell]);
            cell_index++;
        }

        int count_Parents = per_Cell_parents_Stride[num_of_Cells];
        cout << "Number of parents to be simulated: " << count_Parents << endl;

        // cout << "Test " << count_Parents_test << endl;

        // row 0 = indexes
        // row 1 = cells
        int **parent_Indexes_parents_and_their_Cells = create_INT_2D_arrays(2, count_Parents);
        // int *parent_Indexes = (int *)malloc(count_Parents * sizeof(int));
        // int *parents_and_their_Cells = (int *)malloc(count_Parents * sizeof(int));
        cell_index = 0;
        int cell_ID = 0;

        fstream parent_Cells;
        parent_Cells.open(this->cells_of_parents, ios::app);

        for (int cell = start; cell < stop; cell++)
        {
            for (int index_Start = cells_Start_Stop[cell]; index_Start < cells_Start_Stop[cell + 1]; index_Start++)
            {
                parent_Indexes_parents_and_their_Cells[0][cell_index] = parents[index_Start];
                parent_Indexes_parents_and_their_Cells[1][cell_index] = cell_ID;
                parent_Cells << to_string(generation_Current) << "_" << to_string(parents[index_Start]) << "\t" << to_string(generation_Current) << "_" + to_string(processed_Cells + cell_ID) << "\n";
                cell_index++;
            }
            cell_ID++;
        }

        parent_Cells.close();

        // for (size_t i = 0; i < 2; i++)
        // {
        //     for (size_t col = 0; col < count_Parents; col++)
        //     {
        //         cout << parent_Indexes_parents_and_their_Cells[i][col] << " ";
        //     }
        //     cout << endl;
        // }
        // exit(-1);
        // cout << endl;
        // for (size_t i = 0; i < num_of_Cells; i++)
        // {
        //     for (size_t col = per_Cell_parents_Stride[i]; col < per_Cell_parents_Stride[i + 1]; col++)
        //     {
        //         cout << parent_Indexes_parents_and_their_Cells[0][col] << " ";
        //     }
        //     cout << endl;
        // }

        float **current_gen_Parent_data = create_FLOAT_2D_arrays(count_Parents, 1 + (3 * recombination_hotspots));
        float *parent_Proof_reading_probability;
        if (proof_reading_Activate_parent != 0)
        {
            parent_Proof_reading_probability = (float *)malloc(sizeof(float) * count_Parents);
        }

        float **parent_survivability_Probabilities = create_FLOAT_2D_arrays(count_Parents, 1 + recombination_hotspots);

        if (multi_READ == "YES")
        {
            cout << "Initiating multi read" << endl;
            int num_per_Core = count_Parents / this->CPU_cores;
            int remainder = count_Parents % this->CPU_cores;

            for (int core_ID = 0; core_ID < this->CPU_cores; core_ID++)
            {
                int start_Cell = core_ID * num_per_Core;
                int stop_Cell = start_Cell + num_per_Core;

                threads_vec.push_back(thread{&functions_library::read_Profiles_multi_Thread_CELL, this, start_Cell, stop_Cell, parent_Indexes_parents_and_their_Cells, current_gen_Parent_data, parent_Proof_reading_probability, generation_Current, parent_Profiles_Store, parent_survivability_Probabilities});
            }
            if (remainder != 0)
            {
                int start_Cell = count_Parents - remainder;
                int stop_Cell = count_Parents;

                threads_vec.push_back(thread{&functions_library::read_Profiles_multi_Thread_CELL, this, start_Cell, stop_Cell, parent_Indexes_parents_and_their_Cells, current_gen_Parent_data, parent_Proof_reading_probability, generation_Current, parent_Profiles_Store, parent_survivability_Probabilities});
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
        else
        {
            cout << "Initiating single read" << endl;

            for (int parent = 0; parent < count_Parents; parent++)
            {
                string profile_Name = parent_Profiles_Store + "/" + to_string(generation_Current) + "_" + to_string(parent_Indexes_parents_and_their_Cells[0][parent]) + ".profile";
                fstream parent_Profile;
                parent_Profile.open(profile_Name, ios::in);

                string parent_Profile_Line;
                getline(parent_Profile, parent_Profile_Line);
                vector<string> line_Data;
                split(line_Data, parent_Profile_Line, '\t');

                for (int i = 0; i < line_Data.size(); i++)
                {
                    current_gen_Parent_data[parent][i] = stof(line_Data[i]);
                }
                parent_Profile.close();

                if (proof_reading_Activate_parent != 0)
                {
                    fstream parent_Proof_Reading_prob;
                    parent_Proof_Reading_prob.open(profile_Name + "_prob", ios::in);
                    string parent_Proof_Prob_Line;
                    getline(parent_Proof_Reading_prob, parent_Proof_Prob_Line);
                    parent_Proof_reading_probability[parent] = stof(parent_Proof_Prob_Line);
                    parent_Proof_Reading_prob.close();
                }

                fstream parent_Survivability_prob;
                parent_Survivability_prob.open(profile_Name + "_surv", ios::in);
                string parent_survivability_Line;
                getline(parent_Survivability_prob, parent_survivability_Line);
                split(line_Data, parent_survivability_Line, '\t');
                for (int i = 0; i < line_Data.size(); i++)
                {
                    parent_survivability_Probabilities[parent][i] = stof(line_Data[i]);
                }
                parent_Survivability_prob.close();
            }
        }

        // for (size_t i = 0; i < count_Parents; i++)
        // {
        //     for (size_t col = 0; col < (1 + recombination_hotspots); col++)
        //     {
        //         cout << parent_survivability_Probabilities[i][col] << " ";
        //     }
        //     cout << endl;
        //     // cout << " proof: " << parent_Proof_reading_probability[i] << endl;
        // }

        // exit(-1);

        //! clear gpu = DONE
        float **CUDA_current_gen_Parent_data = float_2D_Array_load_to_CUDA(current_gen_Parent_data, count_Parents, 1 + (3 * recombination_hotspots));
        clear_Array_float_CPU(current_gen_Parent_data, count_Parents);

        cout << "Parents loaded to GPU" << endl;

        cout << "Determining progeny numbers" << endl;
        //! clear gpu = done
        int **cuda_Progeny_numbers = create_CUDA_2D_int(recombination_hotspots + 1, count_Parents);

        int full_Rounds = count_Parents / this->gpu_Limit;
        int partial_Rounds = count_Parents % this->gpu_Limit;

        vector<pair<int, int>> start_stops;

        for (int full = 0; full < full_Rounds; full++)
        {
            int start = full * this->gpu_Limit;
            int stop = start + this->gpu_Limit;
            start_stops.push_back(make_pair(start, stop));
        }

        if (partial_Rounds != 0)
        {
            int start = count_Parents - partial_Rounds;
            start_stops.push_back(make_pair(start, count_Parents));
        }

        for (size_t i = 0; i < start_stops.size(); i++)
        {
            int num_of_values_current = start_stops[i].second - start_stops[i].first;

            curandState *state;
            cudaMalloc((void **)&state, num_of_values_current * sizeof(curandState));

            if (progeny_distribution_Type == "Gamma")
            {
                //CUDA_gamma_distribution_PROGENY<<<tot_Blocks, tot_ThreadsperBlock>>>(state, num_of_values_current, cuda_Progeny_numbers, progeny_shape, progeny_scale, start_stops[i].first, recombination_hotspots, CUDA_current_gen_Parent_data);
            }
            else if (progeny_distribution_Type == "Negative binomial")
            {
            }
            cudaError_t err = cudaGetLastError();
            if (err != cudaSuccess)
            {
                printf("CUDA Error: %s\n", cudaGetErrorString(err));

                // Possibly: exit(-1) if program cannot continue....
            }
            cudaDeviceSynchronize();
            // break;

            cudaFree(state);
        }

        cout << "Completed progeny generation via " << start_stops.size() << " GPU rounds" << endl;
        start_stops.clear();

        // for (size_t col = 0; col < count_Parents; col++)
        // {
        //     for (size_t row = 0; row < (recombination_hotspots + 1); row++)
        //     {
        //         cout << Progeny_numbers[row][col] << " ";
        //     }
        //     cout << endl;
        // }

        int sum_Progeny = sum_CUDA(cuda_Progeny_numbers, count_Parents);
        cout << "Sum progeny in cells: " << sum_Progeny << endl;
        // sum_Progeny_in_Generation = sum_Progeny_in_Generation + sum_Progeny;
        //  cout << sum_Progeny_in_Generation << endl;

        // create fitness distribution;
        //! clear gpu = done
        float **CUDA_fitness_distribution;

        //! clear gpu = done
        int *CUDA_per_Cell_parents_Stride;
        cudaMallocManaged(&CUDA_per_Cell_parents_Stride, (num_of_Cells + 1) * sizeof(int));
        cudaMemcpy(CUDA_per_Cell_parents_Stride, per_Cell_parents_Stride, (num_of_Cells + 1) * sizeof(int), cudaMemcpyHostToDevice);
        //! clear gpu = done
        int **CUDA_parent_Indexes_parents_and_their_Cells = int_2D_Array_load_to_CUDA(parent_Indexes_parents_and_their_Cells, 2, count_Parents);

        if (recombination_hotspots != 0)
        {
            cout << "Generating selectivity distributions for putative parents" << endl;
            //! clear gpu = DONE;
            float **CUDA_hotspot_selectivity_Summations = create_CUDA_2D_FLOAT(num_of_Cells, recombination_hotspots);

            start_stops.clear();

            full_Rounds = num_of_Cells / this->gpu_Limit;
            partial_Rounds = num_of_Cells % this->gpu_Limit;

            for (int full = 0; full < full_Rounds; full++)
            {
                int start = full * this->gpu_Limit;
                int stop = start + this->gpu_Limit;
                start_stops.push_back(make_pair(start, stop));
            }

            if (partial_Rounds != 0)
            {
                int start = num_of_Cells - partial_Rounds;
                start_stops.push_back(make_pair(start, num_of_Cells));
            }

            for (size_t i = 0; i < start_stops.size(); i++)
            {
                int num_of_values_current = start_stops[i].second - start_stops[i].first;
                // CUDA_summation_Selectivity_CELLS(int cells, int start_Index, int num_Hotspots, float **CUDA_current_gen_Parent_data, float **CUDA_hotspot_selectivity_Summations, int *CUDA_per_Cell_parents_Stride)
                CUDA_summation_Selectivity_CELLS<<<tot_Blocks, tot_ThreadsperBlock>>>(num_of_values_current, start_stops[i].first, recombination_hotspots, CUDA_current_gen_Parent_data, CUDA_hotspot_selectivity_Summations, CUDA_per_Cell_parents_Stride);

                cudaError_t err = cudaGetLastError();
                if (err != cudaSuccess)
                {
                    printf("CUDA Error: %s\n", cudaGetErrorString(err));

                    // Possibly: exit(-1) if program cannot continue....
                }
                cudaDeviceSynchronize();
            }

            start_stops.clear();

            // float **test_2d = load_to_Host_FLOAT(CUDA_hotspot_selectivity_Summations, num_of_Cells, recombination_hotspots);

            // for (size_t i = 0; i < num_of_Cells; i++)
            // {
            //     for (size_t col = 0; col < recombination_hotspots; col++)
            //     {
            //         cout << test_2d[i][col] << " ";
            //     }
            //     cout << endl;
            // }

            // clear_Array_float_CPU(test_2d, num_of_Cells);
            CUDA_fitness_distribution = create_CUDA_2D_FLOAT(count_Parents, recombination_hotspots);

            cout << "Hotspots' summations done" << endl;

            int total_Cells = count_Parents;
            full_Rounds = total_Cells / this->gpu_Limit;
            partial_Rounds = total_Cells % this->gpu_Limit;

            for (int full = 0; full < full_Rounds; full++)
            {
                int start = full * this->gpu_Limit;
                int stop = start + this->gpu_Limit;
                start_stops.push_back(make_pair(start, stop));
            }

            if (partial_Rounds != 0)
            {
                int start = total_Cells - partial_Rounds;
                start_stops.push_back(make_pair(start, total_Cells));
            }

            for (size_t i = 0; i < start_stops.size(); i++)
            {
                int num_of_values_current = start_stops[i].second - start_stops[i].first;

                CUDA_Selectivity_Distribution_CELLS<<<tot_Blocks, tot_ThreadsperBlock>>>(recombination_hotspots, num_of_values_current, CUDA_current_gen_Parent_data, CUDA_hotspot_selectivity_Summations, CUDA_fitness_distribution, start_stops[i].first, CUDA_parent_Indexes_parents_and_their_Cells);
                cudaError_t err = cudaGetLastError();

                if (err != cudaSuccess)
                {
                    printf("CUDA Error: %s\n", cudaGetErrorString(err));

                    // Possibly: exit(-1) if program cannot continue....
                }
                cudaDeviceSynchronize();
            }

            cout << "Generated selectivity distributions after " << start_stops.size() << " rounds" << endl;
            start_stops.clear();

            // float **test_Distributions = load_to_Host_FLOAT(CUDA_fitness_distribution, count_Parents, recombination_hotspots);

            // for (size_t i = 0; i < count_Parents; i++)
            // {
            //     for (size_t col = 0; col < recombination_hotspots; col++)
            //     {
            //         cout << test_Distributions[i][col] << " ";
            //     }
            //     cout << endl;
            // }

            clear_Array_FLOAT(CUDA_hotspot_selectivity_Summations, num_of_Cells);
        }

        //! create progeny recombination index first.

        cout << "Creating progeny indexes for parents" << endl;

        //! clear properly = done
        int **Progeny_numbers = load_to_Host(cuda_Progeny_numbers, (recombination_hotspots + 1), count_Parents);
        // for (size_t i = 0; i < count_Parents; i++)
        // {
        //     for (size_t row = 0; row < (recombination_hotspots + 1); row++)
        //     {
        //         cout << Progeny_numbers[row][i] << " ";
        //     }
        //     cout << endl;
        // }

        int *progeny_Stride_Index = (int *)malloc((count_Parents + 1) * sizeof(int));
        progeny_Stride_Index[0] = 0;
        //! clear gpu = DONE
        int **progeny_recom_Index_Cuda = create_CUDA_2D_int(sum_Progeny, recombination_hotspots + 1);

        for (int parent = 0; parent < count_Parents; parent++)
        {
            int total_Cells = Progeny_numbers[0][parent];

            full_Rounds = total_Cells / this->gpu_Limit;
            partial_Rounds = total_Cells % this->gpu_Limit;

            start_stops.clear();

            for (int full = 0; full < full_Rounds; full++)
            {
                int start = full * this->gpu_Limit;
                int stop = start + this->gpu_Limit;
                start_stops.push_back(make_pair(start, stop));
            }

            if (partial_Rounds != 0)
            {
                int start = total_Cells - partial_Rounds;
                start_stops.push_back(make_pair(start, total_Cells));
            }

            // float fill_Value = -1;

            for (size_t i = 0; i < start_stops.size(); i++)
            {
                int num_of_values_current = start_stops[i].second - start_stops[i].first;

                progeny_Array_CUDA_CELLS<<<tot_Blocks, tot_ThreadsperBlock>>>(num_of_values_current, progeny_Stride_Index[parent], start_stops[i].first, parent, recombination_hotspots, cuda_Progeny_numbers, progeny_recom_Index_Cuda, -1, CUDA_fitness_distribution, CUDA_per_Cell_parents_Stride, CUDA_parent_Indexes_parents_and_their_Cells);
                cudaError_t err = cudaGetLastError();

                if (err != cudaSuccess)
                {
                    printf("CUDA Error 1: %s\n", cudaGetErrorString(err));
                    exit(-1);
                    // Possibly: exit(-1) if program cannot continue....
                }
                cudaDeviceSynchronize();
            }

            cout << "Created progeny index for parent " << parent + 1 << " after " << start_stops.size() << " round(s)." << endl;

            progeny_Stride_Index[parent + 1] = progeny_Stride_Index[parent] + Progeny_numbers[0][parent];
        }

        cudaFree(CUDA_per_Cell_parents_Stride);
        clear_Array_INT(CUDA_parent_Indexes_parents_and_their_Cells, 2);
        free(per_Cell_parents_Stride);

        // int **test_Load = load_to_Host(progeny_recom_Index_Cuda, sum_Progeny, recombination_hotspots + 1);

        // for (size_t row = 0; row < sum_Progeny; row++)
        // {
        //     for (size_t col = 0; col < (recombination_hotspots + 1); col++)
        //     {
        //         cout << test_Load[row][col] << " ";
        //     }
        //     cout << endl;
        // }
        // clear_Array_int_CPU(test_Load, sum_Progeny);

        // cout << "check" << endl;
        // for (size_t i = 0; i < count_Parents + 1; i++)
        // {
        //     cout << progeny_Stride_Index[i] << endl;
        // }

        // exit(-1);

        start_stops.clear();
        if (recombination_hotspots != 0)
        {
            clear_Array_FLOAT(CUDA_fitness_distribution, count_Parents);

            //! cleared = done
            int *CUDA_progeny_Stride_Index;
            cudaMallocManaged(&CUDA_progeny_Stride_Index, (count_Parents + 1) * sizeof(int));
            cudaMemcpy(CUDA_progeny_Stride_Index, progeny_Stride_Index, (count_Parents + 1) * sizeof(int), cudaMemcpyHostToDevice);

            cout << "Configuring recombinant progenys' parents" << endl;
            int total_Cells = count_Parents;

            full_Rounds = total_Cells / this->gpu_Limit;
            partial_Rounds = total_Cells % this->gpu_Limit;

            for (int full = 0; full < full_Rounds; full++)
            {
                int start = full * this->gpu_Limit;
                int stop = start + this->gpu_Limit;
                start_stops.push_back(make_pair(start, stop));
            }

            if (partial_Rounds != 0)
            {
                int start = total_Cells - partial_Rounds;
                start_stops.push_back(make_pair(start, total_Cells));
            }

            for (size_t i = 0; i < start_stops.size(); i++)
            {
                int num_of_values_current = start_stops[i].second - start_stops[i].first;
                CUDA_progeny_shuffle_CELLs<<<tot_Blocks, tot_ThreadsperBlock>>>(num_of_values_current, start_stops[i].first, recombination_hotspots, progeny_recom_Index_Cuda, CUDA_progeny_Stride_Index);
                cudaError_t err = cudaGetLastError();

                if (err != cudaSuccess)
                {
                    printf("CUDA Error 1: %s\n", cudaGetErrorString(err));
                    exit(-1);
                    // Possibly: exit(-1) if program cannot continue....
                }
                cudaDeviceSynchronize();
            }
            cout << "Completed configuring recombinant progenys' parents after " << start_stops.size() << " round(s)" << endl;
            start_stops.clear();

            cudaFree(CUDA_progeny_Stride_Index);
        }

        free(progeny_Stride_Index);

        clear_Array_INT(cuda_Progeny_numbers, (recombination_hotspots + 1));
        clear_Array_int_CPU(Progeny_numbers, (recombination_hotspots + 1));

        // int **test_Load = load_to_Host(progeny_recom_Index_Cuda, sum_Progeny, recombination_hotspots + 1);

        // for (size_t row = 0; row < 1; row++)
        // {
        //     for (size_t col = 0; col < (recombination_hotspots + 1); col++)
        //     {
        //         cout << test_Load[row][col] << " ";
        //     }
        //     cout << endl;
        // }
        // clear_Array_int_CPU(test_Load, sum_Progeny);

        //! clear gpu = DONE
        float *CUDA_parent_Proof_reading_probability;
        if (proof_reading_Activate_parent != 0)
        {
            CUDA_parent_Proof_reading_probability = copy_1D_to_CUDA_FLOAT(parent_Proof_reading_probability, count_Parents);
            free(parent_Proof_reading_probability);
        }

        //! CLEAR GPU = done
        float **CUDA_parent_survivability_Probabilities;
        CUDA_parent_survivability_Probabilities = float_2D_Array_load_to_CUDA(parent_survivability_Probabilities, count_Parents, recombination_hotspots + 1);
        clear_Array_float_CPU(parent_survivability_Probabilities, count_Parents);

        //! CLEAR GPU DEFINITELY = DONE
        float **CUDA_current_gen_Progeny_data;
        CUDA_current_gen_Progeny_data = create_CUDA_2D_FLOAT(sum_Progeny, 1 + (3 * recombination_hotspots));
        float *CUDA_progeny_Proof_reading_probability;
        if (proof_reading_Activate_parent != 0)
        {
            cudaMallocManaged(&CUDA_progeny_Proof_reading_probability, sum_Progeny * sizeof(float));
        }

        //! CLEAR GPU = done
        float **CUDA_progeny_survivability_Probabilities = create_CUDA_2D_FLOAT(sum_Progeny, 1 + recombination_hotspots);

        // cout << "Surv 2\n";
        // exit(-1);

        int total_Cells = sum_Progeny;

        full_Rounds = total_Cells / this->gpu_Limit;
        partial_Rounds = total_Cells % this->gpu_Limit;

        for (int full = 0; full < full_Rounds; full++)
        {
            int start = full * this->gpu_Limit;
            int stop = start + this->gpu_Limit;
            start_stops.push_back(make_pair(start, stop));
        }

        if (partial_Rounds != 0)
        {
            int start = total_Cells - partial_Rounds;
            start_stops.push_back(make_pair(start, total_Cells));
        }
        cout << "Configuring progeny profiles" << endl;
        for (size_t i = 0; i < start_stops.size(); i++)
        {
            int num_of_values_current = start_stops[i].second - start_stops[i].first;
            CUDA_progeny_Profiles_fill_CELLs<<<tot_Blocks, tot_ThreadsperBlock>>>(num_of_values_current, CUDA_current_gen_Progeny_data, recombination_hotspots, progeny_recom_Index_Cuda, CUDA_current_gen_Parent_data, CUDA_parent_Proof_reading_probability, CUDA_progeny_Proof_reading_probability, start_stops[i].first, proof_reading_Activate_parent, CUDA_parent_survivability_Probabilities, CUDA_progeny_survivability_Probabilities);
            cudaError_t err = cudaGetLastError();

            if (err != cudaSuccess)
            {
                printf("CUDA Error 1: %s\n", cudaGetErrorString(err));
                exit(-1);
                // Possibly: exit(-1) if program cannot continue....
            }
            cudaDeviceSynchronize();
        }

        cout << "Configured progeny profiles after " << start_stops.size() << " rounds" << endl;
        start_stops.clear();

        // exit(-1);

        clear_Array_FLOAT(CUDA_current_gen_Parent_data, count_Parents);
        if (proof_reading_Activate_parent != 0)
        {
            cudaFree(CUDA_parent_Proof_reading_probability);
        }
        clear_Array_FLOAT(CUDA_parent_survivability_Probabilities, count_Parents);

        // float **test_Load_2 = load_to_Host_FLOAT(CUDA_progeny_survivability_Probabilities, sum_Progeny, recombination_hotspots + 1);

        // for (size_t row = 0; row < 1; row++)
        // {
        //     for (size_t col = 0; col < (recombination_hotspots + 1); col++)
        //     {
        //         cout << test_Load_2[row][col] << " ";
        //     }
        //     cout << endl;
        // }

        // clear_Array_float_CPU(test_Load_2, sum_Progeny);

        // exit(-1);

        // check the above
        // float *test_prob = (float *)malloc(sizeof(float) * sum_Progeny);
        // cudaMemcpy(test_prob, CUDA_progeny_Proof_reading_probability, sizeof(float) * sum_Progeny, cudaMemcpyDeviceToHost);
        // float **test_Load_3 = load_to_Host_FLOAT(CUDA_current_gen_Progeny_data, sum_Progeny, 1 + (3 * recombination_hotspots));
        // for (size_t i = 0; i < 1; i++)
        // {
        //     for (size_t c = 0; c < 1 + (3 * recombination_hotspots); c++)
        //     {
        //         cout << test_Load_3[i][c] << " ";
        //     }
        //     cout << endl;
        // }

        // load parents

        cout << "Loading parent sequences" << endl;
        if (multi_READ == "YES")
        {
            cout << "Initiating multi read" << endl;
            for (int parent = 0; parent < count_Parents; parent++)
            {
                sequences.push_back("");
            }
            int num_per_Core = count_Parents / this->CPU_cores;
            int remainder = count_Parents % this->CPU_cores;

            for (int core_ID = 0; core_ID < this->CPU_cores; core_ID++)
            {
                int start_Cell = core_ID * num_per_Core;
                int stop_Cell = start_Cell + num_per_Core;

                threads_vec.push_back(thread{&functions_library::read_nFASTA_multi_Read_CELLS, this, parent_Sequences_Store, genome_Size, start_Cell, stop_Cell, parent_Indexes_parents_and_their_Cells, generation_Current});
            }
            if (remainder != 0)
            {
                int start_Cell = count_Parents - remainder;
                int stop_Cell = count_Parents;

                threads_vec.push_back(thread{&functions_library::read_nFASTA_multi_Read_CELLS, this, parent_Sequences_Store, genome_Size, start_Cell, stop_Cell, parent_Indexes_parents_and_their_Cells, generation_Current});
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
        else
        {
            cout << "Initiating single read" << endl;
            for (int parent = 0; parent < count_Parents; parent++)
            {
                string name_nFASTA = to_string(generation_Current) + "_" + to_string(parent_Indexes_parents_and_their_Cells[0][parent]);
                string location = parent_Sequences_Store + "/" + name_nFASTA + ".nfasta";

                sequences.push_back(read_nFASTA(location, genome_Size));
            }
        }

        cout << "Loading parents to the GPU" << endl;
        //! CLEAR from gpu = DONE
        int **CUDA_all_Parent_Sequences;
        CUDA_all_Parent_Sequences = create_CUDA_2D_int(count_Parents, genome_Size);

        for (int parent = 0; parent < count_Parents; parent++)
        {
            string sequence = sequences[parent];

            char *reference_full, *cuda_reference;
            reference_full = (char *)malloc((sequence.size() + 1) * sizeof(char));
            cudaMallocManaged(&cuda_reference, (sequence.size() + 1) * sizeof(char));
            strcpy(reference_full, sequence.c_str());
            cudaMemcpy(cuda_reference, reference_full, (sequence.size() + 1) * sizeof(char), cudaMemcpyHostToDevice);

            cuda_fill_parent_Sequences<<<tot_Blocks, tot_ThreadsperBlock>>>(genome_Size, CUDA_all_Parent_Sequences, cuda_reference, parent);

            cudaError_t err = cudaGetLastError();

            if (err != cudaSuccess)
            {
                printf("CUDA Error: %s\n", cudaGetErrorString(err));

                // Possibly: exit(-1) if program cannot continue....
            }
            cudaDeviceSynchronize();

            free(reference_full);
        }
        sequences.clear();
        cout << count_Parents << " parents loaded to the GPU" << endl;

        // int **sequences_Test = load_to_Host(CUDA_all_Parent_Sequences, count_Parents, genome_Size);

        // for (size_t i = 0; i < count_Parents; i++)
        // {
        //     for (size_t c = 0; c < genome_Size; c++)
        //     {
        //         cout << sequences_Test[i][c];
        //     }
        //     cout << endl;
        // }

        cout << "Generating progeny sequences" << endl;
        //! CLEAR FROM GPU = DONE
        int **cuda_progeny_Sequences = create_CUDA_2D_int(sum_Progeny, genome_Size);

        total_Cells = sum_Progeny;

        full_Rounds = total_Cells / this->gpu_Limit;
        partial_Rounds = total_Cells % this->gpu_Limit;

        for (int full = 0; full < full_Rounds; full++)
        {
            int start = full * this->gpu_Limit;
            int stop = start + this->gpu_Limit;
            start_stops.push_back(make_pair(start, stop));
        }

        if (partial_Rounds != 0)
        {
            int start = total_Cells - partial_Rounds;
            start_stops.push_back(make_pair(start, total_Cells));
        }

        for (size_t i = 0; i < start_stops.size(); i++)
        {
            int num_of_values_current = start_stops[i].second - start_stops[i].first;
            CUDA_progeny_sequence_generation_CELLS<<<tot_Blocks, tot_ThreadsperBlock>>>(num_of_values_current, start_stops[i].first, CUDA_all_Parent_Sequences, cuda_progeny_Sequences, progeny_recom_Index_Cuda, recombination_hotspots, genome_Size, CUDA_recombination_hotspots_start_stop);
            cudaError_t err = cudaGetLastError();

            if (err != cudaSuccess)
            {
                printf("CUDA Error 1: %s\n", cudaGetErrorString(err));
                exit(-1);
                // Possibly: exit(-1) if program cannot continue....
            }
            cudaDeviceSynchronize();
        }

        clear_Array_INT(CUDA_all_Parent_Sequences, count_Parents);
        int **progeny_recom_Index = load_to_Host(progeny_recom_Index_Cuda, sum_Progeny, recombination_hotspots + 1);
        clear_Array_INT(progeny_recom_Index_Cuda, sum_Progeny);

        // int **sequences_Test = load_to_Host(cuda_progeny_Sequences, sum_Progeny, genome_Size);
        // cout << "Orginal: " << sequences_Test[0][1999] << endl;
        //  for (size_t i = 0; i < 1; i++)
        //  {
        //      for (size_t c = 0; c < genome_Size; c++)
        //      {
        //
        //      }
        //      cout << endl;
        //  }

        if (mutation_hotspots != 0)
        {
            cout << "Mutating progeny sequences" << endl;
            for (size_t i = 0; i < start_stops.size(); i++)
            {
                int num_of_values_current = start_stops[i].second - start_stops[i].first;
                CUDA_mutate_Progeny_CELLS<<<tot_Blocks, tot_ThreadsperBlock>>>(num_of_values_current, mutation_hotspots, generation_Current,
                                                                               CUDA_mutation_rates_Hotspot_generation, CUDA_mutation_Regions_start_stop,
                                                                               cuda_progeny_Sequences, CUDA_progeny_Proof_reading_probability, proof_reading_Activate_parent,
                                                                               CUDA_A_0_mutation, CUDA_T_1_mutation, CUDA_G_2_mutation, CUDA_C_3_mutation,
                                                                               CUDA_sequence_Mutation_tracker, CUDA_current_gen_Progeny_data,
                                                                               CUDA_A_0_fitness, CUDA_T_1_fitness, CUDA_G_2_fitness, CUDA_C_3_fitness,
                                                                               CUDA_A_0_probability_Proof_reading, CUDA_T_1_probability_Proof_reading, CUDA_G_2_probability_Proof_reading, CUDA_C_3_probability_Proof_reading,
                                                                               CUDA_A_0_Recombination, CUDA_T_1_Recombination, CUDA_G_2_Recombination, CUDA_C_3_Recombination,
                                                                               CUDA_stride_Array, start_stops[i].first,
                                                                               CUDA_progeny_survivability_Probabilities,
                                                                               CUDA_A_0_survivability, CUDA_T_1_survivability, CUDA_G_2_survivability, CUDA_C_3_survivability);
                cudaError_t err = cudaGetLastError();

                if (err != cudaSuccess)
                {
                    printf("CUDA Error 1: %s\n", cudaGetErrorString(err));
                    exit(-1);
                    // Possibly: exit(-1) if program cannot continue....
                }
                cudaDeviceSynchronize();
            }
        }

        // int **sequences_Test_change = load_to_Host(cuda_progeny_Sequences, sum_Progeny, genome_Size);
        // cout << "Change: " << sequences_Test_change[0][1999] << endl;

        // float **test_Load_2 = load_to_Host_FLOAT(CUDA_current_gen_Progeny_data, sum_Progeny, 1 + (3 * recombination_hotspots));
        // for (size_t i = 0; i < 1; i++)
        // {
        //     for (size_t c = 0; c < 1 + (3 * recombination_hotspots); c++)
        //     {
        //         cout << test_Load_2[i][c] << " ";
        //     }
        //     cout << endl;
        // }

        // for (size_t i = 0; i < 1; i++)
        // {
        //     for (size_t c = 0; c < genome_Size; c++)
        //     {
        //         cout << sequences_Test[i][c];
        //     }
        //     cout << endl;
        // }

        start_stops.clear();
        cout << "Mutated progeny sequences" << endl;

        // test_Load_2 = load_to_Host_FLOAT(CUDA_progeny_survivability_Probabilities, sum_Progeny, recombination_hotspots + 1);

        // for (size_t row = 0; row < 1; row++)
        // {
        //     for (size_t col = 0; col < (recombination_hotspots + 1); col++)
        //     {
        //         cout << test_Load_2[row][col] << " ";
        //     }
        //     cout << endl;
        // }

        // exit(-1);

        cout << "Writing data to hard disk" << endl;
        float **current_gen_Progeny_data = load_to_Host_FLOAT(CUDA_current_gen_Progeny_data, sum_Progeny, 1 + (3 * recombination_hotspots));
        clear_Array_FLOAT(CUDA_current_gen_Progeny_data, sum_Progeny);

        float **progeny_survivability_Probabilities = load_to_Host_FLOAT(CUDA_progeny_survivability_Probabilities, sum_Progeny, 1 + recombination_hotspots);
        clear_Array_FLOAT(CUDA_progeny_survivability_Probabilities, sum_Progeny);

        float *progeny_Proof_reading_probability;
        if (proof_reading_Activate_parent != 0)
        {
            progeny_Proof_reading_probability = (float *)malloc(sum_Progeny * sizeof(float));
            cudaMemcpy(progeny_Proof_reading_probability, CUDA_progeny_Proof_reading_probability, sum_Progeny * sizeof(float), cudaMemcpyDeviceToHost);
            cudaFree(CUDA_progeny_Proof_reading_probability);
        }

        cout << "Writing progeny profile data" << endl;

        int progeny_ID = sum_Progeny_in_Generation;

        fstream progeny_File_write;
        progeny_File_write.open(progeny_File, ios::app);
        fstream progney_Recombination_write;
        progney_Recombination_write.open(progeny_Recombination_File, ios::app);
        fstream sequence_Profile_write;
        sequence_Profile_write.open(sequence_Profiles, ios::app);

        fstream progeny_Cells;
        progeny_Cells.open(this->cells_of_progeny, ios::app);

        for (int progeny = 0; progeny < sum_Progeny; progeny++)
        {
            float survivability_Total = progeny_survivability_Probabilities[progeny][0];

            string progeny_Name = to_string(generation_Current + 1) + "_" + to_string(progeny_ID);
            string parent_Name = to_string(generation_Current) + "_" + to_string(parent_Indexes_parents_and_their_Cells[0][progeny_recom_Index[progeny][0]]);
            string cell_Name = to_string(generation_Current) + "_" + to_string(parent_Indexes_parents_and_their_Cells[1][progeny_recom_Index[progeny][0]] + processed_Cells);

            progeny_Cells << progeny_Name << "\t" << cell_Name << "\n";

            fstream progeny_Profile_Files;
            progeny_Profile_Files.open(progeny_Profile_Store + "/" + progeny_Name + ".profile", ios::out);

            fstream progeny_survivability_Files;
            progeny_survivability_Files.open(progeny_Profile_Store + "/" + progeny_Name + ".profile_surv", ios::out);

            progeny_File_write << parent_Name << "\t" << progeny_Name << "\tprimary_parent\n";
            progney_Recombination_write << progeny_Name << "\t" << parent_Name;
            sequence_Profile_write << progeny_Name
                                   << "\t" << to_string(generation_Current + 1)
                                   << "\t" << to_string(progeny_survivability_Probabilities[progeny][0])
                                   << "\t" << to_string(current_gen_Progeny_data[progeny][0]);

            progeny_survivability_Files << to_string(progeny_survivability_Probabilities[progeny][0]);
            progeny_Profile_Files << to_string(current_gen_Progeny_data[progeny][0]);

            if (proof_reading_Activate_parent != 0)
            {
                sequence_Profile_write << "\t" << to_string(progeny_Proof_reading_probability[progeny]);
                fstream progeny_Probability_Files;
                progeny_Probability_Files.open(progeny_Profile_Store + "/" + progeny_Name + ".profile_prob", ios::out);
                progeny_Probability_Files << to_string(progeny_Proof_reading_probability[progeny]);
                progeny_Probability_Files.close();
            }
            else
            {
                // sequence_Profile_write << "\tNA";
            }
            for (int recombination_Hotspot = 0; recombination_Hotspot < recombination_hotspots; recombination_Hotspot++)
            {
                if (progeny_recom_Index[progeny][recombination_Hotspot + 1] != -1)
                {
                    string parent_recomb_Name = to_string(generation_Current) + "_" + to_string(parent_Indexes_parents_and_their_Cells[0][progeny_recom_Index[progeny][recombination_Hotspot + 1]]);
                    progeny_File_write << parent_recomb_Name << "\t" << progeny_Name << "\trecombinant_parent_" << to_string(recombination_Hotspot + 1) << "\n";
                    progney_Recombination_write << "\t" << parent_recomb_Name;
                }
                else
                {
                    progney_Recombination_write << "\t" << parent_Name;
                }
                sequence_Profile_write << "\t" << to_string(current_gen_Progeny_data[progeny][(recombination_Hotspot * 3) + 1])
                                       << "\t" << to_string(current_gen_Progeny_data[progeny][(recombination_Hotspot * 3) + 2])
                                       << "\t" << to_string(current_gen_Progeny_data[progeny][(recombination_Hotspot * 3) + 3])
                                       << "\t" << to_string(progeny_survivability_Probabilities[progeny][recombination_Hotspot + 1]);

                progeny_survivability_Files << "\t" << to_string(progeny_survivability_Probabilities[progeny][recombination_Hotspot + 1]);

                survivability_Total = survivability_Total + progeny_survivability_Probabilities[progeny][recombination_Hotspot + 1];

                progeny_Profile_Files << "\t" << to_string(current_gen_Progeny_data[progeny][(recombination_Hotspot * 3) + 1])
                                      << "\t" << to_string(current_gen_Progeny_data[progeny][(recombination_Hotspot * 3) + 2])
                                      << "\t" << to_string(current_gen_Progeny_data[progeny][(recombination_Hotspot * 3) + 3]);
            }

            // Check whether the progeny live or die

            // int survive = 0;

            if (survivability_Total > 1)
            {
                surviving_Progeny.push_back(progeny);
                sequence_Profile_write << "\tYes";
            }
            else if (survivability_Total < 0)
            {
                // survive = 0;
                sequence_Profile_write << "\tNo";
            }
            else
            {
                random_device rd;
                mt19937 gen(rd()); // Mersenne Twister engine

                // Create a Bernoulli distribution with the given probability
                bernoulli_distribution distribution(survivability_Total);

                // Flip the coin and output the result
                bool result = distribution(gen);
                if (result)
                {
                    surviving_Progeny.push_back(progeny);
                    sequence_Profile_write << "\tYes";
                    // cout << result << " Heads" << endl;
                }
                else
                {
                    // survive = 0;
                    sequence_Profile_write << "\tNo";
                    // cout << result << " Tails" << endl;
                }
            }

            // cout << progeny << "\t" << survive << endl;

            progeny_survivability_Files.close();
            progeny_Profile_Files.close();
            progney_Recombination_write << "\n";
            sequence_Profile_write << "\n";

            progeny_ID++;
        }

        progeny_Cells.close();

        progeny_File_write.close();
        progney_Recombination_write.close();
        sequence_Profile_write.close();

        // exit(-1);

        int **progeny_Sequences = load_to_Host(cuda_progeny_Sequences, sum_Progeny, genome_Size);
        clear_Array_INT(cuda_progeny_Sequences, sum_Progeny);

        if (multi_READ == "YES")
        {
            cout << "Initiating multi write FASTA: " << sum_Progeny << endl;

            int num_per_Core = sum_Progeny / this->CPU_cores;
            int remainder = sum_Progeny % this->CPU_cores;

            for (int core_ID = 0; core_ID < this->CPU_cores; core_ID++)
            {
                int start_Cell = core_ID * num_per_Core;
                int stop_Cell = start_Cell + num_per_Core;

                threads_vec.push_back(thread{&functions_library::write_nFASTA_multi_WRITE, this, progeny_Sequences_Store, sum_Progeny_in_Generation, genome_Size, start_Cell, stop_Cell, progeny_Sequences, generation_Current});
            }

            if (remainder != 0)
            {
                int start_Cell = sum_Progeny - remainder;
                int stop_Cell = sum_Progeny;

                threads_vec.push_back(thread{&functions_library::write_nFASTA_multi_WRITE, this, progeny_Sequences_Store, sum_Progeny_in_Generation, genome_Size, start_Cell, stop_Cell, progeny_Sequences, generation_Current});
            }

            for (thread &t : threads_vec)
            {
                if (t.joinable())
                {
                    t.join();
                }
            }

            threads_vec.clear();

            // exit(-1);
        }
        else
        {
            progeny_ID = sum_Progeny_in_Generation;
            cout << "Initiating single write" << endl;

            for (int progeny_Num = 0; progeny_Num < sum_Progeny; progeny_Num++)
            {
                string progeny_Name = to_string(generation_Current + 1) + "_" + to_string(progeny_ID);
                string file_Location = progeny_Sequences_Store + "/" + progeny_Name + ".nfasta";
                config_File_progeny(file_Location);

                fstream nFASTA;
                nFASTA.open(file_Location, ios::app);

                if (nFASTA.is_open())
                {
                    for (int base = 0; base < genome_Size; base++)
                    {
                        nFASTA << progeny_Sequences[progeny_Num][base];
                    }
                    nFASTA.close();
                }
                progeny_ID++;
            }
        }

        sum_Progeny_in_Generation = sum_Progeny_in_Generation + sum_Progeny;
        processed_Cells = processed_Cells + num_of_Cells;
        clear_Array_int_CPU(progeny_Sequences, sum_Progeny);
        clear_Array_int_CPU(progeny_recom_Index, sum_Progeny);
        if (proof_reading_Activate_parent != 0)
        {
            free(progeny_Proof_reading_probability);
        }
        clear_Array_float_CPU(current_gen_Progeny_data, sum_Progeny);
        clear_Array_float_CPU(progeny_survivability_Probabilities, sum_Progeny);

        // exit(-1);
    }
    else
    {
        for (int cell = start; cell < stop; cell++)
        {
            // int parent_Start = cells_Start_Stop[cell];
            // int parent_Stop = cells_Start_Stop[cell + 1];
            cout << "Per Parent mode" << endl;
            cout << "\nCell ID being processed: " << generation_Current + 1 << "_" << cell << endl;
            cout << "Loading parents" << endl;
            int count_Parents = (cells_Start_Stop[cell + 1] - cells_Start_Stop[cell]);
            cout << "Number of infectious units: " << count_Parents << endl;

            float **current_gen_Parent_data = create_FLOAT_2D_arrays(count_Parents, 1 + (3 * recombination_hotspots));
            float *parent_Proof_reading_probability;
            if (proof_reading_Activate_parent != 0)
            {
                parent_Proof_reading_probability = (float *)malloc(sizeof(float) * count_Parents);
            }

            int *parent_Indexes = (int *)malloc(count_Parents * sizeof(int));

            if (multi_READ == "YES")
            {
                cout << "Initiating multi read" << endl;
                int num_per_Core = count_Parents / this->CPU_cores;
                int remainder = count_Parents % this->CPU_cores;

                for (int core_ID = 0; core_ID < this->CPU_cores; core_ID++)
                {
                    int start_Cell = core_ID * num_per_Core;
                    int stop_Cell = start_Cell + num_per_Core;

                    threads_vec.push_back(thread{&functions_library::read_Profiles_multi_Thread, this, start_Cell, stop_Cell, cells_Start_Stop[cell], parents, current_gen_Parent_data, parent_Proof_reading_probability, generation_Current, parent_Profiles_Store, parent_Indexes});
                }
                if (remainder != 0)
                {
                    int start_Cell = count_Parents - remainder;
                    int stop_Cell = count_Parents;

                    threads_vec.push_back(thread{&functions_library::read_Profiles_multi_Thread, this, start_Cell, stop_Cell, cells_Start_Stop[cell], parents, current_gen_Parent_data, parent_Proof_reading_probability, generation_Current, parent_Profiles_Store, parent_Indexes});
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
            else
            {
                cout << "Initiating single read" << endl;
                int parent_Index = 0;
                for (int parent = cells_Start_Stop[cell]; parent < cells_Start_Stop[cell + 1]; parent++)
                {
                    read_Profiles_Single(parent_Profiles_Store + "/" + to_string(generation_Current) + "_" + to_string(parents[parent]) + ".", current_gen_Parent_data, parent_Proof_reading_probability, parent_Index);
                    parent_Indexes[parent_Index] = parents[parent];
                    parent_Index++;
                }
            }

            // cout << "Cell_ID: " << cell << endl;
            // for (int parent = 0; parent < count_Parents; parent++)
            // {
            //     cout << "parent_ID: " << generation_Current << "_" << parent_Indexes[parent] << endl;
            //     for (size_t i = 0; i < (recombination_hotspots * 3) + 1; i++)
            //     {
            //         cout << current_gen_Parent_data[parent][i] << " ";
            //     }
            //     cout << " Prob proof: " << parent_Proof_reading_probability[parent];
            //     cout << endl;
            // }

            float **CUDA_current_gen_Parent_data = float_2D_Array_load_to_CUDA(current_gen_Parent_data, count_Parents, 1 + (3 * recombination_hotspots));
            clear_Array_float_CPU(current_gen_Parent_data, count_Parents);
            float *CUDA_parent_Proof_reading_probability;
            if (proof_reading_Activate_parent != 0)
            {
                CUDA_parent_Proof_reading_probability = copy_1D_to_CUDA_FLOAT(parent_Proof_reading_probability, count_Parents);
                free(parent_Proof_reading_probability);
            }
            cout << "Parents loaded" << endl;

            cout << "Determining progeny numbers" << endl;

            int **cuda_Progeny_numbers = create_CUDA_2D_int(recombination_hotspots + 1, count_Parents);
            // cudaMallocManaged(&cuda_Progeny_numbers, (count_Parents + 1) * (recombination_hotspots + 1) * sizeof(int));

            // int **tmp = (int **)malloc((recombination_hotspots + 1) * sizeof(tmp[0]));

            // for (int i = 0; i < (recombination_hotspots + 1); i++)
            // {
            //     // cout << i << endl;
            //     cudaMalloc((void **)&tmp[i], (count_Parents + 1) * sizeof(tmp[0][0]));
            // }
            // // cout << "run" << endl;
            // cudaMemcpy(cuda_Progeny_numbers, tmp, (recombination_hotspots + 1) * sizeof(int *), cudaMemcpyHostToDevice);
            // free(tmp);

            int full_Rounds = count_Parents / this->gpu_Limit;
            int partial_Rounds = count_Parents % this->gpu_Limit;

            vector<pair<int, int>> start_stops;

            for (int full = 0; full < full_Rounds; full++)
            {
                int start = full * this->gpu_Limit;
                int stop = start + this->gpu_Limit;
                start_stops.push_back(make_pair(start, stop));
            }

            if (partial_Rounds != 0)
            {
                int start = count_Parents - partial_Rounds;
                start_stops.push_back(make_pair(start, count_Parents));
            }

            for (size_t i = 0; i < start_stops.size(); i++)
            {
                int num_of_values_current = start_stops[i].second - start_stops[i].first;

                curandState *state;
                cudaMalloc((void **)&state, num_of_values_current * sizeof(curandState));

                if (progeny_distribution_Type == "Gamma")
                {
                   // CUDA_gamma_distribution_PROGENY<<<tot_Blocks, tot_ThreadsperBlock>>>(state, num_of_values_current, cuda_Progeny_numbers, progeny_shape, progeny_scale, start_stops[i].first, recombination_hotspots, CUDA_current_gen_Parent_data);
                }
                else if (progeny_distribution_Type == "Negative binomial")
                {
                }
                cudaError_t err = cudaGetLastError();
                if (err != cudaSuccess)
                {
                    printf("CUDA Error: %s\n", cudaGetErrorString(err));

                    // Possibly: exit(-1) if program cannot continue....
                }
                cudaDeviceSynchronize();
                // break;

                cudaFree(state);
            }

            cout << "Completed progeny generation via " << start_stops.size() << " GPU rounds" << endl;
            start_stops.clear();

            int **Progeny_numbers = load_to_Host(cuda_Progeny_numbers, (recombination_hotspots + 1), count_Parents);

            // for (size_t col = 0; col < count_Parents; col++)
            // {
            //     for (size_t row = 0; row < (recombination_hotspots + 1); row++)
            //     {
            //         cout << Progeny_numbers[row][col] << " ";
            //     }
            //     cout << endl;
            // }

            int sum_Progeny = sum_CUDA(cuda_Progeny_numbers, count_Parents);
            cout << "Sum progeny in cell: " << sum_Progeny << endl;
            // sum_Progeny_in_Generation = sum_Progeny_in_Generation + sum_Progeny;
            //  cout << sum_Progeny_in_Generation << endl;

            // create fitness distribution;

            float **CUDA_fitness_distribution;
            CUDA_fitness_distribution = create_CUDA_2D_FLOAT(recombination_hotspots, count_Parents);

            if (recombination_hotspots != 0)
            {
                cout << "Generating selectivity distributions for putative parents" << endl;
                float *CUDA_hotspot_selectivity_Summations;
                cudaMallocManaged(&CUDA_hotspot_selectivity_Summations, recombination_hotspots * sizeof(float));

                // run gpu to get summations
                CUDA_summation_Selectivity<<<tot_Blocks, tot_ThreadsperBlock>>>(recombination_hotspots, count_Parents, CUDA_current_gen_Parent_data, CUDA_hotspot_selectivity_Summations);
                cudaError_t err = cudaGetLastError();

                if (err != cudaSuccess)
                {
                    printf("CUDA Error: %s\n", cudaGetErrorString(err));

                    // Possibly: exit(-1) if program cannot continue....
                }
                cudaDeviceSynchronize();

                // float *summations;
                // summations = (float *)malloc(recombination_hotspots * sizeof(float));
                // cudaMemcpy(summations, CUDA_hotspot_selectivity_Summations, recombination_hotspots * sizeof(float), cudaMemcpyDeviceToHost);

                // for (size_t i = 0; i < recombination_hotspots; i++)
                // {
                //     cout << summations[i] << endl;
                // }

                cout << "Hotspots' summations done" << endl;

                int total_Cells = count_Parents;
                full_Rounds = total_Cells / this->gpu_Limit;
                partial_Rounds = total_Cells % this->gpu_Limit;

                for (int full = 0; full < full_Rounds; full++)
                {
                    int start = full * this->gpu_Limit;
                    int stop = start + this->gpu_Limit;
                    start_stops.push_back(make_pair(start, stop));
                }

                if (partial_Rounds != 0)
                {
                    int start = total_Cells - partial_Rounds;
                    start_stops.push_back(make_pair(start, total_Cells));
                }

                for (size_t i = 0; i < start_stops.size(); i++)
                {
                    int num_of_values_current = start_stops[i].second - start_stops[i].first;

                    CUDA_Selectivity_Distribution<<<tot_Blocks, tot_ThreadsperBlock>>>(recombination_hotspots, num_of_values_current, CUDA_current_gen_Parent_data, CUDA_hotspot_selectivity_Summations, CUDA_fitness_distribution, start_stops[i].first);
                    err = cudaGetLastError();

                    if (err != cudaSuccess)
                    {
                        printf("CUDA Error: %s\n", cudaGetErrorString(err));

                        // Possibly: exit(-1) if program cannot continue....
                    }
                    cudaDeviceSynchronize();
                }

                // calculate fitness distribution
                cout << "Generated selectivity distributions after " << start_stops.size() << " rounds" << endl;
                start_stops.clear();
                cudaFree(CUDA_hotspot_selectivity_Summations);
            }

            // float **test_Array = load_to_Host_FLOAT(CUDA_fitness_distribution, recombination_hotspots, count_Parents);
            // for (size_t row = 0; row < recombination_hotspots; row++)
            // {
            //     for (size_t col = 0; col < count_Parents; col++)
            //     {
            //         cout << test_Array[row][col] << " ";
            //     }
            //     cout << endl;
            // }

            // exit(-1);

            //! load parent sequences
            if (multi_READ == "YES")
            {
                cout << "Initiating multi read" << endl;
                for (int parent = 0; parent < count_Parents; parent++)
                {
                    sequences.push_back("");
                }
                int num_per_Core = count_Parents / this->CPU_cores;
                int remainder = count_Parents % this->CPU_cores;

                for (int core_ID = 0; core_ID < this->CPU_cores; core_ID++)
                {
                    int start_Cell = core_ID * num_per_Core;
                    int stop_Cell = start_Cell + num_per_Core;

                    threads_vec.push_back(thread{&functions_library::read_nFASTA_multi_Read, this, parent_Sequences_Store, genome_Size, start_Cell, stop_Cell, parent_Indexes, generation_Current});
                }
                if (remainder != 0)
                {
                    int start_Cell = count_Parents - remainder;
                    int stop_Cell = count_Parents;

                    threads_vec.push_back(thread{&functions_library::read_nFASTA_multi_Read, this, parent_Sequences_Store, genome_Size, start_Cell, stop_Cell, parent_Indexes, generation_Current});
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
            else
            {
                cout << "Initiating single read" << endl;
                for (int parent = 0; parent < count_Parents; parent++)
                {
                    string name_nFASTA = to_string(generation_Current) + "_" + to_string(parent_Indexes[parent]);
                    string location = parent_Sequences_Store + "/" + name_nFASTA + ".nfasta";

                    sequences.push_back(read_nFASTA(location, genome_Size));
                }
            }

            // for (int parent = 0; parent < count_Parents; parent++)
            // {
            //     cout << "File name: " << to_string(generation_Current) + "_" + to_string(parent_Indexes[parent]) << endl;
            //     cout << sequences[parent] << endl;
            //     cout << endl;
            // }

            cout << "Loading parents to the GPU" << endl;
            int **CUDA_all_Parent_Sequences;
            CUDA_all_Parent_Sequences = create_CUDA_2D_int(count_Parents, genome_Size);

            for (int parent = 0; parent < count_Parents; parent++)
            {
                string sequence = sequences[parent];

                char *reference_full, *cuda_reference;
                reference_full = (char *)malloc((sequence.size() + 1) * sizeof(char));
                cudaMallocManaged(&cuda_reference, (sequence.size() + 1) * sizeof(char));
                strcpy(reference_full, sequence.c_str());
                cudaMemcpy(cuda_reference, reference_full, (sequence.size() + 1) * sizeof(char), cudaMemcpyHostToDevice);

                cuda_fill_parent_Sequences<<<tot_Blocks, tot_ThreadsperBlock>>>(genome_Size, CUDA_all_Parent_Sequences, cuda_reference, parent);

                cudaError_t err = cudaGetLastError();

                if (err != cudaSuccess)
                {
                    printf("CUDA Error: %s\n", cudaGetErrorString(err));

                    // Possibly: exit(-1) if program cannot continue....
                }
                cudaDeviceSynchronize();

                free(reference_full);
            }

            sequences.clear();
            cout << count_Parents << " parents loaded to the GPU" << endl;
            // exit(-1);

            for (int parent = 0; parent < count_Parents; parent++)
            {
                cout << "\nProcessing parent " << parent + 1 << " of " << count_Parents << endl;
                int **progeny_recom_Index_Cuda = create_CUDA_2D_int(Progeny_numbers[0][parent], recombination_hotspots + 1);
                int total_Cells = Progeny_numbers[0][parent];

                //  cout << total_Cells << endl;

                full_Rounds = total_Cells / this->gpu_Limit;
                partial_Rounds = total_Cells % this->gpu_Limit;

                start_stops.clear();

                for (int full = 0; full < full_Rounds; full++)
                {
                    int start = full * this->gpu_Limit;
                    int stop = start + this->gpu_Limit;
                    start_stops.push_back(make_pair(start, stop));
                }

                if (partial_Rounds != 0)
                {
                    int start = total_Cells - partial_Rounds;
                    start_stops.push_back(make_pair(start, total_Cells));
                }

                for (size_t i = 0; i < start_stops.size(); i++)
                {
                    int num_of_values_current = start_stops[i].second - start_stops[i].first;
                    // cout << num_of_values_current<< endl;
                    int fill_Value = -1;
                    progeny_Array_CUDA<<<tot_Blocks, tot_ThreadsperBlock>>>(num_of_values_current, parent, recombination_hotspots, cuda_Progeny_numbers, progeny_recom_Index_Cuda, fill_Value, start_stops[i].first, CUDA_fitness_distribution, count_Parents);
                    cudaError_t err = cudaGetLastError();

                    if (err != cudaSuccess)
                    {
                        printf("CUDA Error 1: %s\n", cudaGetErrorString(err));
                        exit(-1);
                        // Possibly: exit(-1) if program cannot continue....
                    }
                    cudaDeviceSynchronize();
                }

                cout << "Created progeny index for parent " << parent + 1 << " after " << start_stops.size() << " round(s)." << endl;

                // int **test_Load = load_to_Host(progeny_recom_Index_Cuda, total_Cells, recombination_hotspots + 1);
                // for (size_t row = 0; row < total_Cells; row++)
                // {
                //     for (size_t col = 0; col < (recombination_hotspots + 1); col++)
                //     {
                //         cout << test_Load[row][col] << " ";
                //     }
                //     cout << endl;
                // }
                // free(test_Load);

                if (recombination_hotspots != 0)
                {
                    cout << "Configuring recombinant progenys' parents" << endl;
                    // randomise
                    CUDA_progeny_shuffle<<<tot_Blocks, tot_ThreadsperBlock>>>(recombination_hotspots, progeny_recom_Index_Cuda, Progeny_numbers[0][parent]);
                    cudaError_t err = cudaGetLastError();

                    if (err != cudaSuccess)
                    {
                        printf("CUDA Error 2: %s\n", cudaGetErrorString(err));
                        exit(-1);
                        // Possibly: exit(-1) if program cannot continue....
                    }
                    cudaDeviceSynchronize();

                    // int** test_Load = load_to_Host(progeny_recom_Index_Cuda, total_Cells, recombination_hotspots + 1);
                    // for (size_t row = 0; row < total_Cells; row++)
                    // {
                    //     for (size_t col = 0; col < (recombination_hotspots + 1); col++)
                    //     {
                    //         cout << test_Load[row][col] << " ";
                    //     }
                    //     cout << endl;
                    // }
                    // free(test_Load);
                }

                // create and write the sequences for this parents progeny a round at a time
                string parent_Name = to_string(generation_Current) + "_" + to_string(parent_Indexes[parent]);
                int **progeny_recom_Index = load_to_Host(progeny_recom_Index_Cuda, Progeny_numbers[0][parent], recombination_hotspots + 1);

                for (size_t i = 0; i < start_stops.size(); i++)
                {
                    int num_of_values_current = start_stops[i].second - start_stops[i].first;

                    cout << "Configuring progeny profiles" << endl;
                    float **CUDA_current_gen_Progeny_data = create_CUDA_2D_FLOAT(num_of_values_current, 1 + (3 * recombination_hotspots));
                    float *CUDA_progeny_Proof_reading_probability;
                    if (proof_reading_Activate_parent != 0)
                    {
                        cudaMallocManaged(&CUDA_progeny_Proof_reading_probability, num_of_values_current * sizeof(float));
                    }
                    CUDA_progeny_Profiles_fill<<<tot_Blocks, tot_ThreadsperBlock>>>(num_of_values_current, CUDA_current_gen_Progeny_data, recombination_hotspots, progeny_recom_Index_Cuda, CUDA_current_gen_Parent_data, CUDA_parent_Proof_reading_probability, CUDA_progeny_Proof_reading_probability, start_stops[i].first, proof_reading_Activate_parent);

                    cudaError_t err = cudaGetLastError();

                    if (err != cudaSuccess)
                    {
                        printf("CUDA Error 3: %s\n", cudaGetErrorString(err));
                        // cout << i << endl;
                        exit(-1);
                        // Possibly: exit(-1) if program cannot continue....
                    }

                    cudaDeviceSynchronize();

                    // for (int row = 0; row < 1; row++)
                    // {
                    //     for (int column = 0; column < (1 + (3 * recombination_hotspots)); column++)
                    //     {
                    //         cout << current_gen_Progeny_data[row][column] << " ";
                    //     }
                    //     cout << endl;
                    // }

                    // free(current_gen_Parent_data);

                    // exit(-1);
                    cout << "Generating progeny sequences" << endl;
                    int **cuda_progeny_Sequences = create_CUDA_2D_int(num_of_values_current, genome_Size);

                    CUDA_progeny_sequence_generation<<<tot_Blocks, tot_ThreadsperBlock>>>(num_of_values_current, start_stops[i].first, CUDA_all_Parent_Sequences, cuda_progeny_Sequences, progeny_recom_Index_Cuda, recombination_hotspots, genome_Size, CUDA_recombination_hotspots_start_stop);

                    err = cudaGetLastError();

                    if (err != cudaSuccess)
                    {
                        printf("CUDA Error 4: %s\n", cudaGetErrorString(err));
                        exit(-1);
                        // Possibly: exit(-1) if program cannot continue....
                    }
                    cudaDeviceSynchronize();

                    // int **progeny_Sequences = load_to_Host(cuda_progeny_Sequences, num_of_values_current, genome_Size);
                    // cout << "Original: " << progeny_Sequences[0][1999] << endl;
                    // cout << "Original: " << progeny_Sequences[0][149] << endl;
                    // float *progeny_Proof_reading_probability = (float *)malloc(num_of_values_current * sizeof(float));
                    // cudaMemcpy(progeny_Proof_reading_probability, CUDA_progeny_Proof_reading_probability, num_of_values_current * sizeof(float), cudaMemcpyDeviceToHost);
                    // cout << progeny_Proof_reading_probability[0] << endl;
                    // free(progeny_Proof_reading_probability);
                    // for (size_t row = 0; row < num_of_values_current; row++)
                    // {
                    //     for (size_t col = 0; col < genome_Size; col++)
                    //     {
                    //         cout << progeny_Sequences[row][col];
                    //     }
                    //     cout << endl;
                    // }
                    // free(progeny_Sequences);
                    // exit(-1);

                    if (mutation_hotspots != 0)
                    {
                        cout << "Mutating progeny sequences" << endl;
                        CUDA_mutate_Progeny<<<tot_Blocks, tot_ThreadsperBlock>>>(num_of_values_current, mutation_hotspots, generation_Current,
                                                                                 CUDA_mutation_rates_Hotspot_generation, CUDA_mutation_Regions_start_stop,
                                                                                 cuda_progeny_Sequences, CUDA_progeny_Proof_reading_probability, proof_reading_Activate_parent,
                                                                                 CUDA_A_0_mutation, CUDA_T_1_mutation, CUDA_G_2_mutation, CUDA_C_3_mutation,
                                                                                 CUDA_sequence_Mutation_tracker, CUDA_current_gen_Progeny_data,
                                                                                 CUDA_A_0_fitness, CUDA_T_1_fitness, CUDA_G_2_fitness, CUDA_C_3_fitness,
                                                                                 CUDA_A_0_probability_Proof_reading, CUDA_T_1_probability_Proof_reading, CUDA_G_2_probability_Proof_reading, CUDA_C_3_probability_Proof_reading,
                                                                                 CUDA_A_0_Recombination, CUDA_T_1_Recombination, CUDA_G_2_Recombination, CUDA_C_3_Recombination,
                                                                                 CUDA_stride_Array);
                        err = cudaGetLastError();

                        if (err != cudaSuccess)
                        {
                            printf("CUDA Error 5: %s\n", cudaGetErrorString(err));
                            exit(-1);
                            // Possibly: exit(-1) if program cannot continue....
                        }
                        cudaDeviceSynchronize();
                    }

                    cout << "Mutated progeny sequences" << endl;

                    // Write sequences and intermediates
                    cout << "Writing data to hard disk" << endl;
                    float **current_gen_Progeny_data = load_to_Host_FLOAT(CUDA_current_gen_Progeny_data, num_of_values_current, 1 + (3 * recombination_hotspots));
                    clear_Array_FLOAT(CUDA_current_gen_Progeny_data, num_of_values_current);
                    // cudaFree(CUDA_current_gen_Progeny_data);

                    float *progeny_Proof_reading_probability;
                    if (proof_reading_Activate_parent != 0)
                    {
                        progeny_Proof_reading_probability = (float *)malloc(num_of_values_current * sizeof(float));
                        cudaMemcpy(progeny_Proof_reading_probability, CUDA_progeny_Proof_reading_probability, num_of_values_current * sizeof(float), cudaMemcpyDeviceToHost);
                        cudaFree(CUDA_progeny_Proof_reading_probability);
                    }

                    cout << "Writing progeny profile data" << endl;

                    int progeny_ID = sum_Progeny_in_Generation;

                    fstream progeny_File_write;
                    progeny_File_write.open(progeny_File, ios::app);
                    fstream progney_Recombination_write;
                    progney_Recombination_write.open(progeny_Recombination_File, ios::app);
                    fstream sequence_Profile_write;
                    sequence_Profile_write.open(sequence_Profiles, ios::app);

                    for (int progeny = 0; progeny < num_of_values_current; progeny++)
                    {

                        string progeny_Name = to_string(generation_Current + 1) + "_" + to_string(progeny_ID);

                        fstream progeny_Profile_Files;
                        progeny_Profile_Files.open(progeny_Profile_Store + "/" + progeny_Name + ".profile", ios::out);

                        progeny_File_write << parent_Name << "\t" << progeny_Name << "\tprimary_parent\n";
                        progney_Recombination_write << progeny_Name << "\t" << parent_Name;
                        sequence_Profile_write << progeny_Name
                                               << "\t" << to_string(generation_Current + 1)
                                               << "\t" << to_string(current_gen_Progeny_data[progeny][0]);
                        progeny_Profile_Files << to_string(current_gen_Progeny_data[progeny][0]);

                        if (proof_reading_Activate_parent != 0)
                        {
                            sequence_Profile_write << "\t" << to_string(progeny_Proof_reading_probability[progeny]);
                            fstream progeny_Probability_Files;
                            progeny_Probability_Files.open(progeny_Profile_Store + "/" + progeny_Name + ".profile_prob", ios::out);
                            progeny_Probability_Files << to_string(progeny_Proof_reading_probability[progeny]);
                            progeny_Probability_Files.close();
                        }
                        else
                        {
                            sequence_Profile_write << "\tNA";
                        }

                        for (int recombination_Hotspot = 0; recombination_Hotspot < recombination_hotspots; recombination_Hotspot++)
                        {
                            if (progeny_recom_Index[progeny][recombination_Hotspot + 1] != -1)
                            {
                                string parent_recomb_Name = to_string(generation_Current) + "_" + to_string(parent_Indexes[progeny_recom_Index[progeny][recombination_Hotspot + 1]]);
                                progeny_File_write << parent_recomb_Name << "\t" << progeny_Name << "\trecombinant_parent_" << to_string(recombination_Hotspot + 1) << "\n";
                                progney_Recombination_write << "\t" << parent_recomb_Name;
                            }
                            else
                            {
                                progney_Recombination_write << "\t" << parent_Name;
                            }
                            sequence_Profile_write << "\t" << to_string(current_gen_Progeny_data[progeny][(recombination_Hotspot * 3) + 1])
                                                   << "\t" << to_string(current_gen_Progeny_data[progeny][(recombination_Hotspot * 3) + 2])
                                                   << "\t" << to_string(current_gen_Progeny_data[progeny][(recombination_Hotspot * 3) + 3]);

                            progeny_Profile_Files << "\t" << to_string(current_gen_Progeny_data[progeny][(recombination_Hotspot * 3) + 1])
                                                  << "\t" << to_string(current_gen_Progeny_data[progeny][(recombination_Hotspot * 3) + 2])
                                                  << "\t" << to_string(current_gen_Progeny_data[progeny][(recombination_Hotspot * 3) + 3]);
                        }
                        progeny_Profile_Files.close();
                        progney_Recombination_write << "\n";
                        sequence_Profile_write << "\n";

                        progeny_ID++;
                    }

                    progeny_File_write.close();
                    progney_Recombination_write.close();
                    sequence_Profile_write.close();

                    int **progeny_Sequences = load_to_Host(cuda_progeny_Sequences, num_of_values_current, genome_Size);
                    clear_Array_INT(cuda_progeny_Sequences, num_of_values_current);
                    // cudaFree(cuda_progeny_Sequences);

                    if (multi_READ == "YES")
                    {
                        cout << "Initiating multi write" << endl;

                        int num_per_Core = num_of_values_current / this->CPU_cores;
                        int remainder = num_of_values_current % this->CPU_cores;

                        for (int core_ID = 0; core_ID < this->CPU_cores; core_ID++)
                        {
                            int start_Cell = core_ID * num_per_Core;
                            int stop_Cell = start_Cell + num_per_Core;

                            threads_vec.push_back(thread{&functions_library::write_nFASTA_multi_WRITE, this, progeny_Sequences_Store, sum_Progeny_in_Generation, genome_Size, start_Cell, stop_Cell, progeny_Sequences, generation_Current});
                        }

                        if (remainder != 0)
                        {
                            int start_Cell = num_of_values_current - remainder;
                            int stop_Cell = num_of_values_current;

                            threads_vec.push_back(thread{&functions_library::write_nFASTA_multi_WRITE, this, progeny_Sequences_Store, sum_Progeny_in_Generation, genome_Size, start_Cell, stop_Cell, progeny_Sequences, generation_Current});
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
                    else
                    {
                        progeny_ID = sum_Progeny_in_Generation;
                        cout << "Initiating single write" << endl;
                        for (int progeny_Num = 0; progeny_Num < num_of_values_current; progeny_Num++)
                        {
                            string progeny_Name = to_string(generation_Current + 1) + "_" + to_string(progeny_ID);
                            string file_Location = progeny_Sequences_Store + "/" + progeny_Name + ".nfasta";
                            config_File_progeny(file_Location);

                            fstream nFASTA;
                            nFASTA.open(file_Location, ios::app);

                            if (nFASTA.is_open())
                            {
                                for (int base = 0; base < genome_Size; base++)
                                {
                                    nFASTA << progeny_Sequences[progeny_Num][base];
                                }
                                nFASTA.close();
                            }
                            progeny_ID++;
                        }
                    }

                    //     // progeny_Sequences = load_to_Host(cuda_progeny_Sequences, num_of_values_current, genome_Size);
                    //     // cout << "Mutated: " << progeny_Sequences[0][1999] << endl;
                    //     // cout << "Mutated: " << progeny_Sequences[0][149] << endl;
                    //     // progeny_Proof_reading_probability = (float *)malloc(num_of_values_current * sizeof(float));
                    //     // cudaMemcpy(progeny_Proof_reading_probability, CUDA_progeny_Proof_reading_probability, num_of_values_current * sizeof(float), cudaMemcpyDeviceToHost);
                    //     // cout << progeny_Proof_reading_probability[0] << endl;
                    //     // free(progeny_Proof_reading_probability);
                    //     // current_gen_Progeny_data = load_to_Host_FLOAT(CUDA_current_gen_Progeny_data, num_of_values_current, 1 + (3 * recombination_hotspots));

                    //     // for (int row = 0; row < 1; row++)
                    //     // {
                    //     //     for (int column = 0; column < (1 + (3 * recombination_hotspots)); column++)
                    //     //     {
                    //     //         cout << current_gen_Progeny_data[row][column] << " ";
                    //     //     }
                    //     //     cout << endl;
                    //     // }
                    //     // exit(-1);
                    //     // after writing all
                    sum_Progeny_in_Generation = sum_Progeny_in_Generation + num_of_values_current;
                }
                // start_stops.clear();
                // cudaFree(progeny_recom_Index_Cuda);
                clear_Array_INT(progeny_recom_Index_Cuda, Progeny_numbers[0][parent]);
                free(progeny_recom_Index);
            }

            // exit(-1);

            // cudaFree(cuda_Progeny_numbers);
            clear_Array_INT(cuda_Progeny_numbers, recombination_hotspots + 1);
            // cudaFree(CUDA_all_Parent_Sequences);
            clear_Array_INT(CUDA_all_Parent_Sequences, count_Parents);
            cudaFree(CUDA_parent_Proof_reading_probability);

            if (recombination_hotspots != 0)
            {
                // cudaFree(CUDA_fitness_distribution);
                clear_Array_FLOAT(CUDA_fitness_distribution, recombination_hotspots);
            }

            if (proof_reading_Activate_parent != 0)
            {
                cudaFree(CUDA_parent_Proof_reading_probability);
            }

            free(Progeny_numbers);

            cout << endl;
        }
    }
}

int **functions_library::create_CUDA_2D_int(int rows, int columns)
{
    cudaSetDevice(CUDA_device_IDs[0]);
    int **CUDA_array;

    cudaMallocManaged(&CUDA_array, rows * sizeof(int *));
    for (int row = 0; row < rows; row++)
    {
        cudaMalloc((void **)&(CUDA_array[row]), (1 + columns) * sizeof(int));
    }

    // cudaMallocManaged(&CUDA_array, (columns + 1) * rows * sizeof(int));
    // int **tmp = (int **)malloc(rows * sizeof(tmp[0]));
    // for (int i = 0; i < rows; i++)
    // {
    //     // cout << i << endl;
    //     cudaMalloc((void **)&tmp[i], (columns + 1) * sizeof(tmp[0][0]));
    // }
    // cudaMemcpy(CUDA_array, tmp, rows * sizeof(int *), cudaMemcpyHostToDevice);
    // free(tmp);

    // for (int i = 0; i < rows; i++)
    // {
    //     free(tmp[i]);
    // }

    return CUDA_array;
}

float **functions_library::create_CUDA_2D_FLOAT(int rows, int columns)
{
    cudaSetDevice(CUDA_device_IDs[0]);
    float **CUDA_array;

    cudaMallocManaged(&CUDA_array, rows * sizeof(float *));
    for (int row = 0; row < rows; row++)
    {
        cudaMalloc((void **)&(CUDA_array[row]), (1 + columns) * sizeof(float));
    }

    // cudaMallocManaged(&CUDA_array, (columns + 1) * rows * sizeof(float));
    // float **tmp = (float **)malloc(rows * sizeof(tmp[0]));
    // for (int i = 0; i < rows; i++)
    // {
    //     // cout << i << endl;
    //     cudaMalloc((void **)&tmp[i], (columns + 1) * sizeof(tmp[0][0]));
    // }
    // cudaMemcpy(CUDA_array, tmp, rows * sizeof(float *), cudaMemcpyHostToDevice);

    // // for (int i = 0; i < rows; i++)
    // // {
    // //     free(tmp[i]);
    // // }

    // free(tmp);

    return CUDA_array;
}

float **functions_library::load_to_Host_FLOAT(float **cuda_2D_Array, int rows, int columns)
{
    cudaSetDevice(CUDA_device_IDs[0]);
    float **Array_host_2D = create_FLOAT_2D_arrays(rows, columns);

    for (int row = 0; row < rows; row++)
    {
        cudaMemcpy(Array_host_2D[row], cuda_2D_Array[row], (columns + 1) * sizeof(cuda_2D_Array[0][0]), cudaMemcpyDeviceToHost);
    }

    return Array_host_2D;
}

int **functions_library::load_to_Host(int **cuda_2D_Array, int rows, int columns)
{
    cudaSetDevice(CUDA_device_IDs[0]);
    int **Array_host_2D = create_INT_2D_arrays(rows, columns);

    for (int row = 0; row < rows; row++)
    {
        cudaMemcpy(Array_host_2D[row], cuda_2D_Array[row], (columns + 1) * sizeof(cuda_2D_Array[0][0]), cudaMemcpyDeviceToHost);
    }

    return Array_host_2D;
}

void functions_library::read_Profiles_Single(string file_Names, float **current_gen_Parent_data, float *parent_Proof_reading_probability, int &parent_Index)
{
    fstream parent_Profile;
    parent_Profile.open(file_Names + "profile", ios::in);
    string parent_Profile_Line;
    getline(parent_Profile, parent_Profile_Line);

    vector<string> line_Data;
    split(line_Data, parent_Profile_Line, '\t');

    for (int i = 0; i < line_Data.size(); i++)
    {
        current_gen_Parent_data[parent_Index][i] = stof(line_Data[i]);
    }

    parent_Profile.close();

    if (proof_reading_Activate_parent != 0)
    {
        fstream parent_Proof_Reading_prob;
        parent_Proof_Reading_prob.open(file_Names + "profile_prob", ios::in);
        string parent_Proof_Prob_Line;
        getline(parent_Proof_Reading_prob, parent_Proof_Prob_Line);
        parent_Proof_reading_probability[parent_Index] = stof(parent_Proof_Prob_Line);

        parent_Proof_Reading_prob.close();
    }
}

void functions_library::read_Profiles_multi_Thread(int start, int stop, int cell_Start, int *parents,
                                                   float **current_gen_Parent_data, float *parent_Proof_reading_probability,
                                                   int generation_Current, string parent_Profiles_Store,
                                                   int *parent_Indexes)
{
    for (int i = start; i < stop; i++)
    {
        string file_Names = parent_Profiles_Store + "/" + to_string(generation_Current) + "_" + to_string(parents[cell_Start + i]) + ".";
        parent_Indexes[i] = parents[cell_Start + i];
        fstream parent_Profile;
        parent_Profile.open(file_Names + "profile", ios::in);
        string parent_Profile_Line;
        getline(parent_Profile, parent_Profile_Line);

        vector<string> line_Data;
        split(line_Data, parent_Profile_Line, '\t');

        for (int col = 0; col < line_Data.size(); col++)
        {
            current_gen_Parent_data[i][col] = stof(line_Data[col]);
            // cout << current_gen_Parent_data[i][col] << endl;
        }

        parent_Profile.close();

        if (proof_reading_Activate_parent != 0)
        {
            fstream parent_Proof_Reading_prob;
            parent_Proof_Reading_prob.open(file_Names + "profile_prob", ios::in);
            string parent_Proof_Prob_Line;
            getline(parent_Proof_Reading_prob, parent_Proof_Prob_Line);
            parent_Proof_reading_probability[i] = stof(parent_Proof_Prob_Line);
            parent_Proof_Reading_prob.close();
        }
    }
}

string functions_library::clean_Line(string &line)
{
    string output_Line = "";

    for (int i = 0; i < line.size(); i++)
    {
        if (isascii(line[i]))
        {
            output_Line = output_Line + line[i];
        }
    }

    return output_Line;
}

string functions_library::clean_Invisible(string line)
{
    string output_Line = "";

    for (int i = 0; i < line.size(); i++)
    {
        char ch = line.at(i);
        if (ch == '\t' || ch == '\r' || ch == '\n' || ch == '\v' ||
            ch == '\f' || ch == '\b' || ch == '\a' || ch == '\0')
        {
            // Skip invisible characters
        }
        else
        {
            output_Line = output_Line + line[i];
        }
    }

    return output_Line;
}

__global__ void cuda_Sequences_to_INT(int num_Sequences, int **sequence_INT, int genome_Length, char *sites)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < num_Sequences)
    {
        int site_Start = tid * genome_Length;
        int site_End = site_Start + genome_Length;

        int bp_Pos = 0;

        for (int site = site_Start; site < site_End; site++)
        {
            if (sites[site] == 'A' || sites[site] == 'a' || sites[site] == '0')
            {
                sequence_INT[tid][bp_Pos] = 0;
            }
            else if (sites[site] == 'T' || sites[site] == 't' || sites[site] == '1')
            {
                sequence_INT[tid][bp_Pos] = 1;
            }
            else if (sites[site] == 'G' || sites[site] == 'g' || sites[site] == '2')
            {
                sequence_INT[tid][bp_Pos] = 2;
            }
            else if (sites[site] == 'C' || sites[site] == 'c' || sites[site] == '3')
            {
                sequence_INT[tid][bp_Pos] = 3;
            }

            bp_Pos++;
        }

        tid += blockDim.x * gridDim.x;
    }
}

int **functions_library::process_Reference_Sequences(vector<string> collect_Sequences, int &genome_Length, int &num_of_Sequences_current)
{
    // int num_of_Sequences_current = start_stops[round].second - start_stops[round].first;

    cout << "Configuring multi gpu distribution of " << num_of_Sequences_current << " sequence(s)\n";

    int standard_num_per_GPU = num_of_Sequences_current / num_Cuda_devices;
    int remainder = num_of_Sequences_current % num_Cuda_devices;

    vector<pair<int, int>> start_stop_Per_GPU;

    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
    {
        int start = gpu * standard_num_per_GPU;
        int stop = start + standard_num_per_GPU;

        start_stop_Per_GPU.push_back(make_pair(start, stop));
    }

    start_stop_Per_GPU[num_Cuda_devices - 1].second = start_stop_Per_GPU[num_Cuda_devices - 1].second + remainder;

    string all_Sequences = "";

    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
    {
        for (int sequence = start_stop_Per_GPU[gpu].first; sequence < start_stop_Per_GPU[gpu].second; sequence++)
        {
            all_Sequences.append(collect_Sequences[sequence]);
        }
    }

    // int site_Index[num_of_Sequences_current + 1];
    // site_Index[0] = 0;

    // for (int i = 0; i < num_of_Sequences_current; i++)
    // {
    //     all_Sequences.append(collect_Sequences[start_stops[round].first + i]);
    //     site_Index[i + 1] = site_Index[i] + collect_Sequences[start_stops[round].first + i].size();
    // }

    char *full_Char;
    full_Char = (char *)malloc((all_Sequences.size() + 1) * sizeof(char));
    // sequence = (int *)malloc((all_Sequences.size() + 1) * sizeof(int));
    strcpy(full_Char, all_Sequences.c_str());

    // cout << "Genome length: " << genome_Length << endl;
    //  exit(-1);

    // for (size_t i = 0; i < num_of_Sequences_current; i++)
    // {
    //     int start = i * genome_Length;
    //     for (size_t s = start; s < start + genome_Length; s++)
    //     {
    //         cout << full_Char[s];
    //     }
    //     cout << "\n\n";
    // }

    cudaStream_t streams[num_Cuda_devices];
    cudaDeviceProp deviceProp;

    char *cuda_full_Char[num_Cuda_devices];
    int **cuda_Sequence[num_Cuda_devices];

    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
    {
        cudaSetDevice(CUDA_device_IDs[gpu]);
        cudaGetDeviceProperties(&deviceProp, gpu);
        cout << "Intializing GPU " << CUDA_device_IDs[gpu] << "'s stream: " << deviceProp.name << endl;

        cudaMalloc(&cuda_full_Char[gpu], (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first) * genome_Length * sizeof(char));
        cudaMemcpy(cuda_full_Char[gpu], full_Char + (start_stop_Per_GPU[gpu].first * genome_Length), (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first) * genome_Length * sizeof(char), cudaMemcpyHostToDevice);

        // cudaMalloc(&cuda_Sequence[gpu], (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first) * genome_Length * sizeof(int));

        cudaMallocManaged(&cuda_Sequence[gpu], (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first) * sizeof(int *));
        for (int row = 0; row < (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first); row++)
        {
            cudaMalloc((void **)&cuda_Sequence[gpu][row], genome_Length * sizeof(int));
        }

        cudaStreamCreate(&streams[gpu]);
    }

    cout << "Loaded " << num_of_Sequences_current << " sequence(s) to the GPU(s)\n";

    // int devID[num_Cuda_devices];

    // cout << "Run\n";

    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
    {
        cudaSetDevice(CUDA_device_IDs[gpu]);
        cuda_Sequences_to_INT<<<tot_Blocks_array[gpu], tot_ThreadsperBlock_array[gpu], 0, streams[gpu]>>>(start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first, cuda_Sequence[gpu], genome_Length, cuda_full_Char[gpu]);
    }

    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
    {
        cudaSetDevice(CUDA_device_IDs[gpu]);
        cudaStreamSynchronize(streams[gpu]);
    }

    cout << "GPU(s) streams completed and synchronized\nCopying data from GPU to Host memory\n";

    int **sequence;
    sequence = create_INT_2D_arrays(num_of_Sequences_current, genome_Length);

    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
    {
        cudaSetDevice(CUDA_device_IDs[gpu]);
        // cudaMemcpy(sequence + (start_stop_Per_GPU[gpu].first * genome_Length), cuda_Sequence[gpu], (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first) * genome_Length * sizeof(int), cudaMemcpyDeviceToHost);
        // cudaMemcpy(sequence[start_stop_Per_GPU[gpu].first], cuda_Sequence[gpu], (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first) * sizeof(int *), cudaMemcpyDeviceToHost);

        for (int row = 0; row < (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first); row++)
        {
            cudaMemcpy(sequence[start_stop_Per_GPU[gpu].first + row], cuda_Sequence[gpu][row], genome_Length * sizeof(int), cudaMemcpyDeviceToHost);
        }
    }

    cout << "Data received by host\n";

    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
    {
        cudaSetDevice(CUDA_device_IDs[gpu]);
        cudaFree(cuda_full_Char[gpu]);
        // cudaFree(cuda_Sequence);

        for (int row = 0; row < (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first); row++)
        {
            cudaFree(cuda_Sequence[gpu][row]);
        }
        cudaFree(cuda_Sequence[gpu]);

        cudaStreamDestroy(streams[gpu]);
    }

    // cudaError_t err = cudaGetLastError();

    // if (err != cudaSuccess)
    // {
    //     printf("CUDA Error: %s\n", cudaGetErrorString(err));
    //     exit(-1);
    //     // Possibly: exit(-1) if program cannot continue....
    // }
    // else
    // {
    //     cout << "OK\n";
    // }

    // exit(-1);

    // for (int row = 0; row < num_of_Sequences_current; row++)
    // {
    //     for (size_t c = 0; c < genome_Length; c++)
    //     {
    //         cout << sequence[row][c];
    //     }
    //     cout << "\n\n";
    // }

    return sequence;

    // for (size_t i = 0; i < num_of_Sequences_current; i++)
    // {
    //     int start = i * genome_Length;
    //     for (size_t s = start; s < start + genome_Length; s++)
    //     {
    //         cout << sequence[s];
    //     }
    //     cout << "\n\n";
    // }
}

vector<string> functions_library::convert_Sequences_Master(int **sequences, int &genome_Length, int &num_of_Sequences_current)
{
    cout << "Converting sequences to strings\n";
    for (int sequence = 0; sequence < num_of_Sequences_current; sequence++)
    {
        all_sequences_String.push_back("");
    }

    int num_per_Core = num_of_Sequences_current / this->CPU_cores;
    int remainder = num_of_Sequences_current % this->CPU_cores;

    vector<thread> threads_vec;

    for (int core_ID = 0; core_ID < this->CPU_cores; core_ID++)
    {
        int start_Cell = core_ID * num_per_Core;
        int stop_Cell = start_Cell + num_per_Core;

        threads_vec.push_back(thread{&functions_library::sequence_to_string_Threads, this, start_Cell, stop_Cell, sequences, genome_Length});
    }

    if (remainder != 0)
    {
        int start_Cell = num_of_Sequences_current - remainder;
        int stop_Cell = num_of_Sequences_current;

        threads_vec.push_back(thread{&functions_library::sequence_to_string_Threads, this, start_Cell, stop_Cell, sequences, genome_Length});
    }

    for (thread &t : threads_vec)
    {
        if (t.joinable())
        {
            t.join();
        }
    }

    threads_vec.clear();

    vector<string> return_Vector = all_sequences_String;
    all_sequences_String.clear();

    return return_Vector;
}

void functions_library::sequence_to_string_Threads(int start, int stop, int **sequences, int genome_Length)
{
    vector<string> sequences_Converted;

    for (int sequence = start; sequence < stop; sequence++)
    {
        string sequence_String = "";
        for (int base = 0; base < genome_Length; base++)
        {
            sequence_String.append(to_string(sequences[sequence][base]));
        }
        sequences_Converted.push_back(sequence_String);
    }

    int index = 0;
    unique_lock<shared_mutex> ul(g_mutex);
    for (int sequence = start; sequence < stop; sequence++)
    {
        all_sequences_String[sequence] = sequences_Converted[index];
        index++;
    }
}

void functions_library::sequence_Write_Configurator(vector<string> &sequence_Write_Store_All, vector<string> &sequence_Write_Store,
                                                    int &max_sequences_per_File, const string &folder_Location, int &last_seq_Num,
                                                    vector<char> &seq_Status)
{
    for (int sequence_Collect = 0; sequence_Collect < sequence_Write_Store.size(); sequence_Collect++)
    {
        sequence_Write_Store_All.push_back(sequence_Write_Store[sequence_Collect]);
    }

    sequence_Write_Store.clear();

    if (sequence_Write_Store_All.size() >= max_sequences_per_File)
    {
        int full_Write_Count = sequence_Write_Store_All.size() / max_sequences_per_File;

        for (int full = 0; full < full_Write_Count; full++)
        {

            string fasta_file_Location = folder_Location + "/" + to_string(last_seq_Num) + "_" + to_string(last_seq_Num + max_sequences_per_File - 1) + ".nfasta";
            fstream fasta_File;
            fasta_File.open(fasta_file_Location, ios::out);

            if (fasta_File.is_open())
            {
                for (int write_Seq = (full * max_sequences_per_File); write_Seq < ((full * max_sequences_per_File) + max_sequences_per_File); write_Seq++)
                {
                    fasta_File << ">" << last_seq_Num;
                    if (seq_Status.size() == 0)
                    {
                        fasta_File << "_A";
                    }
                    else
                    {
                        fasta_File << "_" << seq_Status[write_Seq];
                    }
                    fasta_File << endl;
                    fasta_File << sequence_Write_Store_All[write_Seq] << endl;
                    last_seq_Num++;
                }

                fasta_File.close();
            }
            else
            {
                cout << "ERROR: COULD NOT CREATE NFASTA FILE: " << fasta_file_Location << endl;
                exit(-1);
            }
            // last_seq_Num = last_seq_Num + max_sequences_per_File;
        }

        // int parital_Write_Count = sequence_Write_Store_All.size() % max_sequences_per_File;
        vector<string> sequence_Write_Store_temp;
        vector<char> temp_seq_Status;
        for (int fill = full_Write_Count * max_sequences_per_File; fill < sequence_Write_Store_All.size(); fill++)
        {
            sequence_Write_Store_temp.push_back(sequence_Write_Store_All[fill]);
            if (seq_Status.size() > 0)
            {
                temp_seq_Status.push_back(seq_Status[fill]);
            }
        }

        sequence_Write_Store_All.clear();
        sequence_Write_Store_All = sequence_Write_Store_temp;

        if (seq_Status.size() > 0)
        {
            seq_Status.clear();
            seq_Status = temp_seq_Status;
        }
    }
}

void functions_library::partial_Write_Check(vector<string> &sequence_Write_Store_All,
                                            const string &folder_Location, int &last_seq_Num,
                                            vector<char> &seq_Status)
{
    if (sequence_Write_Store_All.size() > 0)
    {
        string fasta_file_Location = folder_Location + "/" + to_string(last_seq_Num) + "_" + to_string(last_seq_Num + sequence_Write_Store_All.size() - 1) + ".nfasta";
        fstream fasta_File;
        fasta_File.open(fasta_file_Location, ios::out);
        if (fasta_File.is_open())
        {
            for (int write_Seq = 0; write_Seq < sequence_Write_Store_All.size(); write_Seq++)
            {
                fasta_File << ">" << last_seq_Num;
                if (seq_Status.size() == 0)
                {
                    fasta_File << "_A";
                }
                else
                {
                    fasta_File << "_" << seq_Status[write_Seq];
                }
                fasta_File << endl;
                fasta_File << sequence_Write_Store_All[write_Seq] << endl;
                last_seq_Num++;
            }

            fasta_File.close();
        }
        else
        {
            cout << "ERROR: COULD NOT CREATE NFASTA FILE: " << fasta_file_Location << endl;
            exit(-1);
        }
        sequence_Write_Store_All.clear();
    }
}

void functions_library::sequence_Write_Configurator_transfer(vector<string> &sequence_Write_Store_All, vector<string> &sequence_Write_Store,
                                                             int &max_sequences_per_File, const string &folder_Location, int &last_seq_Num,
                                                             vector<char> &seq_Status,
                                                             string sequence_Profiles_Location, string host, string tissue, int current_Generation,
                                                             vector<int> &indexes_Written)
{
    for (int sequence_Collect = 0; sequence_Collect < sequence_Write_Store.size(); sequence_Collect++)
    {
        sequence_Write_Store_All.push_back(sequence_Write_Store[sequence_Collect]);
    }

    sequence_Write_Store.clear();

    if (sequence_Write_Store_All.size() >= max_sequences_per_File)
    {
        int full_Write_Count = sequence_Write_Store_All.size() / max_sequences_per_File;

        for (int full = 0; full < full_Write_Count; full++)
        {

            string fasta_file_Location = folder_Location + "/" + to_string(last_seq_Num) + "_" + to_string(last_seq_Num + max_sequences_per_File - 1) + ".nfasta";
            fstream fasta_File;
            fasta_File.open(fasta_file_Location, ios::out);
            fstream sequence_Profile;
            //"Sequence_ID\tHost\tTissue"
            sequence_Profile.open(sequence_Profiles_Location, ios::app);

            if (fasta_File.is_open())
            {
                for (int write_Seq = (full * max_sequences_per_File); write_Seq < ((full * max_sequences_per_File) + max_sequences_per_File); write_Seq++)
                {
                    fasta_File << ">" << last_seq_Num;
                    indexes_Written.push_back(last_seq_Num);
                    if (seq_Status.size() == 0)
                    {
                        fasta_File << "_A";
                    }
                    else
                    {
                        fasta_File << "_" << seq_Status[write_Seq];
                    }
                    fasta_File << endl;
                    fasta_File << sequence_Write_Store_All[write_Seq] << endl;
                    sequence_Profile << host << "_" << tissue << "_" << current_Generation << "_" << last_seq_Num << "\t" << host << "\t" << tissue << endl;
                    last_seq_Num++;
                }

                fasta_File.close();
                sequence_Profile.close();
            }
            else
            {
                cout << "ERROR: COULD NOT CREATE NFASTA FILE: " << fasta_file_Location << endl;
                exit(-1);
            }
            // last_seq_Num = last_seq_Num + max_sequences_per_File;
        }

        // int parital_Write_Count = sequence_Write_Store_All.size() % max_sequences_per_File;
        vector<string> sequence_Write_Store_temp;
        vector<char> temp_seq_Status;
        for (int fill = full_Write_Count * max_sequences_per_File; fill < sequence_Write_Store_All.size(); fill++)
        {
            sequence_Write_Store_temp.push_back(sequence_Write_Store_All[fill]);
            if (seq_Status.size() > 0)
            {
                temp_seq_Status.push_back(seq_Status[fill]);
            }
        }

        sequence_Write_Store_All.clear();
        sequence_Write_Store_All = sequence_Write_Store_temp;

        if (seq_Status.size() > 0)
        {
            seq_Status.clear();
            seq_Status = temp_seq_Status;
        }
    }
}

void functions_library::partial_Write_Check_transfer(vector<string> &sequence_Write_Store_All,
                                                     const string &folder_Location, int &last_seq_Num,
                                                     vector<char> &seq_Status,
                                                     string sequence_Profiles_Location, string host, string tissue, int current_Generation,
                                                     vector<int> &indexes_Written)
{
    if (sequence_Write_Store_All.size() > 0)
    {
        string fasta_file_Location = folder_Location + "/" + to_string(last_seq_Num) + "_" + to_string(last_seq_Num + sequence_Write_Store_All.size() - 1) + ".nfasta";
        fstream fasta_File;
        fasta_File.open(fasta_file_Location, ios::out);
        fstream sequence_Profile;
        sequence_Profile.open(sequence_Profiles_Location, ios::app);

        if (fasta_File.is_open())
        {
            for (int write_Seq = 0; write_Seq < sequence_Write_Store_All.size(); write_Seq++)
            {
                fasta_File << ">" << last_seq_Num;
                indexes_Written.push_back(last_seq_Num);
                if (seq_Status.size() == 0)
                {
                    fasta_File << "_A";
                }
                else
                {
                    fasta_File << "_" << seq_Status[write_Seq];
                }
                fasta_File << endl;
                fasta_File << sequence_Write_Store_All[write_Seq] << endl;
                sequence_Profile << host << "_" << tissue << "_" << current_Generation << "_" << last_seq_Num << "\t" << host << "\t" << tissue << endl;
                last_seq_Num++;
            }

            fasta_File.close();
            sequence_Profile.close();
        }
        else
        {
            cout << "ERROR: COULD NOT CREATE NFASTA FILE: " << fasta_file_Location << endl;
            exit(-1);
        }
        sequence_Write_Store_All.clear();
    }
}

void functions_library::sequence_Write_Configurator(vector<string> &sequence_Write_Store_All, vector<string> &sequence_Write_Store,
                                                    int &max_sequences_per_File, const string &folder_Location, int &last_seq_Num,
                                                    vector<char> &seq_Status,
                                                    string sequence_Profiles_Location, string host, string tissue, int current_Generation)
{
    for (int sequence_Collect = 0; sequence_Collect < sequence_Write_Store.size(); sequence_Collect++)
    {
        sequence_Write_Store_All.push_back(sequence_Write_Store[sequence_Collect]);
    }

    sequence_Write_Store.clear();

    if (sequence_Write_Store_All.size() >= max_sequences_per_File)
    {
        int full_Write_Count = sequence_Write_Store_All.size() / max_sequences_per_File;

        for (int full = 0; full < full_Write_Count; full++)
        {

            string fasta_file_Location = folder_Location + "/" + to_string(last_seq_Num) + "_" + to_string(last_seq_Num + max_sequences_per_File - 1) + ".nfasta";
            fstream fasta_File;
            fasta_File.open(fasta_file_Location, ios::out);
            fstream sequence_Profile;
            //"Sequence_ID\tHost\tTissue"
            sequence_Profile.open(sequence_Profiles_Location, ios::app);

            if (fasta_File.is_open())
            {
                for (int write_Seq = (full * max_sequences_per_File); write_Seq < ((full * max_sequences_per_File) + max_sequences_per_File); write_Seq++)
                {
                    fasta_File << ">" << last_seq_Num;
                    if (seq_Status.size() == 0)
                    {
                        fasta_File << "_A";
                    }
                    else
                    {
                        fasta_File << "_" << seq_Status[write_Seq];
                    }
                    fasta_File << endl;
                    fasta_File << sequence_Write_Store_All[write_Seq] << endl;
                    sequence_Profile << host << "_" << tissue << "_" << current_Generation << "_" << last_seq_Num << "\t" << host << "\t" << tissue << endl;
                    last_seq_Num++;
                }

                fasta_File.close();
                sequence_Profile.close();
            }
            else
            {
                cout << "ERROR: COULD NOT CREATE NFASTA FILE: " << fasta_file_Location << endl;
                exit(-1);
            }
            // last_seq_Num = last_seq_Num + max_sequences_per_File;
        }

        // int parital_Write_Count = sequence_Write_Store_All.size() % max_sequences_per_File;
        vector<string> sequence_Write_Store_temp;
        vector<char> temp_seq_Status;
        for (int fill = full_Write_Count * max_sequences_per_File; fill < sequence_Write_Store_All.size(); fill++)
        {
            sequence_Write_Store_temp.push_back(sequence_Write_Store_All[fill]);
            if (seq_Status.size() > 0)
            {
                temp_seq_Status.push_back(seq_Status[fill]);
            }
        }

        sequence_Write_Store_All.clear();
        sequence_Write_Store_All = sequence_Write_Store_temp;

        if (seq_Status.size() > 0)
        {
            seq_Status.clear();
            seq_Status = temp_seq_Status;
        }
    }
}

void functions_library::partial_Write_Check(vector<string> &sequence_Write_Store_All,
                                            const string &folder_Location, int &last_seq_Num,
                                            vector<char> &seq_Status,
                                            string sequence_Profiles_Location, string host, string tissue, int current_Generation)
{
    if (sequence_Write_Store_All.size() > 0)
    {
        string fasta_file_Location = folder_Location + "/" + to_string(last_seq_Num) + "_" + to_string(last_seq_Num + sequence_Write_Store_All.size() - 1) + ".nfasta";
        fstream fasta_File;
        fasta_File.open(fasta_file_Location, ios::out);
        fstream sequence_Profile;
        sequence_Profile.open(sequence_Profiles_Location, ios::app);
        if (fasta_File.is_open())
        {
            for (int write_Seq = 0; write_Seq < sequence_Write_Store_All.size(); write_Seq++)
            {
                fasta_File << ">" << last_seq_Num;
                if (seq_Status.size() == 0)
                {
                    fasta_File << "_A";
                }
                else
                {
                    fasta_File << "_" << seq_Status[write_Seq];
                }
                fasta_File << endl;
                fasta_File << sequence_Write_Store_All[write_Seq] << endl;
                sequence_Profile << host << "_" << tissue << "_" << current_Generation << "_" << last_seq_Num << "\t" << host << "\t" << tissue << endl;
                last_seq_Num++;
            }

            fasta_File.close();
            sequence_Profile.close();
        }
        else
        {
            cout << "ERROR: COULD NOT CREATE NFASTA FILE: " << fasta_file_Location << endl;
            exit(-1);
        }
        sequence_Write_Store_All.clear();
    }
}

void functions_library::folder_Delete(string location)
{
    try
    {
        filesystem::remove_all(location);
        cout << "Folder purged: " << location << endl;
    }
    catch (const exception &e)
    {
        cerr << "ERROR PURGING FOLDER: " << e.what() << endl;
        exit(-1);
    }
}

vector<pair<int, int>> functions_library::index_Source_folder(string &source_Target_file_Location, int &tissue_Index, int &current_Generation)
{
    cout << "Initiating indexing source folder: " << source_Target_file_Location << "/" + to_string(tissue_Index) << "/generation_" << to_string(current_Generation) << endl;

    cout << "Indexing tissue: " << tissue_Index + 1 << endl;
    vector<pair<int, int>> indexed_tissue_Folders;
    if (filesystem::exists(source_Target_file_Location + "/" + to_string(tissue_Index) + "/generation_" + to_string(current_Generation)))
    {
        for (const auto &entry : filesystem::directory_iterator(source_Target_file_Location + "/" + to_string(tissue_Index) + "/generation_" + to_string(current_Generation)))
        {
            if (entry.path().extension() == ".nfasta")
            {
                string trim_Extension = entry.path().filename().stem().string();
                vector<string> split_Data;
                split(split_Data, trim_Extension, '_');
                indexed_tissue_Folders.push_back(make_pair(stoi(split_Data[0]), stoi(split_Data[1])));
            }
        }
        sort(indexed_tissue_Folders.begin(), indexed_tissue_Folders.end());
    }
    else
    {
        cout << "ERROR: SOURCE FOLDER INDEXING FAILED. FOLDER DOES NOT EXIST\n";
        exit(-1);
    }

    return indexed_tissue_Folders;
}

vector<pair<int, int>> functions_library::index_Source_folder(string &source_Target_file_Location)
{
    cout << "Initiating indexing source folder: " << source_Target_file_Location << endl;

    // cout << "Indexing tissue: " << tissue_Index + 1 << endl;
    vector<pair<int, int>> indexed_tissue_Folders;
    if (filesystem::exists(source_Target_file_Location))
    {
        for (const auto &entry : filesystem::directory_iterator(source_Target_file_Location))
        {
            if (entry.path().extension() == ".nfasta")
            {
                string trim_Extension = entry.path().filename().stem().string();
                vector<string> split_Data;
                split(split_Data, trim_Extension, '_');
                indexed_tissue_Folders.push_back(make_pair(stoi(split_Data[0]), stoi(split_Data[1])));
            }
        }
        sort(indexed_tissue_Folders.begin(), indexed_tissue_Folders.end());
    }
    else
    {
        cout << "ERROR: SOURCE FOLDER INDEXING FAILED. FOLDER DOES NOT EXIST\n";
        exit(-1);
    }

    return indexed_tissue_Folders;
}

vector<vector<pair<int, int>>> functions_library::index_sequence_Folders(string &source_Target_file_Location, int &num_Tissues, int &current_Generation, string &multi_Read)
{
    cout << "Initiating indexing folder: " << source_Target_file_Location << endl;

    if (multi_Read == "NO")
    {
        cout << "Via single read\n";
        for (int tissue = 0; tissue < num_Tissues; tissue++)
        {
            cout << "Indexing tissue: " << tissue + 1 << endl;
            vector<pair<int, int>> indexed_tissue_Folders;
            if (filesystem::exists(source_Target_file_Location + "/" + to_string(tissue) + "/generation_" + to_string(current_Generation)))
            {
                for (const auto &entry : filesystem::directory_iterator(source_Target_file_Location + "/" + to_string(tissue) + "/generation_" + to_string(current_Generation)))
                {
                    if (entry.path().extension() == ".nfasta")
                    {
                        // filesystem::path file_Name = entry.path().filename().stem().string();
                        string trim_Extension = entry.path().filename().stem().string();
                        // cout << trim_Extension << endl;
                        vector<string> split_Data;
                        split(split_Data, trim_Extension, '_');
                        indexed_tissue_Folders.push_back(make_pair(stoi(split_Data[0]), stoi(split_Data[1])));
                    }
                }
                sort(indexed_tissue_Folders.begin(), indexed_tissue_Folders.end());
            }
            indexed_Source_Folders.push_back(indexed_tissue_Folders);
        }
    }
    else
    {
        cout << "Via multi read\n";

        for (int tissue = 0; tissue < num_Tissues; tissue++)
        {
            vector<pair<int, int>> indexed_tissue_Folders;
            indexed_Source_Folders.push_back(indexed_tissue_Folders);
        }

        int num_per_Core = num_Tissues / this->CPU_cores;
        int remainder = num_Tissues % this->CPU_cores;

        vector<thread> threads_vec;

        for (int core_ID = 0; core_ID < this->CPU_cores; core_ID++)
        {
            int start_Cell = core_ID * num_per_Core;
            int stop_Cell = start_Cell + num_per_Core;

            threads_vec.push_back(thread{&functions_library::thread_Index_sequence_Folders, this, start_Cell, stop_Cell, source_Target_file_Location, current_Generation});
        }

        if (remainder != 0)
        {
            int start_Cell = num_Tissues - remainder;
            int stop_Cell = num_Tissues;

            threads_vec.push_back(thread{&functions_library::thread_Index_sequence_Folders, this, start_Cell, stop_Cell, source_Target_file_Location, current_Generation});
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

    vector<vector<pair<int, int>>> indexed_Source_Folders_Return(indexed_Source_Folders);
    indexed_Source_Folders.clear();

    return indexed_Source_Folders_Return;
}

void functions_library::thread_Index_sequence_Folders(int start, int stop, string source_Target_file_Location, int current_Generation)
{
    vector<vector<pair<int, int>>> indexed_Source_Folders_TEMP;

    for (int tissue = start; tissue < stop; tissue++)
    {
        vector<pair<int, int>> indexed_tissue_Folders;
        if (filesystem::exists(source_Target_file_Location + "/" + to_string(tissue) + "/generation_" + to_string(current_Generation)))
        {
            for (const auto &entry : filesystem::directory_iterator(source_Target_file_Location + "/" + to_string(tissue) + "/generation_" + to_string(current_Generation)))
            {
                if (entry.path().extension() == ".nfasta")
                {
                    // filesystem::path file_Name = entry.path().filename().stem().string();
                    string trim_Extension = entry.path().filename().stem().string();
                    // cout << trim_Extension << endl;
                    vector<string> split_Data;
                    split(split_Data, trim_Extension, '_');
                    indexed_tissue_Folders.push_back(make_pair(stoi(split_Data[0]), stoi(split_Data[1])));
                }
            }
            sort(indexed_tissue_Folders.begin(), indexed_tissue_Folders.end());
        }
        indexed_Source_Folders_TEMP.push_back(indexed_tissue_Folders);
    }

    int index = 0;
    unique_lock<shared_mutex> ul(g_mutex);
    for (int tissue = start; tissue < stop; tissue++)
    {
        indexed_Source_Folders[tissue] = indexed_Source_Folders_TEMP[index];
        index++;
    }
}

vector<string> functions_library::find_Sequences_Master(string &source_Target_file_Location, vector<int> &sequence_List, int &tissue, vector<pair<int, int>> &indexed_Tissue_Folder, int &current_Generation)
{
    int num_Sequences = sequence_List.size();
    cout << "Collecting " << num_Sequences << " sequence(s)\n";
    string folder_Path = source_Target_file_Location + "/" + to_string(tissue) + "/generation_" + to_string(current_Generation);

    int num_per_Core = num_Sequences / this->CPU_cores;
    int remainder = num_Sequences % this->CPU_cores;

    vector<thread> threads_vec;

    for (int core_ID = 0; core_ID < this->CPU_cores; core_ID++)
    {
        int start_Cell = core_ID * num_per_Core;
        int stop_Cell = start_Cell + num_per_Core;

        threads_vec.push_back(thread{&functions_library::thread_find_Files, this, start_Cell, stop_Cell, sequence_List, indexed_Tissue_Folder});
    }

    if (remainder != 0)
    {
        int start_Cell = num_Sequences - remainder;
        int stop_Cell = num_Sequences;

        threads_vec.push_back(thread{&functions_library::thread_find_Files, this, start_Cell, stop_Cell, sequence_List, indexed_Tissue_Folder});
    }

    for (thread &t : threads_vec)
    {
        if (t.joinable())
        {
            t.join();
        }
    }

    threads_vec.clear();

    vector<int> Tissue_files(found_Tissue_Folder_Indexes.begin(), found_Tissue_Folder_Indexes.end());
    found_Tissue_Folder_Indexes.clear();

    cout << Tissue_files.size() << " file(s) identified\n";

    vector<pair<int, int>> sequence_FileIndex_Position_list;
    vector<string> collected_Sequences;
    for (int index = 0; index < sequence_List.size(); index++)
    {
        sequence_FileIndex_Position_list.push_back(make_pair(sequence_List[index], index));
        collected_Sequences.push_back("");
    }

    sort(sequence_FileIndex_Position_list.begin(), sequence_FileIndex_Position_list.end());

    fstream nfasta;
    int index_Files = 0;
    int line_current = 0;

    cout << "Retrieving sequence(s)\n";

    // valid_Sequences = 0;

    nfasta.open(folder_Path + "/" + to_string(indexed_Tissue_Folder[Tissue_files[index_Files]].first) + "_" + to_string(indexed_Tissue_Folder[Tissue_files[index_Files]].second) + ".nfasta", ios::in);

    for (int find = 0; find < sequence_FileIndex_Position_list.size(); find++)
    {
        // cout << "Looking for " << sequence_FileIndex_Position_list[find].first << "\n";

        while ((indexed_Tissue_Folder[Tissue_files[index_Files]].first <= sequence_FileIndex_Position_list[find].first && indexed_Tissue_Folder[Tissue_files[index_Files]].second >= sequence_FileIndex_Position_list[find].first) == 0)
        {
            nfasta.close();
            index_Files++;
            nfasta.open(folder_Path + "/" + to_string(indexed_Tissue_Folder[Tissue_files[index_Files]].first) + "_" + to_string(indexed_Tissue_Folder[Tissue_files[index_Files]].second) + ".nfasta", ios::in);
            line_current = 0;
        }

        if (nfasta.is_open())
        {
            int line_t0_check = (sequence_FileIndex_Position_list[find].first - indexed_Tissue_Folder[Tissue_files[index_Files]].first) * 2;

            string line;
            string sequence = "";

            while (getline(nfasta, line))
            {
                if (line_t0_check == line_current)
                {
                    // cout << line << endl;
                    vector<string> line_Data;
                    split(line_Data, line, '_');
                    // cout << line_Data[0].substr(1) << endl;
                    if (stoi(line_Data[0].substr(1)) == sequence_FileIndex_Position_list[find].first)
                    {
                        // if (line_Data[line_Data.size() - 1].at(0) == 'A')
                        //  {
                        getline(nfasta, line);
                        collected_Sequences[sequence_FileIndex_Position_list[find].second] = line;
                        // valid_Sequences++;
                        line_current++;
                        // }
                    }
                    else
                    {
                        cout << "ERROR: CORRECT SEQUENCE NOT FOUND AT INDEX\n";
                        cout << "Looking for: " << sequence_FileIndex_Position_list[find].first << endl
                             << "Sequence ID at location: " << line << endl
                             << "File: " << folder_Path << "/" << indexed_Tissue_Folder[Tissue_files[index_Files]].first
                             << "_" << indexed_Tissue_Folder[Tissue_files[index_Files]].second << ".nfasta" << endl;
                        exit(-1);
                    }
                    line_current++;
                    break;
                }
                line_current++;
            }
        }
        else
        {
            cout << "ERROR UNABLE TO OPEN NFATSA FILE: " << folder_Path << "/" << indexed_Tissue_Folder[Tissue_files[index_Files]].first << "_" << indexed_Tissue_Folder[Tissue_files[index_Files]].second << ".nfasta" << endl;
            exit(-1);
        }
    }
    nfasta.close();

    // cout << valid_Sequences << " live sequence(s) collected\n";

    return collected_Sequences;
}

vector<string> functions_library::find_Sequences_Master(string &source_Target_file_Location, vector<int> &sequence_List, int &tissue, vector<pair<int, int>> &indexed_Tissue_Folder, int &current_Generation, int &valid_Sequences)
{
    int num_Sequences = sequence_List.size();
    cout << "Collecting " << num_Sequences << " sequence(s)\n";
    string folder_Path = source_Target_file_Location + "/" + to_string(tissue) + "/generation_" + to_string(current_Generation);

    int num_per_Core = num_Sequences / this->CPU_cores;
    int remainder = num_Sequences % this->CPU_cores;

    vector<thread> threads_vec;

    for (int core_ID = 0; core_ID < this->CPU_cores; core_ID++)
    {
        int start_Cell = core_ID * num_per_Core;
        int stop_Cell = start_Cell + num_per_Core;

        threads_vec.push_back(thread{&functions_library::thread_find_Files, this, start_Cell, stop_Cell, sequence_List, indexed_Tissue_Folder});
    }

    if (remainder != 0)
    {
        int start_Cell = num_Sequences - remainder;
        int stop_Cell = num_Sequences;

        threads_vec.push_back(thread{&functions_library::thread_find_Files, this, start_Cell, stop_Cell, sequence_List, indexed_Tissue_Folder});
    }

    for (thread &t : threads_vec)
    {
        if (t.joinable())
        {
            t.join();
        }
    }

    threads_vec.clear();

    vector<int> Tissue_files(found_Tissue_Folder_Indexes.begin(), found_Tissue_Folder_Indexes.end());
    found_Tissue_Folder_Indexes.clear();

    cout << Tissue_files.size() << " file(s) identified\n";

    vector<pair<int, int>> sequence_FileIndex_Position_list;
    vector<string> collected_Sequences;
    for (int index = 0; index < sequence_List.size(); index++)
    {
        sequence_FileIndex_Position_list.push_back(make_pair(sequence_List[index], index));
        collected_Sequences.push_back("");
    }

    sort(sequence_FileIndex_Position_list.begin(), sequence_FileIndex_Position_list.end());

    fstream nfasta;
    int index_Files = 0;
    int line_current = 0;

    cout << "Retrieving sequence(s)\n";

    valid_Sequences = 0;

    nfasta.open(folder_Path + "/" + to_string(indexed_Tissue_Folder[Tissue_files[index_Files]].first) + "_" + to_string(indexed_Tissue_Folder[Tissue_files[index_Files]].second) + ".nfasta", ios::in);

    for (int find = 0; find < sequence_FileIndex_Position_list.size(); find++)
    {
        // cout << "Looking for " << sequence_FileIndex_Position_list[find].first << "\n";

        while ((indexed_Tissue_Folder[Tissue_files[index_Files]].first <= sequence_FileIndex_Position_list[find].first && indexed_Tissue_Folder[Tissue_files[index_Files]].second >= sequence_FileIndex_Position_list[find].first) == 0)
        {
            nfasta.close();
            index_Files++;
            nfasta.open(folder_Path + "/" + to_string(indexed_Tissue_Folder[Tissue_files[index_Files]].first) + "_" + to_string(indexed_Tissue_Folder[Tissue_files[index_Files]].second) + ".nfasta", ios::in);
            line_current = 0;
        }

        if (nfasta.is_open())
        {
            int line_t0_check = (sequence_FileIndex_Position_list[find].first - indexed_Tissue_Folder[Tissue_files[index_Files]].first) * 2;

            string line;
            string sequence = "";

            while (getline(nfasta, line))
            {
                if (line_t0_check == line_current)
                {
                    // cout << line << endl;
                    vector<string> line_Data;
                    split(line_Data, line, '_');
                    // cout << line_Data[0].substr(1) << endl;
                    if (stoi(line_Data[0].substr(1)) == sequence_FileIndex_Position_list[find].first)
                    {
                        if (line_Data[line_Data.size() - 1].at(0) == 'A')
                        {
                            getline(nfasta, line);
                            collected_Sequences[sequence_FileIndex_Position_list[find].second] = line;
                            valid_Sequences++;

                            line_current++;
                        }
                    }
                    else
                    {
                        cout << "ERROR: CORRECT SEQUENCE NOT FOUND AT INDEX\n";
                        cout << "Looking for: " << sequence_FileIndex_Position_list[find].first << endl
                             << "Sequence ID at location: " << line << endl
                             << "File: " << folder_Path << "/" << indexed_Tissue_Folder[Tissue_files[index_Files]].first
                             << "_" << indexed_Tissue_Folder[Tissue_files[index_Files]].second << ".nfasta" << endl;
                        exit(-1);
                    }
                    line_current++;
                    break;
                }
                line_current++;
            }
        }
        else
        {
            cout << "ERROR UNABLE TO OPEN NFATSA FILE: " << folder_Path << "/" << indexed_Tissue_Folder[Tissue_files[index_Files]].first << "_" << indexed_Tissue_Folder[Tissue_files[index_Files]].second << ".nfasta" << endl;
            exit(-1);
        }
    }
    nfasta.close();

    cout << valid_Sequences << " live sequence(s) collected\n";

    return collected_Sequences;
}

void functions_library::thread_find_Files(int start, int stop, vector<int> sequence_List, vector<pair<int, int>> indexed_Tissue_Folder)
{
    vector<int> caught_Indexes;
    for (int sequence = start; sequence < stop; sequence++)
    {
        int top = 0;
        int bottom = indexed_Tissue_Folder.size() - 1;
        int middle = top + ((bottom - top) / 2);

        while (top <= bottom)
        {
            if ((indexed_Tissue_Folder[middle].first <= sequence_List[sequence]) && (indexed_Tissue_Folder[middle].second >= sequence_List[sequence]))
            {
                caught_Indexes.push_back(middle);
                // found = 'Y';
                break;
            }
            else if (indexed_Tissue_Folder[middle].first < sequence_List[sequence])
            {
                top = middle + 1;
            }
            else
            {
                bottom = middle - 1;
            }
            middle = top + ((bottom - top) / 2);
        }
    }

    // int index = 0;
    unique_lock<shared_mutex> ul(g_mutex);
    for (int index = 0; index < caught_Indexes.size(); index++)
    {
        found_Tissue_Folder_Indexes.insert(caught_Indexes[index]);
        //   index++;
    }
}

float functions_library::date_to_Decimal(int year, int month, int day)
{
    cout << "Converting date (yyyy-mm-dd) to decimal\n";
    int daysInMonth[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

    int days_Total = 365;
    if (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0))
    {
        daysInMonth[1] = 29;
        days_Total = 366; // Adjust for leap year
    }

    int days_passed = 0;

    month = month - 1;

    for (int month_i = 0; month_i < month; month_i++)
    {
        days_passed = days_passed + daysInMonth[month_i];
    }

    days_passed = days_passed + day;

    float decimalDate = (float)year + ((float)days_passed / (float)days_Total);

    return decimalDate;
}

void functions_library::decimal_to_Date(float decimal_Date, int &year, int &month, int &day)
{
    cout << "Converting decimal to date (yyyy-mm-dd)\n";

    year = (int)decimal_Date;

    float decimal_Section = decimal_Date - year;

    int daysInMonth[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

    int days_Total = 365;
    if (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0))
    {
        daysInMonth[1] = 29;
        days_Total = 366; // Adjust for leap year
    }

    day = decimal_Section * days_Total;

    // int get_Month = 0;

    int days_passed = 0;
    for (int month_i = 0; month_i < 12; month_i++)
    {
        days_passed = days_passed + daysInMonth[month_i];
        if (days_passed >= day)
        {
            month = month_i + 1;
            days_passed = days_passed - daysInMonth[month_i];
            break;
        }
    }

    day = day - days_passed;

    if (day == 0)
    {
        month = month - 1;
        if (month == 0)
        {
            month = 12;
            year = year - 1;
        }

        day = daysInMonth[month - 1];
    }
}