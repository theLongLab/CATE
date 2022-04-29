#include "test.cuh"

test::test()
{
    cout << "Are we in GPU testing" << endl;
    int nDevices;
    cudaGetDeviceCount(&nDevices);
    for (int i = 0; i < nDevices; i++)
    {
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, i);
        cout << "GPU number\t: " << i << endl;
        cout << "GPU name\t: " << prop.name << endl;
        size_t l_free = 0;
        size_t l_Total = 0;
        cudaError_t error_id = cudaMemGetInfo(&l_free, &l_Total);
        cout << "GPU memory (GB)\t: " << l_Total / (1000 * 1000 * 1000) << endl;
        cout << "GPU number of multiprocessor(s)\t: " << prop.multiProcessorCount << endl;
        cout << "GPU number of blocks per multiprocessor\t: " << prop.maxBlocksPerMultiProcessor << endl;
        cout << "GPU threads per block\t: " << prop.maxThreadsPerBlock << endl;
        this->threads = prop.maxThreadsPerBlock;
        cout << "GPU thread(s) per multiProcessorCount\t: " << prop.maxThreadsPerMultiProcessor << endl
             << endl;
        cout << endl;
    }
    // int device;
    // cudaGetDevice(&device);
    // cout << device << endl;
}

__global__ void cuda_hello(int n, float *x, float *y)
{
    // printf("Yes\n");
    int index = threadIdx.x;
    int stride = blockDim.x;
    for (int i = index; i < n; i = i + stride)
    {
        y[i] = x[i] * y[i];
    }
}

__global__ void testing(int N, int *numbers)
{
    // printf("Yes\n");
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < N)
    {
        numbers[tid] = tid;
        // int c = 0;
        // printf("%d\n", numbers[tid]);
        tid += blockDim.x * gridDim.x;
    }
}

__global__ void testing_2(int N, int *numbers_2)
{
    // printf("Yes\n");
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < N)
    {
        // printf("%d\n", tid);
        numbers_2[tid] = tid + 100;
        // printf("%d\n", numbers_2[tid]);
        tid += blockDim.x * gridDim.x;
    }
}

__global__ void cuda_hello_2(const float *a, float *out, int arraySize)
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

float test::add(int N, float *array)
{
    //   for (int i = 0; i < N; i++)
    // {
    //     //x[i] = 1.5f;
    //     cout<<array[i]<<endl;
    // }

    float *x_Cuda, *y_Cuda;
    float *y_partial = (float *)malloc(sizeof(float));

    cudaMalloc((void **)&x_Cuda, N * sizeof(float));
    cudaMemcpy(x_Cuda, array, N * sizeof(float), cudaMemcpyHostToDevice);

    cudaMalloc((void **)&y_Cuda, sizeof(float));

    cuda_hello_2<<<1, 1024>>>(x_Cuda, y_Cuda, N);

    cudaDeviceSynchronize();
    cudaMemcpy(y_partial, y_Cuda, sizeof(float), cudaMemcpyDeviceToHost);

    // cout << y_partial[0] << endl;
    cudaFree(x_Cuda);
    cudaFree(y_Cuda);

    return y_partial[0];
}

string roundoff(float value, unsigned char prec)
{
    float pow_10 = pow(10.0f, (float)prec);
    return to_string(round(value * pow_10) / pow_10).substr(0, prec + 2);
}

void test::run()
{

    string one = "";

    // char final = new char[one.size()];

    one.append("De\tshan");
    one.append("Perera");
    one.append("Hello");
    // cout << one << endl;

    char *append;
    append = (char *)malloc((one.size() + 1) * sizeof(char));
    strcpy(append, one.c_str());

    for (int i = 0; i < one.size(); i++)
    {
        cout << append[i];
    }

    cout << "H" << endl;

    // convert[0] = line[0].data;
    // line = "Hello";
    // convert[1] = line[0];
    // line = "Mello";
    // convert[2] = line[0];
    // line = "Hello\tTello";
    // convert[3] = line[0];
    // //cout << convert[i] << endl;

    // for (size_t i = 0; i < 4; i++)
    // {
    //     //convert[i] = &line[0];
    //     cout << "i: " << i << ": \t";
    //     cout << convert[i] << endl;
    // }

    // float val = 0.175792507;

    // string round_up =roundoff(val,4);

    // cout << round_up;

    // cout << "endl" << endl;
    // int N = 2000;

    // float *x = (float *)malloc(N * sizeof(float));

    // for (int i = 0; i < N; i++)
    // {
    //     x[i] = 1.25f;
    //     //cout<<x[i]<<endl;
    // }

    // float result = add(N, x);
    // cout << "result " << result << endl;
}

void test::run_2()
{
    auto start = high_resolution_clock::now();
    int N = 1 << 20;

    float *x = new float[N];
    float *y = new float[N];

    cudaMallocManaged(&x, N * sizeof(float));
    cudaMallocManaged(&y, N * sizeof(float));

    for (int i = 0; i < N; i++)
    {
        x[i] = 1.0f;
        y[i] = 2.0f;
    }

    cuda_hello<<<1, 1024>>>(N, x, y);
    cudaDeviceSynchronize();

    cudaFree(x);
    cudaFree(y);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << duration.count() << endl;
}
