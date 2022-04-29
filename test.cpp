#include "test.h"

test::test()
{
    cout << "We are in GPU testing 3" << endl;
    int nDevices;

    cudaGetDeviceCount(&nDevices);
    for (int i = 0; i < nDevices; i++)
    {
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, i);
        cout << "VGA number\t: " << i << endl;
        cout << "VGA name\t: " << prop.name << endl;
    }
}

__global__ void test::cuda_hello()
{

}

void test::run()
{
    cuda_hello<<<1,1>>>();
}
