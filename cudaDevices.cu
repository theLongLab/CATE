#include "cudaDevices.cuh"

cudaDevices::cudaDevices()
{
    /**
     * * Constructor Function 
     * Prints all CUDA devices.
     * No input or return passes.
     **/

    cout << "Listing all CUDA capable devices:" << endl;
    int nDevices;
    cudaGetDeviceCount(&nDevices);
    for (int i = 0; i < nDevices; i++)
    {
        /**
         * Once a CUDA enabled device has been detected its details will be printed.
         * It is advised to  check the printed information with the known details of the device to ensure that
         * NVIDIA's CUDA toolkit is properly communicating with the query device.
         **/

        cout << endl;
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, i);
        cout << "GPU number\t: " << i << endl;
        cout << "GPU name\t: " << prop.name << endl;
        size_t l_free = 0;
        size_t l_Total = 0;
        cudaError_t error_id = cudaMemGetInfo(&l_free, &l_Total);
        cout << "GPU memory (GB)\t: " << l_Total / (1000 * 1000 * 1000) << endl;
        cout << "GPU number of multiprocessor(s)\t: " << prop.multiProcessorCount << endl;
        cout << "GPU block(s) per multiprocessor\t: " << prop.maxBlocksPerMultiProcessor << endl;
        cout << "GPU thread(s) per block\t: " << prop.maxThreadsPerBlock << endl;
        cout << endl;

        /**
         * If there is an ERROR in the execution of the CUDA device the error will be printed.
         **/

        cudaError_t err = cudaGetLastError();

        if (err != cudaSuccess)
        {
            printf("CUDA Error: %s\n", cudaGetErrorString(err));
        }
    }
}