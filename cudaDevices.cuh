#pragma once
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <iostream>
#include <thrust/system_error.h>
#include <thrust/system/cuda/error.h>
#include <sstream>

using namespace std;

class cudaDevices
{
    /**
     * Prints all available CUDA devices present on the current system.
     * User can use this list to determine which CUDA device to be used via the CUDA ID.
     **/

private:
public:
    /**
     *Constructor Function initiates the Function.
     * It also prints all CUDA devices.
     **/
    cudaDevices();
};