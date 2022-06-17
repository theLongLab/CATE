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
private:
public:
    cudaDevices();
};