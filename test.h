#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

using namespace std;

class test
{
public:
    test();
    void run();
    void cuda_hello();
};