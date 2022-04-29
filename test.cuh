#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <chrono>
using namespace std::chrono;
using namespace std;

class test
{
public:
    int threads;
    test();
    float add(int N, float *array);
    void run();
    void run_2();
};