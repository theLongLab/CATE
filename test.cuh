#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <chrono>
#include <thread>
#include <set>
#include <mutex>
#include <shared_mutex>

using namespace std::chrono;
using namespace std;

class test
{
private:
    shared_mutex g_mutex;
    set<int> values;

public:
    int threads;
    test();
    float add(int N, float *array);
    void run();
    void run_2();
    void thread_test();
    void pointer_check(int *x);
    void add_val();
};