![logo](https://user-images.githubusercontent.com/55466094/205796867-c55de996-aa97-415b-963c-9bcdb68a8e20.png)

# CATE (CUDA Accelerated Testing of Evolution)

A fast and scalable CUDA implementation to conduct highly parallelized evolutionary tests on large-scale genomic data.

---

```diff
- NOTE: CATE IS CURRENTLY IN IT'S ALPHA PHASE OF DEVELOPMENT
```

---
#### Description

The CATE software is a CUDA based solution to enable rapid processing of large-scale VCF files to conduct a series of six different tests on evolution.

*Currently the program is under ALPHA testing for accuracy.*

---

#### Prerequisites

1. CUDA capable hardware
2. LINUX or UNIX based kernel
3. NVIDIA's CUDA toolkit (nvcc compiler)
4. C++ compiler (gcc compiler)

---

#### How to INSTALL

To install CATE you may have to compile the code using an nvcc compiler. If so execute the following on the terminal:

Download the repository:
````
git clone "https://github.com/theLongLab/CATE/"
````
````
cd CATE/
````
*cuda 11.3.0 or higher*
````
module load cuda/11.3.0
````

Finally, compile the project:
````
nvcc -std=c++17 *.cu *.cpp -o "CATE_beta"
````
---

#### How to RUN

CATE is a command-line-based software. Its available functions include six different tests on evolution and a series of tools for editing and processing FASTA and VCF files.

The six tests on evolution are:
1. Tajima’s D
2. Fu and Li's D, D*, F, and F \*
3. Fay and Wu’s H and E
4. McDonald–Kreitman test
5. Fixation Index
6. Extended Haplotype Homozygosity

---

Currently, the program's executable is called:  
>Test_Main

To run the software you need a JSON-style parameters file. An example is provided above:

> *parameters.json*.

*The parameters file is used to specify all input and output locations as well as the gene list file locations. Each function's execution can be customized individually using the parameters file.*

The typical syntax for program __execution__ is as follows (example below shows running the Tajima's function):
> program_executable --function parameter_file

> program_executable -f parameter_file

__Example:__

>./Test_Main -t parameters.json

The __HELP__ menu will list all available functions and how each function can be executed. It can be accessed by simply typing -h as the function as shown below:

> ./Test_Main -h

---
#### How to Cite

CATE is currently being submitted for publication. However, if you halready found its framework or the software solution itself useful in your analyses, please CITE the preprint available in [bioRxiv](https://doi.org/10.1101/2023.01.31.526501). 

The details of the citation is listed below:

_Perera, D., Reisenhofer, E., Hussein, S., Higgins, E., Huber, C.D. and Long, Q.   
(2023) CATE: A fast and scalable CUDA implementation to conduct highly parallelized evolutionary tests on large scale genomic data. bioRxiv,   
10.1101/2023.01.31.526501._

---

MIT License

Copyright (c) 2022 The Long Lab

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

---
