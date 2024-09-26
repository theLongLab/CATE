![logo](https://user-images.githubusercontent.com/55466094/205796867-c55de996-aa97-415b-963c-9bcdb68a8e20.png)

# CATE (CUDA Accelerated Testing of Evolution)

A fast and scalable CUDA implementation to conduct highly parallelized evolutionary tests on large-scale genomic data.

---
#### Description

The CATE software is a CUDA based solution to enable rapid processing of large-scale VCF files to conduct a series of six different tests on evolution.

🔵 **Here we have provided only a brief overview of CATE's useability.**   
🟢 **Please refer to [CATE's wiki](https://github.com/theLongLab/CATE/wiki) to obtain a more detailed understanding of its functionality and usability.**

---

#### News

🔴 **CATE is currently under a major update with the integration of APOLLO.**

*Apollo is our high-performance viral epidemic simulation platform powered by CATE's architecture.*

Apollo is already available in CATE for use. Use the **--simulator** or **-sim** command. Documentation and a preprint of the simulation tool, its capabilities, and how to use Apollo are currently being worked on.

The Wiki for Apollo is currently being written and will be completed soon. 

---
#### Prerequisites

1. CUDA capable hardware
2. LINUX or UNIX based kernel
3. NVIDIA's CUDA toolkit (nvcc compiler)
4. C++ compiler (gcc compiler)

---

#### How to INSTALL

![C/C++ CUDA CI](https://github.com/theLongLab/CATE/actions/workflows/c-cpp.yml/badge.svg?event=push)

CATE can be used **on-device** via **Ananconda** or by downloading and building the **GitHub** repo. It can also be used **online** via **Google Colab**.

For the **Google Colab** notebook please follow the link to [CATE on Colab](https://colab.research.google.com/drive/1p8I2umE1U2gEB95eKwg0-fdtOLbgR13-?usp=sharing).

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/theLongLab/CATE/blob/main/Colab_CATE/CATE_on_Colab.ipynb)

To install CATE via [**Anaconda**](https://anaconda.org/deshan_CATE/cate):

[![Anaconda-Server Badge](https://anaconda.org/deshan_cate/cate/badges/version.svg)](https://anaconda.org/deshan_cate/cate)
[![Anaconda-Server Badge](https://anaconda.org/deshan_cate/cate/badges/latest_release_date.svg)](https://anaconda.org/deshan_cate/cate)
[![Anaconda-Server Badge](https://anaconda.org/deshan_cate/cate/badges/platforms.svg)](https://anaconda.org/deshan_cate/cate)

````
conda install deshan_cate::cate
````
To ensure successful installation run the following:
````
CATE -h
````
Else, if you want to install CATE **on-device** using the GitHub repo you might have to compile the code using an nvcc compiler. If so execute the following on the terminal:

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
nvcc -std=c++17 *.cu *.cpp -o "CATE"
````
To ensure successful installation try running:
````
CATE -h
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

CATE comes equipped with Apollo, our viral simulator that spans from network level to individual virion resolution complete with within-host dynamics. Apollo comes with its main simulation function and five additional utility tools.

1. _Apollo simulator_
2. Haplotype retriever
3. Pedigree retriever
4. Segregating sites matcher
5. Base substitution model to JSON
6. Recombination hotspots to JSON

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

CATE has been successfully published in the journal Methods in Ecology and Evolution (MEE). If you find this framework or the software solution useful in your analyses, please CITE the published article available in [MEE, CATE: A fast and scalable CUDA implementation to conduct highly parallelized evolutionary tests on large scale genomic data](https://doi.org/10.1111/2041-210X.14168).

To cite CATE's code please use the Zenodo release:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7987769.svg)](https://doi.org/10.5281/zenodo.7987769)

The details of the citation are listed below:

Perera, D., Reisenhofer, E., Hussein, S., Higgins, E., Huber, C. D., & Long, Q. (2023). 
CATE: A fast and scalable CUDA implementation to conduct highly parallelized evolutionary tests on large scale genomic data. 
Methods in Ecology and Evolution, 00, 1–15. 
[https://doi.org/10.1111/2041-210X.14168](https://doi.org/10.1111/2041-210X.14168).

[![Static Badge](https://img.shields.io/badge/DOI-Methods%20in%20Ecology%20and%20Evolution-%23db0f14)](https://doi.org/10.1111/2041-210X.14168)

---
#### Contact

1. For **CATE** please address your correspondence to:

**Deshan** Perera (duwagedahampriyabala@ucalgary.ca)

Dr. **Quan** Long (quan.long@ucalgary.ca)

Dr. **Christian** D. Huber (cdh5313@psu.edu)

2. For **Apollo** please address your correspondence to:

**Deshan** Perera (duwagedahampriyabala@ucalgary.ca)

Dr. **Quan** Long (quan.long@ucalgary.ca)

Dr. **Alexander** Platt (alexander.platt@pennmedicine.upenn.edu)

---
#### Software License

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
