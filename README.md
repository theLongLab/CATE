# CATE (CUDA Accelerated Testing of Evolution)

A fast and scalable CUDA implementation to conduct highly parallelized evolutionary tests on large scale genomic data.

---

#### <span style="color:red">NOTE: BETA TESTING MODE</span>

---
#### Description

The C.A.T.E (CUDA Accelerated Testing of Evolution) software is CUDA based solution to enable rapid processing of large scale VCF files to conduct a series of six different tests on evolution.

*Currently the program is under BETA testing for accuracy.*

---

#### Prerequisites

1. CUDA capable hardware
2. LINUX or UNIX based kernal
3. NVIDIA's CUDA toolkit (nvcc compiler)
4. C++ compiler (gcc compiler)

---

#### How to RUN

CATE is a command line based software. It's available functions include six different tests on evolution and a series of tools for editing and processing FASTA and VCF files.

The six tests on evolution are:
1. Tajima’s D
2. Fu and Li's D, D*, F and F \*
3. Fay and Wu’s H and E
4. McDonald–Kreitman test
5. Fixation Index
6. Extended Haplotype Homozygosity
