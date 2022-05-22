library(vcfR)
library(PopGenome)
library(pegas)
GENOME.class = readVCF("/work/long_lab/deshan/1000_Genome/Neutrality/testing/chromosomes_Raw/01.vcf.gz",10000,"1",11869,14409)
GENOME.class <- neutrality.stats(GENOME.class, FAST=TRUE)
GENOME.class@Tajima.D
