library(vcfR)
library(PopGenome)
library(pegas)
GENOME.class = readVCF("01.vcf.gz",10000,"1",11869,14409)
GENOME.class <- neutrality.stats(GENOME.class, FAST=TRUE)
GENOME.class@Tajima.D
