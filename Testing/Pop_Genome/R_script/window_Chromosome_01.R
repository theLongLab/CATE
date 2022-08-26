library(vcfR)
library(PopGenome)
library(pegas)

GENOME.class = readVCF("01.vcf.gz",1000,"1",1,10000001)
slide.GENOME.class <- sliding.window.transform(GENOME.class,10000,10000,2)
slide.GENOME.class <- neutrality.stats(slide.GENOME.class, FAST=TRUE)
slide.GENOME.class@Tajima.D

GENOME.class = readVCF("01.vcf.gz",1000,"1",10000001,20000001)
slide.GENOME.class <- sliding.window.transform(GENOME.class,10000,10000,2)
slide.GENOME.class <- neutrality.stats(slide.GENOME.class, FAST=TRUE)
slide.GENOME.class@Tajima.D

GENOME.class = readVCF("01.vcf.gz",1000,"1",20000001,30000001)
slide.GENOME.class <- sliding.window.transform(GENOME.class,10000,10000,2)
slide.GENOME.class <- neutrality.stats(slide.GENOME.class, FAST=TRUE)
slide.GENOME.class@Tajima.D

GENOME.class = readVCF("01.vcf.gz",1000,"1",30000001,40000001)
slide.GENOME.class <- sliding.window.transform(GENOME.class,10000,10000,2)
slide.GENOME.class <- neutrality.stats(slide.GENOME.class, FAST=TRUE)
slide.GENOME.class@Tajima.D

GENOME.class = readVCF("01.vcf.gz",1000,"1",40000001,50000001)
slide.GENOME.class <- sliding.window.transform(GENOME.class,10000,10000,2)
slide.GENOME.class <- neutrality.stats(slide.GENOME.class, FAST=TRUE)
slide.GENOME.class@Tajima.D

GENOME.class = readVCF("01.vcf.gz",1000,"1",50000001,60000001)
slide.GENOME.class <- sliding.window.transform(GENOME.class,10000,10000,2)
slide.GENOME.class <- neutrality.stats(slide.GENOME.class, FAST=TRUE)
slide.GENOME.class@Tajima.D

GENOME.class = readVCF("01.vcf.gz",1000,"1",60000001,70000001)
slide.GENOME.class <- sliding.window.transform(GENOME.class,10000,10000,2)
slide.GENOME.class <- neutrality.stats(slide.GENOME.class, FAST=TRUE)
slide.GENOME.class@Tajima.D

GENOME.class = readVCF("01.vcf.gz",1000,"1",70000001,80000001)
slide.GENOME.class <- sliding.window.transform(GENOME.class,10000,10000,2)
slide.GENOME.class <- neutrality.stats(slide.GENOME.class, FAST=TRUE)
slide.GENOME.class@Tajima.D

GENOME.class = readVCF("01.vcf.gz",1000,"1",80000001,90000001)
slide.GENOME.class <- sliding.window.transform(GENOME.class,10000,10000,2)
slide.GENOME.class <- neutrality.stats(slide.GENOME.class, FAST=TRUE)
slide.GENOME.class@Tajima.D

GENOME.class = readVCF("01.vcf.gz",1000,"1",90000001,100000001)
slide.GENOME.class <- sliding.window.transform(GENOME.class,10000,10000,2)
slide.GENOME.class <- neutrality.stats(slide.GENOME.class, FAST=TRUE)
slide.GENOME.class@Tajima.D

GENOME.class = readVCF("01.vcf.gz",1000,"1",100000001,110000001)
slide.GENOME.class <- sliding.window.transform(GENOME.class,10000,10000,2)
slide.GENOME.class <- neutrality.stats(slide.GENOME.class, FAST=TRUE)
slide.GENOME.class@Tajima.D

GENOME.class = readVCF("01.vcf.gz",1000,"1",110000001,120000001)
slide.GENOME.class <- sliding.window.transform(GENOME.class,10000,10000,2)
slide.GENOME.class <- neutrality.stats(slide.GENOME.class, FAST=TRUE)
slide.GENOME.class@Tajima.D

GENOME.class = readVCF("01.vcf.gz",1000,"1",120000001,130000001)
slide.GENOME.class <- sliding.window.transform(GENOME.class,10000,10000,2)
slide.GENOME.class <- neutrality.stats(slide.GENOME.class, FAST=TRUE)
slide.GENOME.class@Tajima.D

GENOME.class = readVCF("01.vcf.gz",1000,"1",130000001,140000001)
slide.GENOME.class <- sliding.window.transform(GENOME.class,10000,10000,2)
slide.GENOME.class <- neutrality.stats(slide.GENOME.class, FAST=TRUE)
slide.GENOME.class@Tajima.D

GENOME.class = readVCF("01.vcf.gz",1000,"1",140000001,150000001)
slide.GENOME.class <- sliding.window.transform(GENOME.class,10000,10000,2)
slide.GENOME.class <- neutrality.stats(slide.GENOME.class, FAST=TRUE)
slide.GENOME.class@Tajima.D

GENOME.class = readVCF("01.vcf.gz",1000,"1",150000001,160000001)
slide.GENOME.class <- sliding.window.transform(GENOME.class,10000,10000,2)
slide.GENOME.class <- neutrality.stats(slide.GENOME.class, FAST=TRUE)
slide.GENOME.class@Tajima.D

GENOME.class = readVCF("01.vcf.gz",1000,"1",160000001,170000001)
slide.GENOME.class <- sliding.window.transform(GENOME.class,10000,10000,2)
slide.GENOME.class <- neutrality.stats(slide.GENOME.class, FAST=TRUE)
slide.GENOME.class@Tajima.D

GENOME.class = readVCF("01.vcf.gz",1000,"1",170000001,180000001)
slide.GENOME.class <- sliding.window.transform(GENOME.class,10000,10000,2)
slide.GENOME.class <- neutrality.stats(slide.GENOME.class, FAST=TRUE)
slide.GENOME.class@Tajima.D

GENOME.class = readVCF("01.vcf.gz",1000,"1",180000001,190000001)
slide.GENOME.class <- sliding.window.transform(GENOME.class,10000,10000,2)
slide.GENOME.class <- neutrality.stats(slide.GENOME.class, FAST=TRUE)
slide.GENOME.class@Tajima.D

GENOME.class = readVCF("01.vcf.gz",1000,"1",190000001,200000001)
slide.GENOME.class <- sliding.window.transform(GENOME.class,10000,10000,2)
slide.GENOME.class <- neutrality.stats(slide.GENOME.class, FAST=TRUE)
slide.GENOME.class@Tajima.D

GENOME.class = readVCF("01.vcf.gz",1000,"1",200000001,210000001)
slide.GENOME.class <- sliding.window.transform(GENOME.class,10000,10000,2)
slide.GENOME.class <- neutrality.stats(slide.GENOME.class, FAST=TRUE)
slide.GENOME.class@Tajima.D

GENOME.class = readVCF("01.vcf.gz",1000,"1",210000001,220000001)
slide.GENOME.class <- sliding.window.transform(GENOME.class,10000,10000,2)
slide.GENOME.class <- neutrality.stats(slide.GENOME.class, FAST=TRUE)
slide.GENOME.class@Tajima.D

GENOME.class = readVCF("01.vcf.gz",1000,"1",220000001,230000001)
slide.GENOME.class <- sliding.window.transform(GENOME.class,10000,10000,2)
slide.GENOME.class <- neutrality.stats(slide.GENOME.class, FAST=TRUE)
slide.GENOME.class@Tajima.D

GENOME.class = readVCF("01.vcf.gz",1000,"1",230000001,240000001)
slide.GENOME.class <- sliding.window.transform(GENOME.class,10000,10000,2)
slide.GENOME.class <- neutrality.stats(slide.GENOME.class, FAST=TRUE)
slide.GENOME.class@Tajima.D

GENOME.class = readVCF("01.vcf.gz",1000,"1",240000001,249250000)
slide.GENOME.class <- sliding.window.transform(GENOME.class,10000,10000,2)
slide.GENOME.class <- neutrality.stats(slide.GENOME.class, FAST=TRUE)
slide.GENOME.class@Tajima.D
