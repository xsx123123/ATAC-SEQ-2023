#!/usr/bin/env bash
# author:zhang jian
# date:2023.3.7
# version:1.0v
# e-maiil:zhangjian199567@outlook.com
# this is building index scripts by bowtie2
genome=/data/genome/new_mus/ensembl_mm39/genome/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa
index_dir=/data/genome/new_mus/ensembl_mm39/index/bowtie2_index_genome
num_threads=10
# create index dir 
mkdir -p $index_dir
# use bowtie2 build index
bowtie2-build --threads $num_threads \
              -f $genome  \
              $index_dir/index
echo "index build finish!!!!!"

