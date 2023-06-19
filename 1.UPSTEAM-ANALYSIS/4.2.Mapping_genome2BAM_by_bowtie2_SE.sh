#!/usr/bin/env bash
# author:zhang jian
# date:2023.3.8
# version:1.0v
# e-maiil:zhangjian199567@outlook.com
# this is mapping sam2bam sort index & states scripts
# set run parameter
index=/data/genome/new_mus/ensembl_mm39/index/bowtie2_index_genome/index
num_thread=20
# change dir
cd $1
# create folder
mkdir -p  4.mapping_result/mapping_report
mkdir -p  5.bam_file/mapping_result_state
# mapping sam2bam sort index & state
for i in ./3.rm_adapter/*.fastq.gz
do
file_name=`basseame $i .fastq.gz`
echo "mapping" $file_name
# mapping by hisat2
bowtie2 -k 1 \
        -x $index \
        --mm \
        --threads $num_thread \
        --local \
        -U ./3.rm_adapter/$file_name.fastq.gz \
        -S ./4.mapping_result/$file_name.sam 2> /4.mapping_result/mapping_report/$file_name.bowtie2.log
# sam2bam by samtools
samtools sort -@ $num_thread \
              -b ./4.mapping_result/$file_name.sam \
              -o ./5.bam_file/$file_name.bam
# sort by samtools
samtools sort -@ $num_thread \
              ./5.bam_file/$file_name.bam \
              -o ./5.bam_file/$file_name.sort.bam
# index by samtools
samtools index -@ $num_thread \
               ./5.bam_file/$file_name.sort.bam
# states by samtools
samtools flagstat -@ $num_thread \
                  ./5.bam_file/$file_name.sort.bam > ./5.bam_file/mapping_result_state/$file_name.flagstat
echo $file_name "mapping sam2bam sort index & state finish!!!!!!!"
done
