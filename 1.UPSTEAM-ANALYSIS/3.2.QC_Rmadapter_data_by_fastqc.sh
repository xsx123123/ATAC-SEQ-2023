#!/usr/bin/env bash
# author:zhang jian
# date:2023.6.16
# version:1.0v
# e-maiil:zhangjian199567@outlook.com
# this is fastqc scripts for remove adapter clean data
# set run parameter
threads=20
# change dir
cd $1 
# create dir
mkdir -p  2.fastqc_report/remove_adapter_data_fastqc_report
# use fastqc qc  fastq.gz  file
fastqc -t $threads -q -o 2.fastqc_report/remove_adapter_data_fastqc_report \
       ./3.rm_adapter/*fastq.gz
# meige fastqc report
multiqc ./2.fastqc_report/remove_adapter_data_fastqc_report \
        -o ./2.fastqc_report/remove_adapter_data_fastqc_report
# how to use this scripts
# 2.fastqc_report.sh <fastq.gz file path> 
