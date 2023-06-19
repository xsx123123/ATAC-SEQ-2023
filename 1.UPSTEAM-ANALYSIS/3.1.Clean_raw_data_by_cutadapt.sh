# !/usr/bin/env bash
# author:zhang jian
# date:2023.6.16
# version:1.0v
# e-maiil:zhangjian199567@outlook.com
# this is cutadapt scripts for fastq.gz file
# cahnge dir 
cd $1
# set run parameter 
threads=20
save_dir=$1/3.rm_adapter
# create dir 
mkdir -p $save_dir
# use adapter by cutadapt
for each in ./1.raw_data/*.R1.fastq.gz
do
file_name=`basename $each .R1.fastq.gz`
echo $file_name "is being do"
cutadapt -j $threads \
         -q 25 \
         -A CTGTCTCTTATA \
         -a CTGTCTCTTATA \
         -m 54 \
         -o ./3.rm_adapter/$file_name.R1.fastq.gz \
         -p ./3.rm_adapter/$file_name.R2.fastq.gz \
         -O 5 \
         -e 0.2 \
         ./1.raw_data/$file_name.R1.fastq.gz \
         ./1.raw_data/$file_name.R2.fastq.gz
echo $file_name "is done!!!!!"
done

# read min length = (max length read - fastqc(bad ATCG rate length)) *0.4
# read min length(54) = (max length read(150) - fastqc(bad ATCG rate length(15))) *0.4