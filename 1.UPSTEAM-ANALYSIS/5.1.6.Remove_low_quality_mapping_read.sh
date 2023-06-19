#!/usr/bin/env bash
# author:zhang jian
# date:2023.3.8
# version:1.0v
# e-maiil:zhangjian199567@outlook.com
# this is remove pcr dup read for sort bam file
# set run parameter
num_thread=20
# change dir 
cd $1
# create dir
mkdir -p 11.rm_low_quality_bam
# use samtools remove pcr dup
for i in ./10.rm_more_mapping_read_bam/*.rmMulti.bam
do
file_name=`basename $i .rmMulti.bam`
samtools view -@ $num_thread -F 1804 -f 2  -b ./10.rm_more_mapping_read_bam/$file_name.rmMulti.bam > ./11.rm_low_quality_bam/$file_name.fin.bam
done
# -F 1804 -f 2 can repmve low mapping quality read
# Retain properly paired reads -f 2
