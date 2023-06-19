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
mkdir -p 10.rm_more_mapping_read_bam
# use samtools remove pcr dup
for i in ./9.rm_dup_bam/*.redup.bam
do
file_name=`basename $i .redup.bam`
samtools view -@ $num_thread \
              -h \
              -b \
              -q 30 \
              ./9.rm_dup_bam/$file_name.redup.bam > ./10.rm_more_mapping_read_bam/$file_name.rmMulti.bam
done
# Remove multi-mapped reads (i.e. those with MAPQ < 30, using -q in SAMtools)

