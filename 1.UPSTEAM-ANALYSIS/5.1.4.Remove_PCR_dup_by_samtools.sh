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
mkdir -p 9.rm_dup_bam
# use samtools remove pcr dup
for i in ./7.add_dup_tag_bam/*.marked.bam
do
file_name=`basename $i .marked.bam`
samtools view -@ $num_thread \
              -h \
              -b \
              -F 1024 \
              ./7.add_dup_tag_bam/$file_name.marked.bam > ./9.rm_dup_bam/$file_name.redup.bam
done
# MarkDuplicates will add a FALG 1024 to duplicate reads, we can remove them using samtools:

