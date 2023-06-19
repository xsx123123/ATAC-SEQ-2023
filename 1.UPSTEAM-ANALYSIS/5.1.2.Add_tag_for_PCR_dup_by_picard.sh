#!/usr/bin/env bash
# author:zhang jian
# date:2023.3.8
# version:1.0v
# e-maiil:zhangjian199567@outlook.com
# this is remove pcr dup read for sort bam file
# set run parameter
num_thread=20
picard=/home/jian/biosoftware/picard/build/libs/picard.jar
# change dir
cd $1
# create folder
mkdir -p  7.add_dup_tag_bam/metrics
# mapping sam2bam sort index & state
for i in ./6.rm_mt_read/*.rechrMt.bam
do
file_name=`basename $i .rechrMt.bam`
echo $file_name "is do!!!!"
# remove pcr dup by picard MarkDuplicates
java -XX:ParallelGCThreads=$num_thread \
     -jar $picard MarkDuplicates \
     QUIET=true \
     INPUT=./6.rm_mt_read/$file_name.rechrMt.bam \
     OUTPUT=./7.add_dup_tag_bam/$file_name.marked.bam \
     METRICS_FILE=./7.add_dup_tag_bam/metrics/$file_name.sorted.metrics \
     REMOVE_DUPLICATES=false \
     CREATE_INDEX=true 
echo $file_name "is DONE!!!!"
done
# --TMP_DIR <File> One or more directories with space available to be used by this program for temporary storage of working files  
# VALIDATION_STRINGENCY Validation stringency for all SAM files read by this program

