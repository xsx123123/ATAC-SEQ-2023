#!/usr/bin/env bash
# author:zhang jian
# date:2023.3.8
# version:1.0v
# e-maiil:zhangjian199567@outlook.com
# this is remove pcr dup read for sort bam file
# set run parameter
bed=/data/genome/new_mus/ensembl_mm39/ATAC-seq_promoters/mm39.bed.gz
genome=/data/genome/new_mus/ensembl_mm39/genome/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa
Rscript_stats=/home/jian/biosoftware/alfred-0.2.7/scripts/stats.R 
# change dir
cd $1
# create folder
mkdir -p  8.states_unfiltered_BAM
# mapping sam2bam sort index & state
for i in ./7.add_dup_tag_bam/*.marked.bam
do
file_name=`basename $i .marked.bam`
# Run stats on unfiltered BAM
alfred qc -b  $bed \
          -r  $genome \
          -o ./8.states_unfiltered_BAM/$file_name.bamStats.unfiltered.tsv.gz\
          -j ./8.states_unfiltered_BAM/$file_name.bamStats.unfiltered.json.gz \
          ./7.add_dup_tag_bam/$file_name.marked.bam
done
# run Rscript stats.r to draw ATAC-seq QC result
for i in ./8.states_unfiltered_BAM/*.bamStats.unfiltered.tsv.gz
do
file_name=`basename $i .bamStats.unfiltered.tsv.gz`
Rscript $Rscript_stats ./8.states_unfiltered_BAM/$file_name.bamStats.unfiltered.tsv.gz
done

