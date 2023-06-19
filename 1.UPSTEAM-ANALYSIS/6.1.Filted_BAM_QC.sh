# !/usr/bin/env bash
# author:zhang jian
# date:2023.3.8
# version:1.0v
# e-maiil:zhangjian199567@outlook.com
# this is final filted bam file QC scripts
# set run parameter
bed=/data/genome/new_mus/ensembl_mm39/ATAC-seq_promoters/mm39.bed.gz
genome=/data/genome/new_mus/ensembl_mm39/genome/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa
Rscript_stats=/home/jian/biosoftware/alfred-0.2.7/scripts/stats.R 
# change dir
cd $1
# create folder
mkdir -p  12.final_filted_bam_QC
# mapping sam2bam sort index & state
for i in ./11.rm_low_quality_bam/*.fin.bam
do
file_name=`basename $i .fin.bam`
# Run stats on unfiltered BAM
alfred qc -b  $bed \
          -r  $genome \
          -o ./12.final_filted_bam_QC/$file_name.bamStats.fin.filtered.tsv.gz\
          -j ./12.final_filted_bam_QC/$file_name.bamStats.fin.unfiltered.json.gz \
          ./11.rm_low_quality_bam/$file_name.fin.bam
done
# run Rscript stats.r to draw ATAC-seq QC result
for i in ./11.final_filted_bam_QC/*.bamStats.fin.filtered.tsv.gz
do
file_name=`basename $i .bamStats.fin.filtered.tsv.gz`
Rscript $Rscript_stats ./12.final_filted_bam_QC/$file_name.bamStats.fin.filtered.tsv.gz
done