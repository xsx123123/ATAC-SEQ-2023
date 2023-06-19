#!/usr/bin/env bash
# author:zhang jian
# date:2023.3.8
# version:1.0v
# e-maiil:zhangjian199567@outlook.com
# this is an gtf convert to bed scripts for alfred software
# set run parameter
gtf_dir=/data/genome/new_mus/ensembl_mm39/annotation/Mus_musculus.GRCm39.109.chr.gtf
bed_dir=/data/genome/new_mus/ensembl_mm39/ATAC-seq_promoters/mm39.bed

# convert GTF to BED by sed & awk
sed 's/"/\t/g' $gtf_dir | \
     awk 'BEGIN{OFS=FS="\t"}{if($3=="gene") \
          {if($7=="+") {start=$4-1000; end=$4+1000;} \
          else {if($7=="-") start=$5-1000; end=$5+1000; } \
          if(start<0) start=0; print $1,start,end,"promoter";}}' > $bed_dir

# compress .bed file
pigz -p 10 $bed_dir

