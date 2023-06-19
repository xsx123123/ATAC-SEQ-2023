#!/usr/bin/env bash
# author:zhang jian
# date:2023.3.8
# version:1.0v
# e-maiil:zhangjian199567@outlook.com
# this is remove mapping to MT chr read
# set run parameter
MT_ID=MT
thread=20
#change dir
cd $1
# create dir
mkdir -p  6.rm_mt_read/remove_result
# use samtools grep remove MT read
for i in ./5.bam_file/*.sort.bam
do
file_name=`basename $i .sort.bam`
samtools view -@ $thread \
              -h ./5.bam_file/$file_name.sort.bam | \
               grep -v $MT_ID | samtools sort -@ $thread \
               -O bam -o ./6.rm_mt_read/$file_name.rechrMt.bam
# states by samtools
samtools flagstat -@ $thread ./6.rm_mt_read/$file_name.rechrMt.bam > ./6.rm_mt_read/remove_result/$file_name.rechrMt.flagstat
done

# -h, --with-header          Include header in SAM output
# this is mapping result sample, you should certain mito chr NAME "MT" or "chrM"
# A00838:783:HFVNHDSX5:1:1472:28736:7983	161	1	3023671	1	26M	MT	6088	0	
# A00838:783:HFVNHDSX5:1:1356:1380:26991	417	1	3023683	1	24M	MT	6118	0	
# A00838:783:HFVNHDSX5:1:2615:10538:31328	433	1	4428526	1	29M	MT	6090	0	


