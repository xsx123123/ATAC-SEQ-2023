#!/usr/bin/env bash
# author:zhang jian
# date:2023.1.7
# version:1.0v
# e-maiil:zhangjian199567@outlook.com
# this is download data form embi-ebi by acsp 

#!/usr/bin/env bash
#Name: ascp_ENA.sh
#NOTE: You can only run this script on the linux system with bash 
#Usage: bash $0  <filereport_read_run_*_tsv.txt>  <output_path>
## input
input_file=$1
output_path=$2
threads=10
usage="Usage:\nbash $0  <filereport_read_run_*_tsv.txt>  <output_path>"

## check parameters
if [ $# != 2 ];then
    echo -e $usage
    exit 1
fi

## check the ascp
###NOTE: ascp must be install && add it to the envionment variable
ascp_install='\nThe ascp command must be installed locally and add environment variables!\n\nwget -c https://d3gcli72yxqn2z.cloudfront.net/connect_latest/v4/bin/ibm-aspera-connect_4.1.1.73_linux.tar.gz\ntar -zxvf ibm-aspera-connect_4.1.1.73_linux.tar.gz\nsh ibm-aspera-connect_4.1.1.73_linux.sh\necho "PATH=$PATH:$HOME/.aspera/connect/bin/" >> $HOME/.bashrc\nsource $HOME/.bashrc\n'
ascp -h > /dev/null
if [ $? == 0 ];then
    echo "Command ascp is ready !"
else
    echo -e ${ascp_install}
    exit 1
fi

## check the input file
if [ -r $1 ];then
    echo "Input file $1 is ready!"
else
    echo "ERROR: check the input file $1 !"
    echo -e $usage
    exit 1
fi

## check the output dirctory
if [ -d $2  ];then
    echo "Output dirctory $2 is ready !"
else 
    echo "ERROR: check the output dirctory $2 !"
    echo -e $usage
    exit 1
fi

## download
sleep 1
echo "Trying to download fastq files......"

lims=0
grep -o 'ftp.*gz' ${input_file} |sed -e 's/;/\n/g' -e 's/ftp.sra.ebi.ac.uk//g' | while read fq_ftp 
do
    all=$(echo ${fq_ftp} |wc -l)
    if [ $all == 0  ];then
        echo "The $1 don't have effective ftp address. Please check it!"
    fi
    ((lims++)) 
    ascp -QT -l 300m -P33001 -i /home/jian/miniconda3/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:${fq_ftp}   ${output_path} &  #这里必须放入后台，不然多线程就是空谈。
    #Example: ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR949/SRR949627/SRR949627_1.fastq.gz ./

    if [ ${lims} -ge ${threads} ];then  #-ge 大于等于
        wait  # wait不加PID就是等待所有子进程全部结束才继续。
        lims=0  #重置限制,继续循环上面的命令
    fi
done

#bash ascp_for_ENA.sh  <filereport_read_run_*_tsv.txt>  <output path>  