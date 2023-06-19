#!/usr/bin/env bash
# author:zhang jian
# date:2023.6.16
# version:2.0v
# e-maiil:zhangjian199567@outlook.com
# this is an download scripts for huawei OBS data 
echo -e '\n\n'
echo "#######################download###########################################"
echo "#           this scripts only can use doenload huawei OBS date           #"
echo "#######################download###########################################"
echo -e '\n\n'
echo "pleace check the download url:$1"
echo "pleace check the key of download:$2"
echo -e '\n\n'
read -t 10 -p "Hit ENTER or wait ten seconds" 

figlet -m 60 -f big ZHANGJIANÍ
# set download dir 
save_dir=/data/data_2022
# use obsutil download data
obsutil share-cp \
           $1 \
           $save_dir \
           -ac=$2 \
           -r \
           -f \
           -u \
           -vlength \
           -vmd5 
echo "download is finish!!!!!!!"
echo "download date save in :" $save_dir