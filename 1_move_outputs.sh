#!/bin/bash

# This script aims to move the outputs from scratch here and group the outputs according to datasets

cd /fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac_atac

data_file=/fs/ess/PCON0022/liyang/STREAM/benchmarking/dataset_list.txt
from_dir=/fs/ess/scratch/PCON0022/liyang/stream/benchmarking/Signac_atac


# Enumerate the directories for different parameter settings
dir_ar=`ls -l -I backup $from_dir | grep -E '^d' | awk '{print $9}'`
#dir_ar=unset ${dir_ar[@]}
#dir_ar=("${dir_ar[@]:1}")

#for param_dir in `ls -l $from_dir | grep -E '^d' | awk '{print $9}'`;
#do
#    if [ "$param_dir" = "backup" ];
#    then
#        continue
#    fi
#    
#    echo $param_dir
#done


while read line;
do
    ln_ar=(${line// / })
    data=${ln_ar[0]}
    echo $data
    mkdir $data
    
    for param in $dir_ar;
    do
        #echo $param
        cp ${from_dir}/${param}/${data}_dir/regulons.qsave ${data}/${param}.qsave
        #echo ${from_dir}/${param}/${data}_dir/regulons.qsave
        #echo ${data}/${param}.qsave
    done
    
    #break
done < $data_file
