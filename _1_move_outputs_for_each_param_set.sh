#!/bin/bash

cd /fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac_atac
from_dir=/fs/ess/scratch/PCON0022/liyang/stream/benchmarking/Signac_atac

for dir in `ls -l $from_dir | grep -E '^d' | awk '{print $9}'`;
do
    if [[ $dir == backup ]]
    then
        continue
    fi
    
    mkdir $dir
    echo "----> ${dir}"
    for data in  `ls -l ${from_dir}/${dir} | grep -E '^d' | awk '{print $9}'`;
    do
        echo  "--------> ${dir}/${data//_dir/.qsave}"
        cp ${from_dir}/${dir}/${data}/regulons.qsave ${dir}/${data//_dir/.qsave}
    done
    #break
done
