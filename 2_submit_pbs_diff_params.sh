#!/bin/bash

cd /fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac_atac

dir_list=$(ls -l ./ | awk '/^d/ {print $NF}')

for dir in $dir_list
do
    if [ $dir = "backup" ];
    then
        continue
    fi
    
    cd $dir
    #echo $dir

    for pbs in *.pbs
    do
        echo $pbs
        sbatch $pbs
        sleep 0.1s
    done
    
    cd ..
done
