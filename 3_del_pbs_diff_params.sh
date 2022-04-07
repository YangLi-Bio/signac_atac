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
    echo $dir

    rm Signac_*.out
    rm -r *_dir/
    
    cd ..
done
