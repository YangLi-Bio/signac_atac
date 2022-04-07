#!/bin/bash

cd /fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac_atac

data_list="../dataset_list.txt"

cat $data_list | while read line
do
    array=(${line})
    dir=${array[0]}_dir
#    echo "Directory: $dir"
    
    job=${dir/_dir/}
    echo "Job: $job"
    
#    rds=${dir/_dir/.RDS}
    
#    org=${array[1]}
#    echo "Organism: $org"
    
#    cd $dir
    sbatch ${job}_Signac.pbs
#    cd ..
done
