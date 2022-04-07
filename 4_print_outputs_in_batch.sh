#!/bin/bash

cd /fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac

data_list="../dataset_list.txt"
output_file="Signac_batch_output.txt"

rm $output_file
touch $output_file

cat $data_list | while read line
do
    array=(${line})
    dir=${array[0]}_dir
    echo "Directory: $dir"
#    cd $dir
    
    echo $dir >> ${output_file}
    tail -n 5 Signac_${array[0]}*.out >> ${output_file}
    echo "" >> ${output_file}
    echo "" >> ${output_file}
    echo "" >> ${output_file}

#    cd ..
done
