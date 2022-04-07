#!/bin/bash

cd /fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac_atac/
data_list="../../dataset_list.txt"


# Parameters to tune
reso_array=(0.4 0.8) # resolution for clustering
fc_array=(0 0.25) # cutoff of log-fold changes
p_array=(0.01 0.05) # cutoff of adjusted p-values
var_array=(0 0.25) # cutoff of covariance of Cicero


# Run Signac-ATAC under different parameter settings
for((i=0;i<2;i++));
do
    for((j=0;j<2;j++));
    do
        for((k=0;k<2;k++));
        do
            for((l=0;l<2;l++));
            do
                if [[ $i == 1 && $j == 0 && $k == 0 && $l == 0 ]] # already run
                then
                    continue
                fi
                
                reso=${reso_array[$i]}
                fc=${fc_array[$j]}
                padj=${p_array[$k]}
                covar=${var_array[$l]}
                root_dir=reso_${reso}_fc_${fc}_padj_${padj}_var_${covar}
                echo $root_dir
                mkdir $root_dir
                
                cd $root_dir
                cat $data_list | while read line
                do
			            array=(${line})
		        	    dir=${array[0]}_dir
	
		        	    job=${dir/_dir/}
		        	    #echo "Job: $job"
			    
			            rds=/fs/ess/PCON0022/liyang/STREAM/benchmarking/${dir/_dir/.RDS}
			    
			           org=${array[1]}
			           #echo "Organism: $org"
			           echo -e "#!/bin/bash\n#SBATCH --job-name=Signac_atac_${job}\n#SBATCH --time=23:50:59\n#SBATCH --output="Signac_atac_${job}.out"\n#SBATCH --account=PCON0022\n#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=20\n#SBATCH --mem=200GB\n#SBATCH --gpus-per-node=1\n\nset -e\n\nmodule load R/4.1.0-gnu9.1\n\ncd /fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac_atac/${root_dir}/\nstart=\$(date +%s)\nsleep 5;\n\nRscript /fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac_atac/Signac_atac.R $rds $org ./ ${reso} ${fc} ${padj} ${covar}\n\nend=\$(date +%s)\ntake=\$(( end - start ))\necho \${take}" > "${job}_Signac_atac.pbs"
		      	done
		      	cd ..
		      	
		    done
		 done
    done
done

# -gnu9.1
