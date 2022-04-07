#!/bin/bash

# This script aims to run pathway enrichment for each directory

cd /fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac_atac
data_file=/fs/ess/PCON0022/liyang/STREAM/benchmarking/dataset_list.txt


cat $data_file | while read line
do
    ln_ar=(${line// / })
    data=${ln_ar[0]}
    org=${ln_ar[1]}
    #echo "#!/bin/bash"
    echo '#!/bin/bash' > "${data}_Signac_enrichr.pbs"
    echo -e "#SBATCH --job-name=Signac_enrichr_${data}\n#SBATCH --time=11:50:59\n#SBATCH --output=Signac_enrichr_${data}.out\n#SBATCH --account=PCON0022\n#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=4\n#SBATCH --mem=40GB\n#SBATCH --gpus-per-node=1\n\nset -e\n\nmodule load R/4.1.0-gnu9.1\n\ncd /fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac_atac/\nstart=\$(date +%s)\nsleep 5;\n\nRscript 2_pathway_enrichment_per_dir.R /fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac_atac/${data} $org\n\nend=\$(date +%s)\ntake=\$(( end - start ))\necho \${take}" >> "${data}_Signac_enrichr.pbs"
done
