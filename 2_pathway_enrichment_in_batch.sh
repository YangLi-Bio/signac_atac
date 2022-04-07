#!/bin/bash
#SBATCH --job-name=Pathway_Enrichment
#SBATCH --time=23:50:59
#SBATCH --account=PCON0022
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=50GB
#SBATCH --gpus-per-node=1

set -e

module load R/4.0.2-gnu9.1

cd /fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac_atac
start=$(date +%s)
sleep 5;

module load R/4.1.0-gnu9.1
Rscript 2_pathway_enrichment.R

end=$(date +%s)
take=$(( end - start ))
echo ${take}