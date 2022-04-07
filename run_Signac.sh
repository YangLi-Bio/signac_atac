#!/usr/bin/bash
#SBATCH --account PCON0022
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=40GB
#SBATCH --gpus-per-node=1

cd /fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac_atac

module load R/4.1.0-gnu9.1
Rscript Signac_atac.R /fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac_atac/example.RDS hg38 ./ 0.4 0.25 0.01 0
