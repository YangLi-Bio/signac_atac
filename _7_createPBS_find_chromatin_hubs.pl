#!/usr/bin/perl

use strict;
use warnings;

my $datalist = "/fs/ess/PCON0022/liyang/STREAM/benchmarking/dataset_list.txt";

open (IN, $datalist);
while (<IN>) {
    chomp;
    my @ln_array = split(/\s+/, $_);
    my $data = $ln_array[0];
    my $org = $ln_array[1];
    my $dir = "/fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac_atac\/$data";
    my @dir_list = </fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac_atac/$data/*.qsave>;
    my $rds = "/fs/ess/PCON0022/liyang/STREAM/benchmarking/".${data}.".RDS";
    
    if ($#dir_list < 1) {
        next;
    }
    
    print $#dir_list + 1, " files in dir: ", $dir, "\n";
    #chdir $dir;
    
    foreach my $param(@dir_list)
    {
        #print $param, "\n";
        my @param_ar = split(/\//, $param);
        my $job = "$param_ar[$#param_ar - 1]\-$param_ar[$#param_ar]\.Signac_find_hubs";
        $job =~ s/\.links\.qsave//g;
        
        my $printfile = $job;
        $printfile = $job."\.out";
        
        pop @param_ar;
        pop @param_ar;
        push @param_ar, $job;
        #my $pbs = join("\.", @param_ar);
        #$pbs .= ".kegg.pbs";
        my $pbs = "find_chromatin_hubs_tasks\/".$job."\.pbs";

        print $pbs, "\n";
        
        open(PBS, ">", $pbs) or die "Error: failed in openning ", $job, "\n";
        print PBS "#!/bin/bash\n";
        print PBS "#SBATCH --time=03:50:59\n";
        print PBS "#SBATCH --output=", $printfile, "\n";
        print PBS "#SBATCH --account=PCON0022\n";
        print PBS "#SBATCH --nodes=1\n";
        print PBS "#SBATCH --ntasks-per-node=2\n";
        print PBS "#SBATCH --mem=16GB\n";
        print PBS "#SBATCH --gpus-per-node=1\n\n";
        print PBS "set \-e\n\n";
        print PBS "module load R/4.1.0-gnu9.1\n\n";
        print PBS "cd /fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac_atac/\n";
        print PBS "Rscript /fs/ess/PCON0022/liyang/STREAM/comparison/2_chrom_hubs/2_find_chrom_hub.R $param\n\n";
        close PBS;
    }
    
    #chdir "/fs/ess/PCON0022/liyang/Signac_atac/benchmarking/Signac_joint_regulons/";
}
close IN;
