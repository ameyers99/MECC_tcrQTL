#!/bin/bash
#SBATCH --partition=sapphire
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=05:00:00
#SBATCH --job-name=cdr3_aa

# job command(s) below
dd="/data/gpfs/projects/punim2842/ImmunoXcan/data/MECC/immuneACCESS" # Directory with QC'd ImmuneACCESS sample-level .tsv TCRseq data

./cdr3-trbv_seq.sh 1 100 "$dd" # captures all CDR3 lengths

##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s