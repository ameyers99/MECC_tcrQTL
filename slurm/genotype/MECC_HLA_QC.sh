#!/bin/bash
#SBATCH --partition=sapphire
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=30:00
#SBATCH --job-name=MECC_HLA_QC

# job command(s) below
# RUN WITHIN THE 'MECC_tcrQTL' PROJECT DIRECTORY

PROJ_DIR="<PATH TO 'MECC_tcrQTL' PROJECT DIRECTORY>"
BFILE="<PATH TO OUTPUT SNP2HLA PLINK1 FILES>"

cd $PROJ_DIR

Rscript --vanilla ./code/genotype/MECC_HLA_QC.R \
--output_dir ./data/genotype \
--tmp ./data/genotype/tmp/tmp \
--bfile $BFILE \
--plink plink \
--af_filter 0.005 \
--call_filter 0.90 \
--ref ./data/genotype/HLA_DICTIONARY_AA.hg19.imgt3320.AA_tf.in_ref.rds \
--verbose 1

##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s
