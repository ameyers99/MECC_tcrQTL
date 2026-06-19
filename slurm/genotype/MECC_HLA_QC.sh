#!/bin/bash
#SBATCH --partition=sapphire
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=30:00
#SBATCH --job-name=MECC_HLA_QC
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

# job command(s) below
PROJ_DIR="<PATH TO CLONED 'MECC_tcrQTL' PROJECT DIRECTORY>"
BFILE="<PATH TO SNP2HLA PLINK v1.9 FILE PREFIX (NO .bed/.bim/.fam)>"
APPTAINER_MODULE="<APPTAINER MODULE CALL, TESTED ON 'Apptainer/1.1.8'"

module purge
module load $APPTAINER_MODULE

cd $PROJ_DIR

LOG_DIR="$PROJ_DIR/slurm/genotype/logs"
mkdir -p "$LOG_DIR"
exec > "$LOG_DIR/MECC_HLA_QC_${SLURM_JOB_ID}.out" 2> "$LOG_DIR/MECC_HLA_QC_${SLURM_JOB_ID}.err"

Apptainer exec ./env/tcr_qtl_1.0.sif \
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
