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
pheno="cdr3-trbv_aa.freq" 

# Create site-length-specific CDR3 amino acid frequency call files
for length in {12..18}; do # Covers 95% of all productive CDR3 sequences in MECC
    for pos in $(seq 1 $length); do
        ./cdr3-trbv_aa.sh "$length" "$pos" "$dd"
    done
done

# Index CDR3 amino acid phenotype names with site and length
cd $dd/$pheno/
for d in */; do
    d=${d%/}
    f="$d/cdr3-trbv_aa.txt.gz"
    zcat "$f" | awk -v prefix="$d" '
    BEGIN{OFS=" "}
    NR==1 {print; next}
    {
        sub(/\+/, "+" prefix "_", $2)
        print
    }' | gzip > "${f%.txt.gz}.tmp.gz"
    mv "${f%.txt.gz}.tmp.gz" "$f"
done

# Concatonate all site-length-specific CDR3 amino acid frequency call files
cd $dd/$pheno/
set -euo pipefail
OUTDIR="combined" # Create output directory
mkdir -p "${OUTDIR}"
OUTFILE="${OUTDIR}/cdr3-trbv_aa.txt" # Output file
first=1
for dir in L*_P*; do # Loop through all subdirectories
    [ -d "$dir" ] || continue # Skip if not a directory
    file=""
    if [ -f "${dir}/cdr3-trbv_aa.txt.gz" ]; then
        file="${dir}/cdr3-trbv_aa.txt.gz"
        reader="zcat"
    else
        echo "No cdr3-trbv_aa file found in ${dir}"
        continue
    fi
    echo "Processing ${file}"
    if [ $first -eq 1 ]; then
        $reader "$file" > "$OUTFILE" # Keep header from first file
        first=0
    else
        $reader "$file" | tail -n +2 >> "$OUTFILE" # Skip header for subsequent files
    fi
done
gzip -f "$OUTFILE"
echo "Finished:"
echo "${OUTFILE}.gz"

##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s