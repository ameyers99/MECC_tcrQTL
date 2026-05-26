# TCR phenotypes used as input for tcrQTL mapping analysis
All files are long format with columns:
-  Sample: Unique patient ID
-  tcr: TCR ID
-  rate: Number of templates with given phenotype, divided by all productive templates within the given phenotype class

## Quantitative phenotypes:
### CDR3 position-specific amino acid frequencies, for receptor length 12-18
- ./cdr3_aa.txt.gz
- TCR ID is in the format CDR3-length-position-residue
### TRBV family frequency  
- ./trbv.txt.gz
### Joint TRBV family + CDR3 position-specific amino acid frequency, for receptor length 12-18
- ./cdr3-trbv_aa_public.txt.gz
- Pre-filtered for expression prevalence in ≥10% of all samples
- TCR ID is in the format TRBV+CDR3-length-position-residue

## Binary phenotypes
### CDR3 amino acid sequence
- ./cdr3_seq.txt.gz
### Joint TRBV family + CDR3 amino acid sequence
- ./cdr3-trbv_seq.txt.gz
### CDR3 specificity clusters (GLIPH), without conditioning clustering on TRBV family
- ./cdr3_gliph.txt.gz
### CDR3 specificity clusters (GLIPH), conditioning clustering on TRBV family
- ./cdr3-trbv_gliph.txt.gz
