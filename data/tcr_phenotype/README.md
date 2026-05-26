# TCR phenotypes used as input for tcrQTL mapping analysis
## 

## Quantitative phenotypes:
### CDR3 position-specific amino acid frequencies, for receptor length 12-18
- ./cdr3_aa.txt.gz
### TRBV family frequency  
- ./trbv.txt.gz
### Joint TRBV family + CDR3 position-specific amino acid frequency, for receptor length 12-18
- ./cdr3-trbv_aa_public.txt.gz (pre-filtered for expression prevalence in ≥10% of all samples)

## Binary phenotypes
### CDR3 amino acid sequence
- ./cdr3_seq.txt.gz
### Joint TRBV family + CDR3 amino acid sequence
- ./cdr3-trbv_seq.txt.gz
### CDR3 specificity clusters (GLIPH), without conditioning clustering on TRBV family
- ./cdr3_gliph.txt.gz
### CDR3 specificity clusters (GLIPH), conditioning clustering on TRBV family
- ./cdr3-trbv_gliph.txt.gz
