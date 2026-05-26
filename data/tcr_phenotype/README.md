# TCR phenotypes used for tcrQTL mapping
This directory contains TCR phenotype matrices used as input for tcrQTL mapping analyses.
All files are provided in **long format** with the following columns:
- **Sample**: Unique patient identifier  
- **tcr**: TCR feature identifier (see file-specific definitions below)  
- **rate**: Proportion of productive templates within the relevant phenotype class  

## Quantitative phenotypes
### CDR3 position-specific amino acid frequencies (length 12–18)
- `./cdr3_aa.txt.gz`

CDR3 amino acid frequencies at each position, stratified by receptor length (12–18 amino acids).

**TCR ID format:** `CDR3-length-position-residue`

---

### TRBV family usage
- `./trbv.txt.gz`

TRBV family frequency across samples.

**TCR ID format:** `TRBV-family`

---

### Joint TRBV family + CDR3 position-specific amino acid frequencies (length 12–18)
- `./cdr3-trbv_aa_public.txt.gz`

Combined TRBV family and CDR3 position-specific amino acid features.

**Notes:**
- Pre-filtered for features present in **≥10% of samples**

**TCR ID format:** `TRBV-family + CDR3-length-position-residue`

---

## Binary phenotypes

### CDR3 amino acid sequences
- `./cdr3_seq.txt.gz`

Binary presence/absence of full CDR3 amino acid sequences.

---

### Joint TRBV family + CDR3 amino acid sequences
- `./cdr3-trbv_seq.txt.gz`

Binary features combining TRBV family and full CDR3 amino acid sequences.

---

### CDR3 specificity clusters (GLIPH; TRBV-agnostic)
- `./cdr3_gliph.txt.gz`

GLIPH-derived CDR3 specificity clusters without conditioning on TRBV family.

---

### CDR3 specificity clusters (GLIPH; TRBV-conditioned)
- `./cdr3-trbv_gliph.txt.gz`

GLIPH-derived specificity clusters conditioned on TRBV family usage.

---

## Notes

- All files are gzipped, tab-delimited text files.
- The **rate** column represents the fraction of productive templates within each phenotype class.
- These matrices are designed for downstream tcrQTL association analyses.
