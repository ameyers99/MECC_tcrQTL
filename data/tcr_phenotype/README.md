# TCR phenotype matrices for tcrQTL mapping

Files provided in **long format** with the columns:

- **Sample**: Unique patient identifier  
- **tcr**: TCR feature identifier (see file-specific definitions below)  
- **rate**: Proportion of productive templates within the relevant phenotype class  

---

## Quantitative phenotypes

### CDR3 position-specific amino acid frequencies (length 12–18)
- `./cdr3_aa.txt.gz`

CDR3 amino acid frequencies at each position, stratified by receptor length (12–18 amino acids).

**TCR ID format:** `CDR3-length-position-residue`

---

### TRBV family usage
- `./trbv.txt.gz`

TRBV family frequency across samples.

---

### Joint TRBV family + CDR3 position-specific amino acid frequencies (length 12–18)
- `./cdr3-trbv_aa_public.txt.gz`

Combined TRBV family and CDR3 position-specific amino acid features.

**Notes:**
- Pre-filtered for features present in **≥10% of samples**

**TCR ID format:** `TRBV family + CDR3-length-position-residue`

---

## Binary phenotypes
These phenotypes were binarised into presence/absence and filtered for ≥5% prevalence before analysis.

### CDR3 amino acid sequences
- `./cdr3_seq.txt.gz`

---

### Joint TRBV family + CDR3 amino acid sequences
- `./cdr3-trbv_seq.txt.gz`

---

### CDR3 specificity clusters (GLIPH; TRBV agnostic)
- `./cdr3_gliph.txt.gz`

---

### CDR3 specificity clusters (GLIPH; conditioned on TRBV family)
- `./cdr3-trbv_gliph.txt.gz`
