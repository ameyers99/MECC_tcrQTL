# TCR phenotype matrices for tcrQTL mapping

Files provided in **long format** with columns:

- **Sample**: Unique patient identifier  
- **tcr**: TCR feature identifier (see file-specific definitions below)  
- **rate**: Proportion of productive templates within the relevant phenotype class. Templates were used to ensure comparable scales across batches, reducing bias from PCR amplification and sequencing depth differences. Further, template-based normalisation better reflects underlying clonotype abundance compared to reads.

**Notes:**
- Many clonotypes had unresolved TRBV genes; therefore, analyses including TRBV features were restricted to family-level resolution
- A seconday analysis was conducted with TRBV gene usage as a quantitative phenotype to evaluate HLA associations with finer-resolution repertoire structure


---

## Quantitative phenotypes

Residualised for covariates, rank-normalised, and standardised before analysis

### CDR3 position-specific amino acid frequencies
- `./cdr3_aa.txt.gz`

**Notes:**
- Stratified by CDR3 length;
- Restricted to length 12-18 amino acids (captures 95% of receptors in the cohort)

**TCR ID format:** `CDR3-length-position-residue`

---

### TRBV gene usage
- `./trbv_gene.txt.gz`

**Notes:**
- Restricted to receptors with resolved TRBV gene

---

### TRBV family usage
- `./trbv.txt.gz`

**Notes:**
- Restricted to receptors with resolved TRBV family

---

### Joint TRBV family + CDR3 position-specific amino acid frequencies
- `./cdr3-trbv_aa_public.txt.gz`

**Notes:**
- Stratified by CDR3 length;
- Restricted to length 12-18 amino acids (captures 95% of receptors in the cohort)
- Pre-filtered for features present in **≥10% of samples** before upload

**TCR ID format:** `TRBV family + CDR3-length-position-residue`

---

## Binary phenotypes
Binarised into presence/absence and filtered for ≥k% prevalence before analysis.

### CDR3 amino acid sequences
- `./cdr3_seq.txt.gz`

---

### Joint TRBV family + CDR3 amino acid sequences
- `./cdr3-trbv_seq.txt.gz`

**Notes:**
- Restricted to receptors with resolved TRBV family

---

### CDR3 specificity clusters (TRBV agnostic)
- `./cdr3_gliph.txt.gz`

**Notes:**
- Pre-filtered for public sequences present in **n≥2 samples** before clustering with GLIPH1
- Includes sequences of 'singleton clusters' in `./cdr3_seq.txt.gz`
- Restricted to receptors with CDR3 length >6 amino acids

---

### CDR3 specificity clusters (conditioned on TRBV family)
- `./cdr3-trbv_gliph.txt.gz`

**Notes:**
- Pre-filtered for public sequences present in **n≥2 samples** before clustering with GLIPH1
- Includes sequences of 'singleton clusters' in `./cdr3-trbv_seq.txt.gz`
- Restricted to receptors with resolved TRBV family and CDR3 length >6 amino acids
