# TCR phenotype matrices for tcrQTL mapping

Files provided in **long format** with columns:

- **Sample**: Unique patient identifier  
- **tcr**: TCR feature identifier (see file-specific definitions below)  
- **rate**: Proportion of productive templates within the relevant phenotype class. Templates were used to ensure comparable scales across batches, reducing bias from PCR amplification and sequencing depth differences. Further, template-based normalisation better reflects underlying clonotype abundance compared to reads.

**Notes:**
- Many clonotypes had unresolved TRBV genes; therefore, analyses including TRBV features were restricted to family-level resolution
- **rate** values for GLIPH2 antigen specificity clusters are computed as the sum of rate values for TCR sequences contributing to the cluster


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

### TRBV family usage
- `./trbv.txt.gz`

**Notes:**
- Restricted to receptors with resolved TRBV family
- Total number of productive templates used to normalise feature frequency, irrespecitve of TRBV resolution, to minimise sample-specific biases in sequencing quality

---

## Binary phenotypes
Binarised into presence/absence and filtered for ≥1% prevalence before analysis.

### CDR3 amino acid sequences
- `./cdr3_seq.txt.gz`

---

### Joint TRBV family + CDR3 amino acid sequences
- `./cdr3-trbv_seq.txt.gz`

**Notes:**
- Restricted to receptors with resolved TRBV family
- Total number of productive templates used to normalise feature frequency, irrespecitve of TRBV resolution, to minimise sample-specific biases in sequencing quality

---

### CDR3 global similarity clusters (TRBV agnostic)
- `./cdr3_gliph2-global.txt.gz`

**Notes:**
- Restricted to receptors with CDR3 length >6 amino acids
- Sequences with a Hamming distance of 1 whose different amino acid has a BLOSUM62 score ≥0 are evaluated as global similarities
- Filtered for clusters with ≥3 contributing TCRs

---

### CDR3 global similarity clusters (conditioned on TRBV family)
- `./cdr3-trbv_gliph2-global.txt.gz`

**Notes:**
- Restricted to receptors with resolved TRBV family and CDR3 length >6 amino acids
- Filtered for clusters with ≥3 contributing TCRs
- Sequences with a Hamming distance of 1 whose different amino acid has a BLOSUM62 score ≥0 are evaluated as global similarities

---

### CDR3 local motifs (TRBV agnostic)
- `./cdr3_gliph2-local.txt.gz`

**Notes:**
- Restricted to receptors with CDR3 length >6 amino acids
- Restricted 2mer, 3mer, and 4mers, including discontinuous motifs 
- Filtered for motifs with ≥3 contributing TCRs, as well as Fisher score <0.05 and OvE ≥10-fold enrichment compared to GLIPH naive reference repertoire
