# Overview
We tested HLA allomorphs and amino acid polymorhisms for association with colorectal tumour-infiltrating T-cell receptor (TCR) composition. We refer to these as TCR quantitative trait loci (tcrQTL) analyses. Further, we tested for heterogeneity by tumour subtype, including mismatch repair deficiency (MMRd), and conducted HLA colocalisation (HLAcoloc) analyses with colorectal cancer risk.

## Main dataset
We analysed controlled-access germline genotype data and tumour TCR sequencing profiles from 2,750 incident, population-based colorectal cancer cases from the Molecular Epidemiology of Colorectal Cancer (MECC) study, generated using the immunoSEQ platform. Raw TCR sequencing data are available via the Adaptive Biotechnologies immuneACCESS repository: https://clients.adaptivebiotech.com/pub/schmit-2025-bmcg. 
- Reference: Schmit SL, Tsai YY, Bonner JD, Sanz-Pamplona R, Joshi AD, Ugai T, et al. Germline genetic regulation of the colorectal tumor immune microenvironment. BMC Genomics. 2024;25(1):409.

## Purpose of the analysis
- HLA variants influence sporadic colorectal cancer risk, and recent fine-mapping studies have identified HLA allomorphs and amino acid polymorphisms modifying risk in Lynch syndrome carriers, who are predisposed to immunogenic MMRd tumours.
- HLA genotypes are established quantitative trait loci (QTL) for TCR repertoire features.
- Risk-modifying HLA variants are hypothesised to influence colorectal cancer risk by altering tumour neoantigen presentation to hypervariable TCRs.
- We therefore investigated whether HLA allomorphs (2-field alleles, class II heterodimers) and amino acid polymorphisms influence tumour-infiltrating TCR composition. These analyses support emerging genetic causal inference frameworks, including ImmunoXcan (https://github.com/ameyers99/ImmunoXcan) and HLAcoloc (https://github.com/DrGBL/hlacoloc), to identify putative HLA–TCR-cancer relationships for functional studies and neoantigen discovery.

## tcrQTL analysis
**See:** `./vignettes/QC/MECC_immunoSEQ_QC.ipynb` **for descriptive statistics and QC analysis of data in** `./data/tcr_phenotype/`

### To balance data sparcity and biological interpretability, we prioritised features as follows:
#### Primary analysis
- CDR3 position-specific amino acid frequencies
- CDR3 enriched motifs
#### Secondary analysis
- CDR3 amino acid sequences
- CDR3 global similarity clusters (TRBV agnostic)
- TRBV family usage
#### Exploratory/sensitivity analysis
- Joint TRBV family + CDR3 amino acid sequences
- CDR3 global similarity clusters (conditioned on TRBV family)

---

### Quantitative phenotypes
- TCR features were quantified as the number of templates with a given feature, normalised by the total productive templates per sample
- Features were residualised for covariates, rank-normalised, and standardised before analysis
- Univariate association analyses used marginal and regularised linear regression with cross-validation, including LASSO, elastic net, GBLUP, and BSLMM: `./code/tcrQTL/ImmunoXcan_compute-weights.R`
- Multivariate association analyses modelled TCR features jointly using regularised multivariate methods with cross-validation, including MRCE, multivariate elastic net, sparse PLS, sparse group LASSO, multi-task LASSO, joinet stacked elastic net, and super learner stacking: `./code/tcrQTL/ImmunoXcan_compute-weights_multivariate.R`
- Heritability was estimated using GREML with a modified multi-allelic GRM for HLA loci: `./code/tcrQTL/GRM_multiallelic.R`
- Regularised regression and multivariate omnibus likelihood ratio tests excluded the most common allele/residue at each locus as reference

**Features:**
- CDR3 position-specific amino acid frequencies: `./data/tcr_phenotype/cdr3_aa.txt.gz`
- TRBV family usage: `./data/tcr_phenotype/trbv.txt.gz`

---

### Binary phenotypes
- TCR features were quantified as above, and binarised into presence/absence traits before analysis
- Marginal association analyses used Firth-type penalised logistic regression, adjusted for covariates: `./code/tcrQTL/ImmunoXcan_logit_sumstats.R`
- Heritability was estimated using GREML and transformed to the underlying liability scale using the observed feature prevalence. We used a modified multi-allelic GRM for HLA loci: `./code/tcrQTL/GRM_multiallelic.R`
- Multivariate omnibus likelihood ratio tests excluded the most common allele/residue at each locus as reference

**Features:**
- CDR3 amino acid sequences: `./data/tcr_phenotype/cdr3_seq.txt.gz`
- Joint TRBV family + CDR3 amino acid sequences: `./data/tcr_phenotype/cdr3-trbv_seq.txt.gz`
- CDR3 global similarity clusters (TRBV agnostic): `./data/tcr_phenotype/cdr3_gliph2-global.txt.gz`
- CDR3 global similarity clusters (conditioned on TRBV family): `./data/tcr_phenotype/cdr3-trbv_gliph2-global.txt.gz`
- CDR3 enriched motifs: `./data/tcr_phenotype/cdr3_gliph2-local.txt.gz`

## System requirements
This project is designed to run in a containerised environment using Apptainer/Singularity on HPC systems.
### Container environment
All software dependencies are defined in the `env/` directory:
- `Dockerfile`: defines system libraries and base R environment
- `renv.lock`: specifies exact R package versions
The environment can be built using Docker and executed via Apptainer/Singularity on HPC systems.
### HPC runtime requirements
To execute workflows on HPC systems:
- Apptainer ≥ 1.1.8 (for running container images)
- SLURM (job scheduler)
- Optional module tools (if not using containers directly)
### Core external tools (inside container)
- GCTA (v1.94.1)
- PLINK (v1.9b_6.21)
- GEMMA (v0.98.5)
- OpenMPI (v4.1.4)
### R environment
R (v4.2.2) and all package dependencies are managed via `renv` and defined in:
- `env/renv.lock`
Restore with:
```r
renv::restore()
