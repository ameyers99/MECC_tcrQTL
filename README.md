# Overview
We tested HLA allomorphs and amino acid polymorhisms for association with colorectal tumour-infiltrating T-cell receptor (TCR) composition. We refer to these as TCR quantitative trait loci (tcrQTL) analyses. Further, we tested for heterogeneity by tumour subtype, including mismatch repair deficiency (MMRd), and conducted HLA colocalisation analyses with colorectal cancer risk.

## Main dataset
We analysed controlled-access germline genotype data and tumour TCR sequencing profiles from 2,750 incident, population-based colorectal cancer cases from the Molecular Epidemiology of Colorectal Cancer (MECC) study, generated using the immunoSEQ platform. Raw TCR sequencing data are available via the Adaptive Biotechnologies immuneACCESS repository: https://clients.adaptivebiotech.com/pub/schmit-2025-bmcg. 
- Reference: Schmit SL, Tsai YY, Bonner JD, Sanz-Pamplona R, Joshi AD, Ugai T, et al. Germline genetic regulation of the colorectal tumor immune microenvironment. BMC Genomics. 2024;25(1):409.

## Purpose of the analysis
- HLA variants influence sporadic colorectal cancer risk, and recent fine-mapping studies have identified HLA allomorphs and amino acid polymorphisms modifying risk in Lynch syndrome carriers, who are predisposed to immunogenic MMRd tumours.
- HLA genotypes are established quantitative trait loci (QTL) for TCR repertoire features.
- Risk-modifying HLA variants are hypothesised to influence colorectal cancer risk by altering tumour neoantigen presentation to hypervariable TCRs.
- We therefore investigated whether HLA allomorphs (2-field alleles, class II heterodimers) and amino acid polymorphisms influence tumour-infiltrating TCR composition. These analyses support emerging genetic causal inference frameworks, including ImmunoXcan (https://github.com/ameyers99/ImmunoXcan) and HLAcoloc (https://github.com/DrGBL/hlacoloc), to identify putative HLA–TCR-cancer relationships underlying disease aetiology.

## tcrQTL analysis
### Quantitative phenotypes

**Features:**
- CDR3 position-specific amino acid frequencies: 
- Joint TRBV family + CDR3 position-specific amino acid frequencies:
- TRBV family usage (primary V gene trait):
- TRBV gene usage (secondary V gene trait):
- 


- Residualised, rank normalised, standardised
- Amino acid frequency at each CDR3 position for each length, plus joint frequency of TRBV and amino acid frequency at each CDR3 position of each length
TRBV gene/family usage
Number of templates for the given feature normalised by the total number of productive templates
Linear regression
Univariate regularised (LASSO, BSLMM, GBLUP, elastic net) (reference script) and multivariate regularised regression (multivariate regression with covariance estimation, multivariate elastic net, joinet stacked elastic net, sparse PLS, sparse group lasso, multi-task lasso / L21 regularization, super learner stacking) (reference script)
GREML heritability, modified multi-allelic GRM (reference script)
V family justification


### Binary phenotypes
GREML heritability, liability threshold model (reference script)
Presence/absence
CDR3, TRBV+CDR3; and further clustered using GLIPH
Firth penalised logistic regression

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
