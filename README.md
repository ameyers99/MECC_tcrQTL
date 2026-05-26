# Overview
We tested associations between HLA genotypes and colorectal tumour-infiltrating T-cell receptor (TCR) composition. We refer to these as TCR quantitative trait loci (tcrQTL) analyses. We further assessed whether associations differed by tumour mismatch repair deficiency (MMRd) status, and conducted HLA colocalisation analyses with colorectal cancer risk.

## Main dataset
We analysed controlled-access germline genotype data and tumour TCR sequencing profiles from 2,750 incident, population-based colorectal cancer cases from the Molecular Epidemiology of Colorectal Cancer (MECC) study, generated using the immunoSEQ platform. Raw TCR sequencing data are available via the Adaptive Biotechnologies immuneACCESS repository: https://clients.adaptivebiotech.com/pub/schmit-2025-bmcg. 
- Reference: Schmit SL, Tsai YY, Bonner JD, Sanz-Pamplona R, Joshi AD, Ugai T, et al. Germline genetic regulation of the colorectal tumor immune microenvironment. BMC Genomics. 2024;25(1):409.

## Purpose of the analysis
- Polymorphisms in HLA genes influence sporadic colorectal cancer risk, and recent HLA fine-mapping studies have identified amino acid variants associated with risk in Lynch syndrome carriers, who are predisposed to immunogenic MMRd tumours.
- HLA genotypes are established quantitative trait loci (QTL) for circulating TCR repertoire features.
- Risk-modifying HLA variants are hypothesised to influence colorectal cancer susceptibility by altering tumour neoantigen presentation to hypervariable TCRs.
- We therefore investigated whether HLA alleles and amino acid polymorphisms regulate tumour-infiltrating TCR composition. These analyses support emerging genetic causal inference frameworks, including ImmunoXcan (https://github.com/ameyers99/ImmunoXcan) and HLAcoloc (https://github.com/DrGBL/hlacoloc), to identify putative HLA–TCR–cancer relationships stratified by MMRd status.

## Quantitative TCR phenotypes
### CDR3 amino acid composition
- ...
### TRBV usage 
- ...
### tcrQTL analysis
- ...

## Binary TCR phenotypes
### Clonotypes
- ...
### Antigen specificity clusters
- ...
### tcrQTL analysis
- ...

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
