## Core scripts to run HLA-tcrQTL mapping
Outputs .RDat files for each TCR phenotype. Depending on the pipeline, output includes:
- TCR feature prevalence
- Marginal HLA association summary statistics (standardised beta, SE, p, N): bi-allelic-encoded genotypes and optionally cis/trans DPA1-DPB1 and DQA1-DQB1 heterodimer interaction terms 
- Omnibus tests for each locus (LRT, rsq, p, N): protein positions for polymorphic residues, or HLA genes for 2-field alleles
- TCR feature prediction models (standardised weights)
- Cross-validated feature prediction performance (rsq, p)
- Multi-allelic-adapted GREML heritability (hsq, SE, p)
- Genotype coordinates and allele frequency
