## R script to prepare imputed HLA genotypes for tcrQTL mapping:
- Takes as input PLINK1 binary (.bed/.bim/.fam) files with SNP2HLA encoding;
- Assumes P/A alleles and no field separator (e.g., HLA_A_0201, HLA_DRB1_0301)
- Outputs PLINK1 binary files for (a) 2-field alleles (field-separated names) and (b) amino acid polymorphisms, with T/A (presence/absence) allele encoding irrespective of MAF – these files are used as direct input for tcrQTL analysis
- Optionally filters on MAF and call rate
