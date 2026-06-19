## R script to prepare imputed HLA genotypes for tcrQTL mapping
- Takes as input PLINK1 binary (.bed/.bim/.fam) files with SNP2HLA encoding;
- Assumes P/A alleles and no field separator (e.g., HLA_A_0201, HLA_DRB1_0301)
- Outputs an LD matrix and PLINK1 binary files for (a) 2-field alleles (field-separated names; e.g., HLA_A_02_01, HLA_DRB1_03_01) and (b) amino acid polymorphisms, both called using IMGT reference data `./data/genotype/HLA_DICTIONARY_AA.hg19.imgt3320.AA_tf.in_ref.rds`;
- The resulting files use T/A (presence/absence) allele encoding, irrespective of MAF, and are used as input for the tcrQTL pipelines
- Hard-calls probabilistic genotype dosages (≥0.5 posterior)
- Optionally filters on MAF and call rate
