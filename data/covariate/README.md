## Covariate files
- Technical ImmunoSEQ covariates: `immunoSEQ.cov`
- HLA-tcrQTL mapping covariates: `tcrQTL.cov`
  - `tcrQTL.cov` should contain at minimum: sex, age, ancestral PCs, ImmunoSEQ primer set; genotyping centre should **ONLY** be included if PCs were not computed separately for each batch (to avoid collinearity)
  - `tcrQTL.cov` covariates must be numeric coded (0/1 for binary), and one-hot encoded for multi-level (e.g., genotyping centre) with the reference group (i.e., most common genotyping centre) excluded as the reference.
  - `./mock_data/example_tcrQTL.cov` gives an example for four genotyping centres
