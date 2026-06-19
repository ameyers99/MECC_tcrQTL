## Environment specification (Dockerfile) to build a reproducible container for local development and HPC execution via Apptainer/Singularity.
For example, to create a .sif file for running the following pipeline:
- `docker build -t tcr_qtl:1.0 --file tcr_qtl.dockerfile .`
- `docker save tcr_qtl:1.0 | gzip > tcr_qtl.tar.gz`
- `gunzip -c tcr_qtl.tar.gz > tcr_qtl.tar`
- `apptainer build tcr_qtl_1.0.sif docker-archive://tcr_qtl.tar`
