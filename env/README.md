## Environment specification (Dockerfile) to build a reproducible container for local development and HPC execution via Apptainer/Singularity.
### For example, to create a .sif file for running analyses with Apptainer, run the following on a machine with docker insya:
- `docker build -t tcr_qtl:1.0 --file tcr_qtl.dockerfile .`
### Followed by:
- `docker save tcr_qtl:1.0 > tcr_qtl.tar`
- `apptainer build tcr_qtl_1.0.sif docker-archive://tcr_qtl.tar`
### or simply:
- `apptainer build tcr_qtl_1.0.sif docker-daemon://tcr_qtl:1.0`
