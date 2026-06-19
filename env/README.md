## Environment specification (Dockerfile) to build a reproducible container for local development and HPC execution via Apptainer
### Build Docker image (local machine):
- `docker build -t tcr_qtl:1.0 --file tcr_qtl.dockerfile .`
- `docker save tcr_qtl:1.0 > tcr_qtl.tar`
### Build Apptainer image (HPC or local)
- `apptainer build tcr_qtl_1.0.sif docker-archive://tcr_qtl.tar`

---

## Slurm scripts for usage are in `./slurm`
