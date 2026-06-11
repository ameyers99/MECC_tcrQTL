# Script to collapse sample-level CDR3 sequence frequencies into antigen-specificity clusters with GLIPH2

rm(list = ls())

suppressPackageStartupMessages({
  library(turboGliph)
  library(data.table)
  library(tidyverse)
})

# --------------------------
# Run GLIPH2
# --------------------------

# Read CDR3 sequences
system("zcat ./MECC_tcrQTL/data/tcr_phenotype/cdr3_seq.txt.gz | awk 'NR>1 {print $2}' > ./MECC_tcrQTL/data/tcr_phenotype/cdr3_seqs_all.txt")
df <- fread(
  "./MECC_tcrQTL/data/tcr_phenotype/cdr3_seqs_all.txt",
  header = FALSE,
  nThread = 32
)
setnames(df, "V1", "CDR3b")

# Run GLIPH2
gliph2(
  df,
  result_folder = "./MECC_tcrQTL/data/tcr_phenotype/GLIPH2/CDR3/",
  refdb_beta = "gliph_reference",
  v_usage_freq = NULL,
  cdr3_length_freq = NULL,
  ref_cluster_size = "original",
  sim_depth = 1000,
  lcminp = 0.05, 
  lcminove = c(10, 10, 10),
  motif_distance_cutoff = 3,
  kmer_mindepth = 3,
  accept_sequences_with_C_F_start_end = TRUE,
  min_seq_length = 6,
  structboundaries = TRUE,
  boundary_size = 3,
  motif_length = base::c(2, 3, 4),
  discontinuous_motifs = TRUE,
  local_similarities = TRUE,
  global_similarities = TRUE,
  global_vgene = FALSE, # Swap to TRUE when conditioning on TRBV family usage
  all_aa_interchangeable = FALSE,
  boost_local_significance = TRUE,
  cluster_min_size = 3,
  hla_cutoff = 0.1,
  n_cores = 32
)

# --------------------------
# Convert to frequencies
# --------------------------

# CDR3 expression calls
pid <- fread(
  "./MECC_tcrQTL/data/tcr_phenotype/cdr3_seq.txt.gz",
  sep = " ",
  header = TRUE
)
pid <- pid %>% rename(cdr3=tcr)

# Local motif calls
local <- read.delim("./MECC_tcrQTL/data/tcr_phenotype/GLIPH2/CDR3/local_similarities_minp_0.05_minove_10_10_10_kmer_mindepth_3.txt.gz")
local$motif <- paste0(local$motif,"_",local$start,"_",local$stop)
local <- subset(local, num_in_sample >= 3) 
setDT(local)
local <- local[
  ,
  .(
    cdr3 = unique(unlist(strsplit(members, " ")))
  ),
  by = .(motif, num_in_sample, num_in_ref,
         fisher.score, num_fold, start, stop)
]
local <- merge(pid,local,by=c("cdr3"))
local <- local[,c("Sample","motif","rate")]
local <- local %>% rename(tcr=motif)
local <- local[
  ,
  .(rate = sum(rate)),
  by = .(Sample, tcr)
]
write.table(
  local,
  file = gzfile("./MECC_tcrQTL/data/tcr_phenotype/cdr3_gliph2-local.txt.gz"),
  sep = " ",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

# Global similarity calls
global <- read.delim("./MECC_tcrQTL/data/tcr_phenotype/GLIPH2/CDR3/global_similarities.txt.gz")
global <- subset(global, cluster_size >= 3)
global <- subset(global, cluster_tag != "%")
global$cluster_tag <- paste0(global$cluster_tag,"_",global$aa_at_position)
setDT(global)
global <- global[
  ,
  .(
    cdr3 = unique(unlist(strsplit(CDR3b, " ")))
  ),
  by = .(cluster_tag, cluster_size, unique_CDR3b, num_in_ref, 
         fisher.score, aa_at_position, TRBV)
]
global <- merge(pid,global,by=c("cdr3"))
global <- global[,c("Sample","cluster_tag","rate")]
global <- global %>% rename(tcr=cluster_tag)
global <- global[
  ,
  .(rate = sum(rate)),
  by = .(Sample, tcr)
]
write.table(
  global,
  file = gzfile("./MECC_tcrQTL/data/tcr_phenotype/cdr3_gliph2-global.txt.gz"),
  sep = " ",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
