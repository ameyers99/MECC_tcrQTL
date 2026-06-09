###################################################
# Rscript: GRM_multiallelic.R
# Purpose: Estimate GRM for multi-allelic HLA loci
# Author: Aaron Meyers
# Date: 2nd Oct 2025
###################################################

# Reference 
# - Formulae described in HLA-tcrQTL mapping paper PMID:27479906

library(plink2R)
library(tidyverse)
library(optparse)

# Command options 
option_list = list(
  make_option("--bfile", action="store", default=NA, type='character',
              help="Path to PLINK binary input file prefix (minus bed/bim/fam) [required]"),
  make_option("--tmp", action="store", default=NA, type='character',
              help="Path to temporary files [required]"),
  make_option("--PATH_gcta", action="store", default="gcta64", type='character',
              help="Path to GCTA executable [%default]"),
  make_option("--variant", action="store", default="allele", type='character',
              help="Variant type for GRM estimation, either 'allele' or 'amino_acid' [default: %default]"), 
  make_option("--verbose", action="store", default=1, type="integer",
              help="How much chatter to print: 0=nothing; 1=minimal; 2=all [default: %default]")
  )

opt = parse_args(OptionParser(option_list=option_list))

if ( opt$verbose == 2 ) {
  SYS_PRINT = F
} else {
  SYS_PRINT = T
}


# Read in PLINK binary files 
dat <- read_plink(opt$bfile, impute="avg") 
rownames(dat$bed) <- sub(".*:", "", rownames(dat$bed))
n <- nrow(dat$bed)


# Store locus-level variant IDs  
if (opt$variant == "allele") { 
  loci <- unique(sub("_$", "", sub("^(([^_]+_){2}).*", "\\1", colnames(dat$bed)))) # change if different naming format
} else {
  if (opt$variant == "amino_acid") {
    loci <- unique(sub("_$", "", sub("^(([^_]+_){2}).*", "\\1", colnames(dat$bed)))) # change if different naming format
  }
}
  

# Compute GRM
GRM_list <- list()
for (prefix in loci) {
  # - Condition on locus
  cols <- grep(paste0("^", prefix), colnames(dat$bed))
  x <- dat$bed[, cols, drop = FALSE]
  m <- ncol(x)
  # - GRM 
  p_hat <- apply(x, 2, sum) / (2*n)
  w <- apply(rbind(x, p_hat), 2, function(z) (z - 2*z[length(z)]) / sqrt(2*z[length(z)]*(1 - z[length(z)])))[1:n, , drop = FALSE]
  A <- w %*% t(w) / m
  GRM_list[[prefix]] <- A
}
  # - Average over loci 
  avg_GRM <- Reduce("+", GRM_list) / length(GRM_list)


# Save lower triangle elements to grm-gz format 
# - .grm.gz file 
total_m <- ncol(dat$bed)
idx <- which(lower.tri(avg_GRM, diag = TRUE), arr.ind = TRUE)
df_GRM <- data.frame(
  ind1 = idx[, 1],
  ind2 = idx[, 2],
  m = total_m,
  GRM = avg_GRM[idx]
)
df_GRM <- df_GRM[order(df_GRM$ind1, df_GRM$ind2), ]
outfile <- paste0(opt$tmp,".grm.gz")
write.table(
  df_GRM,
  file = outfile,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)
# - .grm.id file 
dat$fam <- dat$fam[match(rownames(avg_GRM), dat$fam$V1), ]
dat$fam[,3:6] <- NULL
outfile <- paste0(opt$tmp,".grm.id")
write.table(
  dat$fam,
  file = outfile,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)


# Convert to GCTA binary format 
arg = paste(opt$PATH_gcta," --grm-gz ",opt$tmp," --make-grm-bin --out ",opt$tmp, sep='')
system(arg, ignore.stdout=SYS_PRINT, ignore.stderr=SYS_PRINT)