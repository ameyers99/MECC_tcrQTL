#######################################################
# Rscript: MECC_HLA_QC.R
# Purpose: Call HLA amino acid haplotypes and prepare \
#          SNP2HLA genotypes for HLA-tcrQTL mapping
# Author: Aaron Meyers
# Date: 05 April 2026
#######################################################

set.seed(1999)

# NOTE
# - Assumes SNP2HLA allele coding (P/A) and naming (e.g., HLA_B_1402)

# --------------------------------------------------------
# SET UP ENVIRONMENT
# --------------------------------------------------------

# Load packages
suppressMessages({
  suppressWarnings({
    library("optparse", quietly=TRUE, warn.conflicts=FALSE)
    library("plink2R", quietly=TRUE, warn.conflicts=FALSE)
    library("tidyverse", quietly=TRUE, warn.conflicts=FALSE)
    library("data.table", quietly=TRUE, warn.conflicts=FALSE)
  })
})

# Command options
option_list = list(
  make_option("--output_dir", action="store", default=NA, type='character', 
              help="Path to output directory"), 
  make_option("--bfile", action="store", default=NA, type='character',
              help="Path to PLINK binary input file prefix (minus bed/bim/fam)"),
  make_option("--tmp", action="store", default=NA, type='character',
              help="Path to temporary files"),
  make_option("--plink", action="store", default="plink", type='character',
              help="Path to plink executable [%default]"),         
  make_option("--af_filter", action="store", default=0.005, type='double', 
              help="Minor allele frequency threshold [default: %default]"), 
  make_option("--call_filter", action="store", default=0.90, type='double', 
              help="Genotype call rate threshold [default: %default]"), 
  make_option("--ref", action="store", default=NA, type='character',
              help="Path to HLA info data HLA_DICTIONARY_AA.hg19.imgt3320.AA_tf.in_ref.rds"),
  make_option("--noclean", action="store_true", default=FALSE,
              help="Do not delete any temporary files (for debugging) [default: %default]"),
  make_option("--verbose", action="store", default=1, type="integer",
              help="How much chatter to print: 0=nothing; 1=minimal; 2=all [default: %default]")
)
opt = parse_args(OptionParser(option_list=option_list))

if ( opt$verbose == 2 ) {
  SYS_PRINT = F
 } else {
  SYS_PRINT = T
 }

# --------------------------------------------------------
# CLEANUP FUNCTION
# --------------------------------------------------------

cleanup = function() {
  if ( !opt$noclean ) {
    parent_dir <- dirname(opt$tmp)
    arg <- paste("rm -rf", parent_dir)
    system(arg)
  }
}

# --------------------------------------------------------
# IMPORT HLA ALLELE CALLS
# --------------------------------------------------------

# Make temporary directory for this run
tmp_dir <- dirname(opt$tmp)
arg <- paste("mkdir -p", tmp_dir)
system(arg)
   
# Extract HLA allele calls
arg <- paste0("awk '$2 ~ /^HLA/ {print $2}' ", opt$bfile, ".bim > ", opt$tmp,".txt")
system(arg)
arg <- paste0(opt$plink, " --bfile ", opt$bfile, " --extract ", opt$tmp,".txt", " --make-bed --out ", opt$tmp,"1")
system(arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)

# Force T/A encoding (SNP2HLA presence/absence, P/A) irrespectie of MAF
arg <- paste0("awk '$2 ~ /^HLA_/ { ",
  "a1 = ($5==\"P\" ? \"T\" : ($5==\"A\" ? \"A\" : \"\")); ",
  "a2 = ($6==\"P\" ? \"T\" : ($6==\"A\" ? \"A\" : \"\")); ",
  "print $2 \"\\t\" $5 \"\\t\" $6 \"\\t\" a1 \"\\t\" a2 ",
  "}' ",
  opt$tmp,"1.bim > ", opt$tmp,".txt")
system(arg)
arg <- paste0(opt$plink, " --bfile ", opt$tmp,"1", " --update-alleles ", opt$tmp,".txt", " --make-bed --out ", opt$tmp,"2")
system(arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)

arg <- paste0("awk 'BEGIN{OFS=\"\\t\"} {print $0, \"T\"}' ",opt$tmp, "2.bim > ", opt$tmp, ".txt")
system(arg)
arg <- paste0(opt$plink, " --bfile ", opt$tmp,"2", " --a1-allele ", opt$tmp,".txt 7 2", " --make-bed --out ", opt$tmp,"3")
system(arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)

# Read genotypes
genos <- paste0(opt$tmp,"3")
genos = read_plink(genos,impute="none")
if (opt$verbose > 0) {
  print(paste0("Number of HLA alleles: ",nrow(genos$bim)))
 }

# --------------------------------------------------------
# MAP AND RETAIN 2-FIELD (4-DIGIT) ALLELE CODES
# --------------------------------------------------------

# Plausible field separators
x <- genos$bim$V2
prefix <- sub("^(([^_]*_){2}).*", "\\1", x)
rem <- sub("^([^_]*_){2}", "", x)
genos$bim$V2.1 <- paste0(
  prefix,
  substr(rem, 1, 2), "_",
  substr(rem, 3, nchar(rem))
)
genos$bim$V2.2 <- paste0(
  prefix,
  substr(rem, 1, 3), "_",
  substr(rem, 4, nchar(rem)))
genos$bim <- genos$bim[,c("V2","V2.1","V2.2")]

# Referent HLA info
info <- readRDS(opt$ref)

# Map 2-field allele codes to reference
ref_vec <- unique(info$hla)
m1 <- genos$bim$V2.1 %in% ref_vec
m2 <- genos$bim$V2.2 %in% ref_vec
genos$bim$hla <- NA_character_
genos$bim$hla[m1 & !m2] <- genos$bim$V2.1[m1 & !m2]
genos$bim$hla[m2 & !m1] <- genos$bim$V2.2[m2 & !m1]
ambiguous <- which(m1 & m2)
if (length(ambiguous) > 0) {
  for (i in ambiguous) {
    message(sprintf("%s 2-field could not be resolved", genos$bim$V2[i]))
  }
  genos$bim$hla[ambiguous] <- NA_character_
}

# Update resolved 2-field allele names
map <- genos$bim[!is.na(genos$bim$hla), c("V2", "hla")]
map <- map[!duplicated(map$V2), ]
rename_vec <- setNames(map$hla, map$V2)
colnames(genos$bed) <- ifelse(
  colnames(genos$bed) %in% names(rename_vec),
  rename_vec[colnames(genos$bed)],
  colnames(genos$bed)
)

# Keep resolved 2-field alleles
hla_resolved <- unique(na.omit(genos$bim$hla))
keep_cols <- colnames(genos$bed) %in% hla_resolved
genos$bed <- genos$bed[, keep_cols, drop = FALSE]

# --------------------------------------------------------
# CALL AMINO ACID POLYMORPHISMS
# --------------------------------------------------------

info$pos <- gsub("-", ".", info$pos)
info$var <- paste0(info$gene,"_",info$pos,"_",info$AA)
allvar <- sort( unique(info$var) )
allelenames <- colnames(genos$bed)
adopted<-NULL
geno_df <- as.data.frame(genos$bed)
for(thisvar in allvar){
  hlas<-as.character(subset(info,var==thisvar)$hla)
  hlas<-hlas[hlas %in% allelenames]
  if(length(hlas)>0){
    geno_df$thisvar <- rowSums(geno_df[hlas])
    colnames(geno_df)[ncol(geno_df)] <- thisvar
    adopted<-c(adopted,thisvar)
  }
}

# --------------------------------------------------------
# CALL RATE AND MAF FILTER
# --------------------------------------------------------

# Hard-call genotypes 
geno_df[] <- lapply(geno_df, function(x) round(as.numeric(x)))

# Call rate filter
if (!is.na(opt$call_filter)) {
  vars <- colnames(geno_df)
  for (var in vars) {
    rate <- mean(!is.na(geno_df[[var]]))
    if (rate < opt$call_filter) {
      if (opt$verbose >= 1) {
        message(sprintf("Removing %s | call rate = %.4f", var, rate))
      }
      geno_df[[var]] <- NULL
    }
  }
}

# MAF filter
if (!is.na(opt$af_filter)) {
  vars <- colnames(geno_df)
  for (var in vars) {
    af <- sum(geno_df[[var]], na.rm = TRUE) / (2*(sum(!is.na(geno_df[[var]]))))
    if (af < opt$af_filter || af > (1-opt$af_filter)) {
      if (opt$verbose >= 1) {
        message(sprintf("Removing %s | AF = %.4f", var, af))
      }
      geno_df[[var]] <- NULL
    }
  }
}

# --------------------------------------------------------
# LD MATRIX
# --------------------------------------------------------

ld <- cor(geno_df)
ld_path <- paste0(opt$output_dir,"/HLA_LD.rds")
saveRDS(ld, ld_path)

# --------------------------------------------------------
# OUTPUT PED/MAP
# --------------------------------------------------------

# Split into allele and residue sets
allele_cols <- grep("^HLA_", colnames(geno_df), value = TRUE)
residue_cols <- setdiff(colnames(geno_df), allele_cols)
alleles <- geno_df[, c(allele_cols), drop = FALSE]
residues <- geno_df[, c(residue_cols), drop = FALSE]

# Cleanup IIDs
add_ids_from_rownames <- function(df) {
  rn <- rownames(df)
  spl <- strsplit(rn, ":", fixed = TRUE)
  df <- cbind(
    FID = vapply(spl, `[`, "", 2),
    IID = vapply(spl, `[`, "", 2),
    df
  )
  df
}
alleles   <- add_ids_from_rownames(alleles)
residues  <- add_ids_from_rownames(residues)

# Output
data_frames <- list(
  alleles = alleles,
  residues = residues
)
convert_dosage_to_geno <- function(x) {
  if (is.na(x)) {
    c("0", "0")
  } else if (x == 0) {
    c("A", "A")
  } else if (x == 1) {
    c("A", "T")
  } else if (x == 2) {
    c("T", "T")
  } else {
    c("0", "0") # fallback
  }
}
for (name in names(data_frames)) {
  df <- data_frames[[name]]
  geno_vars <- setdiff(names(df), c("FID", "IID"))
  ped_header <- df[, c("FID", "IID")]
  ped_header$PID <- 0
  ped_header$MID <- 0
  ped_header$Sex <- 0
  ped_header$Pheno <- -9
  geno_mat <- do.call(cbind, lapply(df[, geno_vars], function(col) {
    t(sapply(col, convert_dosage_to_geno))
  }))
  ped <- cbind(ped_header, geno_mat)
  write.table(ped, file = paste0(opt$output_dir, "/", name, ".ped"), quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
  map <- data.frame(
    CHR = 6,
    SNP = geno_vars,
    GENPOS = 0,
    BP = seq_along(geno_vars)
  )
  write.table(map, file = paste0(opt$output_dir, "/", name, ".map"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

# --------------------------------------------------------
# WRITE PLINK BINARY FILES
# --------------------------------------------------------

for (dat in names(data_frames)) {
  arg <- paste0(opt$plink, " --file ", opt$output_dir,"/",dat, " --keep-allele-order --make-bed --out ", opt$output_dir,"/",dat,"_tmp")
  system(arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)
}

# Force T/A encoding irrespective of MAF
for (dat in names(data_frames)) {
  arg <- paste0("awk 'BEGIN{OFS=\"\\t\"} {print $0, \"T\"}' ", opt$output_dir,"/",dat,"_tmp.bim > ", opt$output_dir,"/",dat,"_tmp.txt")
  system(arg)

  arg <- paste0(opt$plink, " --bfile ", opt$output_dir,"/",dat,"_tmp", " --a1-allele ", opt$output_dir,"/",dat,"_tmp.txt 7 2", " --make-bed --out ", opt$output_dir,"/",dat)
  system(arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)
  
  arg <- paste0("rm ",opt$output_dir,"/",dat,"_tmp*")
  system(arg)
}

# --------------------------------------------------------
# CLEANUP
# --------------------------------------------------------

if ( opt$verbose >= 1 ) cat("Cleaning up\n")
cleanup()