########################################################
# Rscript: ImmunoXcan_logit_sumstats.R
# Purpose: Compute marginal HLA association summary \
#          statistics for binary TCR expression vector
# Author: Aaron Meyers
# Date: 2nd Oct 2025
########################################################

# References
# - FUSION TWAS (PMID:26854917)

set.seed(1999)

# Load packages
suppressMessages({
  suppressWarnings({
    library("optparse", quietly=TRUE, warn.conflicts=FALSE)
    library("plink2R", quietly=TRUE, warn.conflicts=FALSE)
    library("methods", quietly=TRUE, warn.conflicts=FALSE)
    library("tidyverse", quietly=TRUE, warn.conflicts=FALSE)
    library("data.table", quietly=TRUE, warn.conflicts=FALSE)
    library("logistf", quietly=TRUE, warn.conflicts=FALSE)
  })
})

# Command options
option_list = list(
  make_option("--input_file", action="store", default=NA, type='character', 
              help="Path to full TCR call file (space-separated, long format; columns: Sample, *pheno_name*, rate)"),
  make_option("--output_dir", action="store", default=NA, type='character', 
              help="Path to output directory"), 
  make_option("--pheno_name", action="store", default=NA, type='character', 
              help="Name of molecular phenotype in input_file"), 
  make_option("--bfile", action="store", default=NA, type='character',
              help="Path to PLINK binary input file prefix (minus bed/bim/fam)"),
  make_option("--tmp", action="store", default=NA, type='character',
              help="Path to temporary files"),
  make_option("--PATH_plink", action="store", default="plink", type='character',
              help="Path to plink executable [%default]"),
  make_option("--PATH_gcta", action="store", default="gcta_nr_robust", type='character',
              help="Path to GCTA executable [%default]"),
  make_option("--PATH_multiallelic_grm", action="store", default="GRM_multiallelic.R", type='character',
              help="Path to multiallelic GRM R script [%default]"),
  make_option("--covar", action="store", default=NA, type='character',
              help="Path to quantitative covariates (PLINK format) [optional]"),          
  make_option("--pe_filter", action="store", default=0.10, type='double', 
              help="Minimum expression prevalence for which to compute weights [default: %default]"), 
  make_option("--hsq_set", action="store", default=NA, type='double',
              help="Skip heritability estimation and set hsq estimate to this value [optional]"),              
  make_option("--verbose", action="store", default=1, type="integer",
              help="How much chatter to print: 0=nothing; 1=minimal; 2=all [default: %default]"),
  make_option("--noclean", action="store_true", default=FALSE,
              help="Do not delete any temporary files (for debugging) [default: %default]"), 
  make_option("--variant", action="store", default="allele", type='character', 
              help="Variant type, either 'allele' or 'amino_acid' [default: %default]"),
  make_option("--heterodimer", action="store", default=NA, type='character',
              help="Optional heterodimer mode: DP, DQ, or both (comma-separated: 'DP,DQ') [default: none]"),
  make_option("--heterodimer_maf", action="store", default=0.01, type='double',
              help="Minimum MAF for heterodimer features (only used if --heterodimer is set) [default: 0.01]"),
  make_option("--omnibus", action="store_true", default=FALSE,
              help="Optional omnibus likelihood test for locus (gene-level for alleles, position-level for amino acids) [default: %default]"),
  make_option("--prefix", action="store", default=NA, type='character',
              help="Prefix to add to output files indicating e.g. tissue type or source"),
  make_option("--subset_file", action="store", default=NA, type='character', 
              help="Partitioned file with list of TCRs to process in this job [optional]") 
)

 opt = parse_args(OptionParser(option_list=option_list))

 opt$heterodimer <- if (!is.na(opt$heterodimer)) {
  strsplit(opt$heterodimer, ",")[[1]]
} else {
  character(0)
}

 if ( opt$verbose == 2 ) {
  SYS_PRINT = F
 } else {
  SYS_PRINT = T
 }


# ------------------------------------------------------------------------
# --- CLEANUP FUNCTION
# ------------------------------------------------------------------------

cleanup = function() {
  if ( !opt$noclean ) {
    parent_dir <- dirname(opt$tmp)
    arg <- paste("rm -rf", parent_dir)
    system(arg)
  }
}

# --------------------------------------------------------
# META-DATA
# --------------------------------------------------------

 # TCR phenotypes 
if (!is.na(opt$subset_file)) {
    tcr_list <- scan(opt$subset_file, what=character()) # Read in TCR names from subset file or; 
} else {
    tcr_list <- unique(fread(opt$input_file, header = TRUE)[[opt$pheno_name]]) # Read in all unique TCRs from input file of TCR calls
}

# ------------------------------------------------------------------------
# --- PREPARE PHENOTYPES AND COMMENCE WEIGHT ESTIMATION
# ------------------------------------------------------------------------

# I/O checks
if ( system( paste(opt$PATH_plink,"--help") , ignore.stdout=T,ignore.stderr=T ) != 0 ) {
	cat( "ERROR: plink could not be executed, set with --PATH_plink\n" , sep='', file=stderr() )
	cleanup()
	q()
}
if ( !is.na(opt$hsq_set) && system( opt$PATH_gcta , ignore.stdout=T,ignore.stderr=T ) != 0 ){
        cat( "ERROR: gcta could not be executed, set with --PATH_gcta\n" , sep='', file=stderr() )
        cleanup()
        q()
}


# Load fam file
fam = read.table(paste(opt$bfile,".fam",sep=''),as.is=T)
all_samples <- data.table(Sample = unique(fam$V2))


# Loop of TCR phenotypes
 for (tcr_pheno in tcr_list) {


  # Make temporary directory for this run
  tmp_dir <- dirname(opt$tmp)
  arg <- paste("mkdir -p", tmp_dir)
  system(arg)
   

  # Read in rows for phenotype i
  df_sub <- fread(opt$input_file, header = TRUE)[ get(opt$pheno_name) == tcr_pheno, .(Sample, rate) ] 
  df_sub$rate <- 1 # enforces a binary phenotype

  # Merge with all samples, fill missing phenotype expression rates with 0
  df_sub <- merge(all_samples, df_sub, by = "Sample", all.x = TRUE) 
  df_sub[is.na(rate), rate := 0]


  # Check missingness (rate == 0) and filter on expression prevalence
  tcr.pe <- mean(df_sub$rate > 0, na.rm = TRUE)
  if (tcr.pe < opt$pe_filter | tcr.pe > 1-opt$pe_filter) { 
     if (opt$verbose >= 1) {
         message(sprintf("Skipping %s - %.1f%% of samples have non-zero rate",
                         tcr_pheno, tcr.pe * 100))
     }
     next
  }


  # Signpost
  if ( opt$verbose >= 1 ) print(paste("Commencing ImmunoXcan weight estimation for ",tcr_pheno)) 


  # Prepare .pheno format
  df_sub[, IID := Sample]
  df_sub[, FID := IID]
  df_sub <- df_sub[, .(FID, IID, rate)]


  # Create safe filename
  input_dir <- dirname(opt$input_file)
  tmp_file <- file.path(input_dir, paste0(gsub("-", "_", tcr_pheno), ".pheno"))


  # Write temp file
  fwrite(df_sub, file = tmp_file, sep = "\t", quote = FALSE, col.names = FALSE)
  

  # Output path
  out_name <- paste0(opt$prefix,"_",tcr_pheno)
  out_path <- file.path(opt$output_dir,out_name)


  # Perform i/o checks here:
  files = paste(opt$bfile,c(".bed",".bim",".fam"),sep='')
  if ( !is.na(tmp_file) ) files = c(files,tmp_file)
  if ( !is.na(opt$covar) ) files = c(files,opt$covar)

  for ( f in files ) {
	if ( !file.exists(f) ){
		cat( "ERROR: ", f , " input file does not exist\n" , sep='', file=stderr() )
		cleanup()
		q()
	}
}


# Make/fetch the phenotype file
if ( !is.na(tmp_file) ) {
	pheno.file = tmp_file
	pheno = read.table(pheno.file,as.is=T)
	# Match up data
	m = match( paste(fam[,1],fam[,2]) , paste(pheno[,1],pheno[,2]) )
	m.keep = !is.na(m)
	fam = fam[m.keep,]
	m = m[m.keep]
	pheno = pheno[m,]
} else {
	pheno.file = paste(opt$tmp,".pheno",sep='')
	pheno = fam[,c(1,2,6)]
	write.table(pheno,quote=F,row.names=F,col.names=F,file=pheno.file)
}

# Load in the covariates if needed 
if ( !is.na(opt$covar) ) {
	covar = ( read.table(opt$covar,as.is=T,head=T) )
	if ( opt$verbose >= 1 ) cat( "Loaded",ncol(covar)-2,"covariates\n")
	# Match up data
	m = match( paste(fam[,1],fam[,2]) , paste(covar[,1],covar[,2]) )
	m.keep = !is.na(m)
	fam = fam[m.keep,]
	pheno = pheno[m.keep,]
	m = m[m.keep]
	covar = covar[m,]
}

geno.file = opt$tmp
# recode to the intersection of samples and new phenotype
arg = paste( opt$PATH_plink ," --allow-no-sex --bfile ",opt$bfile," --pheno ",pheno.file," --keep ",pheno.file," --keep-allele-order --make-bed --out ",geno.file,sep='')
system(arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)


# ------------------------------------------------------------------------
# --- HERITABILITY ANALYSIS
# ------------------------------------------------------------------------

if ( is.na(opt$hsq_set) ) {
if ( opt$verbose >= 1 ) cat("### Estimating heritability\n")

# 1. Generate GRM
arg = paste( "Rscript --vanilla ", opt$PATH_multiallelic_grm," --bfile ",geno.file," --tmp ",opt$tmp," --PATH_gcta ",opt$PATH_gcta," --variant ",opt$variant," --verbose ",opt$verbose,sep='' )
system(arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)

# 2. Estimate heritability
if ( !is.na(opt$covar) ) {
arg = paste( opt$PATH_gcta ," --grm ",opt$tmp," --pheno ",pheno.file," --qcovar ",opt$covar," --prevalence ",tcr.pe," --out ",opt$tmp," --reml --reml-no-constrain --reml-lrt 1",sep='')
} else {
arg = paste( opt$PATH_gcta ," --grm ",opt$tmp," --pheno ",pheno.file," --prevalence ",tcr.pe," --out ",opt$tmp," --reml --reml-no-constrain --reml-lrt 1",sep='')
}
system(arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)

# 3. Evaluate LRT and V(G)/Vp
if ( !file.exists( paste(opt$tmp,".hsq",sep='') ) ) {
	cat(opt$tmp,"does not exist, likely GCTA could not converge\n",file=stderr())
        hsq <- c(NA, NA)   # hsq and SE
        hsq.pv <- 1
    } else {
    hsq.file = read.table(file=paste(opt$tmp,".hsq",sep=''),as.is=T,fill=T) 
    hsq = as.numeric(unlist(hsq.file[hsq.file[,1] == "V(G)/Vp_L",2:3])) 
    hsq.pv = as.numeric(unlist(hsq.file[hsq.file[,1] == "Pval",2])) 
}


if ( opt$verbose >= 1 ) cat("Heritability (se):",hsq,"LRT P-value:",hsq.pv,'\n')

} else {
if ( opt$verbose >= 1 ) cat("### Skipping heritability estimate\n")
hsq = opt$hsq_set
hsq.pv = NA
}

# ------------------------------------------------------------------------
# --- MARGINAL SUMMARY STATISTICS 
# ------------------------------------------------------------------------

# Read genotype data
genos = read_plink(geno.file,impute="avg")
N.tot = nrow(genos$bed)

# Call heterodimers 
if (opt$variant == "allele") {

  if (length(opt$heterodimer) > 0 && is.null(opt$heterodimer_maf)) {
    stop("--heterodimer_maf must be set when --heterodimer is used")
  }

  if ("DP" %in% opt$heterodimer) { # DPA1-DPB1
    new_vars <- character()
    DPB1 <- grep("^HLA_DPB1_", colnames(genos$bed), value = TRUE)
    DPA1 <- grep("^HLA_DPA1_", colnames(genos$bed), value = TRUE)
    new_cols <- list()
    for (a in DPA1) {
      for (b in DPB1) {
        new_name <- paste(a, b, sep = "__")
        new_val  <- genos$bed[, a] * genos$bed[, b]
        freq <- sum(new_val, na.rm = TRUE) /
          (sum(!is.na(new_val)) * 2)
        if (freq >= opt$heterodimer_maf) {
          new_cols[[new_name]] <- new_val
          new_vars <- c(new_vars, new_name)
        }
      }
    }
    if (length(new_cols) > 0) {
      genos$bed <- cbind(genos$bed, do.call(cbind, new_cols))
    }
  }
if ("DQ" %in% opt$heterodimer) { # DQA1-DQB1
    new_vars <- character()
    DQB1 <- grep("^HLA_DQB1_", colnames(genos$bed), value = TRUE)
    DQA1 <- grep("^HLA_DQA1_", colnames(genos$bed), value = TRUE)
    new_cols <- list()
    for (a in DQA1) {
      for (b in DQB1) {
        new_name <- paste(a, b, sep = "__")
        new_val  <- genos$bed[, a] * genos$bed[, b]
        freq <- sum(new_val, na.rm = TRUE) /
          (sum(!is.na(new_val)) * 2)
        if (freq >= opt$heterodimer_maf) {
          new_cols[[new_name]] <- new_val
          new_vars <- c(new_vars, new_name)
        }
      }
    }
    if (length(new_cols) > 0) {
      genos$bed <- cbind(genos$bed, do.call(cbind, new_cols))
    }
  }
}

sds = apply(genos$bed,2,sd)
# important : genotypes are standardised and scaled here:
genos$bed = scale(genos$bed)

# Marker names
hlas <- as.character(colnames(genos$bed))

# Phenotype
pheno = read.table(pheno.file,as.is=T,header=F)
colnames(pheno) <- c("V1", "V2", "outcome")

if ( opt$verbose >= 1 ) cat(nrow(pheno),"phenotyped samples, ",nrow(genos$bed),"genotyped samples, ",ncol(genos$bed)," markers\n")

# Convert genotype matrix to data frame
geno_df <- as.data.frame(genos$bed)
tmp <- do.call(rbind, strsplit(rownames(geno_df), ":", fixed = TRUE))
geno_df$V1 <- tmp[, 1]
geno_df$V2 <- tmp[, 2]
geno_df <- merge(pheno,geno_df,by=c("V1","V2"), sort=F)

# Add covariates once
covs <- NULL
if (!is.na(opt$covar)) {
  covar <- read.table(opt$covar, as.is=T, header=T)
  colnames(covar)[1:2] <- c("V1", "V2")
  geno_df <- merge(geno_df, covar, by = c("V1", "V2"), sort=F)
  covs <- names(covar)[-(1:2)]
}

# Output object
sum.stats <- data.frame(
  tcr  = character(0),
  geno = character(0),
  beta = numeric(0),
  se   = numeric(0),
  pval = numeric(0),
  N    = integer(0),
  stringsAsFactors = FALSE
)

# Compute marginal associations
if (opt$verbose >= 1) cat("### Computing marginal associations\n")

tcr_name <- tcr_pheno

for (hla in hlas) {
  # Build formula
  if (length(covs) > 0) {
    form <- as.formula(paste("outcome ~", paste(c(hla, covs), collapse = " + ")))
  } else {
    form <- as.formula(paste("outcome ~", hla))
    }
   lr2_fit <- tryCatch(
    logistf(
      formula = form,
      data = geno_df,
      na.action = na.omit
    ),
    error = function(e) {
      message("Model failed for ", hla, ": ", e$message)
      return(NULL)
    }
  )

  beta <- NA
  se   <- NA
  pval <- NA
  n_samples <- NA

  if (!is.null(lr2_fit)) {

    n_samples <- nobs(lr2_fit)

    if (hla %in% names(lr2_fit$coefficients)) {
      i <- which(names(lr2_fit$coefficients) == hla)
      beta <- lr2_fit$coefficients[i]
      se   <- sqrt(lr2_fit$var[i, i])
      pval <- lr2_fit$prob[i]
    }
  }
  
  # Store results
  sum.stats <- rbind(
    sum.stats,
    data.frame(
      tcr  = tcr_name,
      geno = hla,
      beta = beta,
      se   = se,
      pval = pval,
      N    = n_samples,
      stringsAsFactors = FALSE
    )
  )
}

# ------------------------------------------------------------------------
# --- COMPUTE MAF
# ------------------------------------------------------------------------

arg = paste( opt$PATH_plink ," --allow-no-sex --bfile ",geno.file," --freq --keep-allele-order "," --out ",geno.file,"_maf",sep='')
system(arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)
frq_file <- paste0(geno.file,"_maf.frq")
sample_af <- read.table(frq_file,header=T)

# ------------------------------------------------------------------------
# --- OPTIONAL: LOCUS-LEVEL OMNIBUS TEST
# ------------------------------------------------------------------------

# Identify referent allele
ref <- sample_af %>% # NB: assumes plink files are consistently bi-allelic (presence/absence) encoded, without multi-residue variants 
  mutate(
    group = if (opt$variant == "allele") { 
      str_extract(SNP, "^[^_]+_[^_]+") # change if different naming format
    } else if (opt$variant == "amino_acid") { 
      str_extract(SNP, "^[^_]+_[^_]+") # change if different naming format
    }
  ) %>%
  group_by(group) %>%
  filter(
    if (opt$variant == "amino_acid") {
      dplyr::n() > 1
    } else {
      TRUE
    }
  ) %>%
  slice_max(order_by = MAF, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-group)

if (opt$verbose >= 1) {
  cat("Referent alleles and frequency:\n")  
  for (i in 1:nrow(ref)) {
    cat(ref$SNP[i], ref$MAF[i], "\n")
  }
}


# Omnibus test
if (opt$omnibus) {
if ( opt$verbose >= 1 ) cat("### Running locus-level omnibus tests\n")

# Drop heterodimers if defined
if (opt$variant == "allele") { 
  if (length(opt$heterodimer) > 0) { 
    drop_cols <- grepl("DQA1", names(geno_df)) & grepl("DQB1", names(geno_df)) |
             grepl("DPA1", names(geno_df)) & grepl("DPB1", names(geno_df))
    geno_df <- geno_df[, !drop_cols, drop = FALSE]
  }
}
# Drop referent allele
geno_df <- geno_df[, !(colnames(geno_df) %in% ref$SNP), drop = FALSE]

# Alternate alleles per locus
sel <- colnames(geno_df) %in% hlas
prefix <- sub("^([^_]+_[^_]+)_.*$", "\\1", colnames(geno_df)[sel])
groups <- split(which(sel), prefix)

# Output object
omnibus <- data.frame(
  tcr  = character(0),
  locus = character(0),
  pval   = numeric(0),
  df     = integer(0),
  plr_chisq = numeric(0), # penalised likelihood ratio test statistic
  N    = integer(0),
  stringsAsFactors = FALSE
)

# Loop through loci
for (g in names(groups)) {
  idx <- groups[[g]]
  idx_names <- colnames(geno_df)[idx]
  keep_cols <- unique(c(idx_names, covs, "outcome"))
  keep_cols <- intersect(keep_cols, colnames(geno_df))
  X <- geno_df[, keep_cols, drop = FALSE]
  
  if (length(idx_names) == 0) next
  
  # Build formula
  preds <- setdiff(colnames(X), "outcome")
  if (length(covs) > 0) {
  form1 <- as.formula(
    paste("outcome ~", paste(c(preds), collapse = " + "))
  )
  form2 <- as.formula(
    paste("outcome ~", paste(covs, collapse = " + "))
  )
  } else {
  form1 <- as.formula(
    paste("outcome ~", paste(preds, collapse = " + "))
  )
  form2 <- as.formula("outcome ~ 1")
  }

   lr2_fit1 <- tryCatch(
    logistf(
      formula = form1,
      data = X,
      na.action = na.omit
    ),
    error = function(e) {
      message("Model failed for full model for ", g, ": ", e$message)
      return(NULL)
    }
  )
  lr2_fit2 <- tryCatch(
    logistf(
      formula = form2,
      data = X,
      na.action = na.omit
    ),
    error = function(e) {
      message("Model failed for null model for ", g, ": ", e$message)
      return(NULL)
    }
  )

  pval   = NA
  df     = NA
  dev = NA
  n_samples    = NA

  if (!is.null(lr2_fit1)) {
    n_samples <- nobs(lr2_fit1)
    lrt <- anova(lr2_fit2, lr2_fit1)
    pval <- lrt$pval
    df   <- lrt$df
    chisq <- lrt$chisq # delta deviance not make sense here
  }

  omnibus <- rbind(omnibus, data.frame(
    tcr = tcr_name,
    locus = g,
    df     = df, 
    plr_chisq = chisq, # penalised likelihood ratio test statistic
    pval   = pval,
    N    = n_samples,
    stringsAsFactors = FALSE
  ))
  }
  rownames(omnibus) <- NULL
}


# ------------------------------------------------------------------------
# --- OUTPUT WRITING
# ------------------------------------------------------------------------

# save output
genotypes = genos$bim 
if (opt$omnibus) {
save(sum.stats, omnibus, sample_af, genotypes, hsq, hsq.pv, N.tot, tcr.pe,
     file = paste0(out_path, ".wgt.RDat"))
} else {
save(sum.stats, sample_af, genotypes, hsq, hsq.pv, N.tot, tcr.pe,
     file = paste0(out_path, ".wgt.RDat"))
}

# cleanup
if ( opt$verbose >= 1 ) cat("Cleaning up\n")
cleanup()

# Remove temporary files 
file.remove(tmp_file)

}
