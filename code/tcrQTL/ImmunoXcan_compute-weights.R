#####################################################
# Rscript: ImmunoXcan_compute-weights.R
# Purpose: Compute HLA-based prediction models \
#          and marginal summary statistics for \
#          univariate TCR expression vector
# Author: Aaron Meyers
# Date: 2nd Oct 2025
#####################################################

# References
# - FUSION TWAS (PMID:26854917)

# --------------------------------------------------------
# Load packages and define options
# --------------------------------------------------------

set.seed(1999)

# Load packages
suppressMessages({
  suppressWarnings({
    library("optparse", quietly=TRUE, warn.conflicts=FALSE)
    library("plink2R", quietly=TRUE, warn.conflicts=FALSE)
    library("glmnet", quietly=TRUE, warn.conflicts=FALSE)
    library("methods", quietly=TRUE, warn.conflicts=FALSE)
    library("tidyverse", quietly=TRUE, warn.conflicts=FALSE)
    library("data.table", quietly=TRUE, warn.conflicts=FALSE)
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
  make_option("--PATH_gemma", action="store", default="gemma", type='character',
              help="Path to GEMMA executable [%default]"),
  make_option("--covar", action="store", default=NA, type='character',
              help="Path to quantitative covariates (PLINK format) [optional]"),
  make_option("--resid", action="store_true", default=FALSE,
              help="Also regress the covariates out of the genotypes [default: %default]"),              
  make_option("--pe_filter", action="store", default=0.20, type='double', 
              help="Minimum expression prevalence for which to compute weights [default: %default]"), 
  make_option("--hsq_set", action="store", default=NA, type='double',
              help="Skip heritability estimation and set hsq estimate to this value [optional]"),
  make_option("--crossval", action="store", default=5, type='double',
              help="How many folds of cross-validation, 0 to skip [default: %default]"),
  make_option("--verbose", action="store", default=1, type="integer",
              help="How much chatter to print: 0=nothing; 1=minimal; 2=all [default: %default]"),
  make_option("--noclean", action="store_true", default=FALSE,
              help="Do not delete any temporary files (for debugging) [default: %default]"),
  make_option("--rn", action="store_true", default=FALSE,
              help="Rank-normalise the phenotype after all QC: [default: %default]"),		  
  make_option("--variant", action="store", default="allele", type='character', 
              help="Variant type, either 'allele' or 'amino_acid' [default: %default]"),
  make_option("--weight_prefix", action="store", default=NA, type='character',
              help="Prefix to add to weight output files indicating e.g. tissue type or source"),
  make_option("--heterodimer", action="store", default=NA, type='character',
              help="Optional heterodimer mode: DP, DQ, or both (comma-separated: 'DP,DQ') [default: none]"),
  make_option("--heterodimer_maf", action="store", default=0.01, type='double',
              help="Minimum MAF for heterodimer features (only used if --heterodimer is set) [default: 0.01]"),
  make_option("--omnibus", action="store_true", default=FALSE,
              help="Optional omnibus likelihood test for locus (gene-level for alleles, position-level for amino acids) [default: %default]"),              
  make_option("--models", action="store", default="blup,lasso,top1,enet", type='character',
              help="Comma-separated list of prediction models [default: %default]\n
					Available models:\n
					top1:\tTop eQTL (standard marginal eQTL Z-scores always computed and stored)\n
					blup:\t Best Unbiased Linear Predictor (dual of ridge regression)\n
					bslmm:\t Bayesian Sparse Linear Model (spike/slab MCMC)\n
					lasso:\t LASSO regression (with mixing parameter of 1)\n
					enet:\t Elastic-net regression (with mixing parameter of 0.5)\n"),	
  make_option("--subset_file", action="store", default=NA, type='character', 
              help="Partitioned file with list of TCRs to process in this job [optional]") 
)

 opt = parse_args(OptionParser(option_list=option_list))

 opt$heterodimer <- if (!is.na(opt$heterodimer)) {
  strsplit(opt$heterodimer, ",")[[1]]
  } else {
  character(0)
 }

 models = unique( c(unlist(strsplit(opt$models,',')),"top1") ) # top1 always included
 M = length(models)

 if ( ! all(models %in% c('blup', 'lasso', 'top1', 'enet', 'bslmm')) ) {
	cat( "ERROR: --models flag included invalid models\n" , sep='' , file=stderr() )
	q()
 }

 if ( opt$verbose == 2 ) {
  SYS_PRINT = F
 } else {
  SYS_PRINT = T
 }


# --------------------------------------------------------
# STATS FUNCTIONS
# --------------------------------------------------------

# GenABEL standardisation 
"ztransform" <- 
function(formula,data,family=gaussian) {
	if (missing(data)) {
		if(is(formula,"formula")) 
			data <- environment(formula)
		else  
			data <- environment()
#		wasdata <- 0
	} else {
		if (is(data,"gwaa.data")) {
			data <- data@phdata
		} 
		else if (!is(data,"data.frame")) {
			stop("data argument should be of gwaa.data or data.frame class")
		}
#		attach(data,pos=2,warn.conflicts=FALSE)
#		wasdata <- 1
	}
	
	if (is.character(family)) 
           family <- get(family, mode = "function", envir = parent.frame())
	if (is.function(family)) 
           family <- family()
	if (is.null(family$family)) {
           print(family)
           stop("'family' not recognized")
	}
	
	if ( is(try(formula,silent=TRUE),"try-error") ) { 
		formula <- data[[as(match.call()[["formula"]],"character")]] 
	}
	
	if (is(formula,"formula")) {
#		mf <- model.frame(formula,data,na.action=na.omit,drop.unused.levels=TRUE)
		mf <- model.frame(formula,data,na.action=na.pass,drop.unused.levels=TRUE)
		mids <- complete.cases(mf)
		mf <- mf[mids,]
		y <- model.response(mf)
		desmat <- model.matrix(formula,mf)
		lmf <- glm.fit(desmat,y,family=family)
#		if (wasdata) 
#			mids <- rownames(data) %in% rownames(mf)
#		else 
		resid <- lmf$resid
#		print(formula)
	} else if (is(formula,"numeric") || is(formula,"integer") || is(formula,"double")) {
		y <- formula
		mids <- (!is.na(y))
		y <- y[mids]
		resid <- y
		if (length(unique(resid))==1) stop("trait is monomorphic")
		if (length(unique(resid))==2) stop("trait is binary")
	} else {
		stop("formula argument must be a formula or one of (numeric, integer, double)")
	}
	y <- (resid-mean(resid))/sd(resid)
#	if (wasdata==1) detach(data)
	tmeas <- as.logical(mids)
	out <- rep(NA,length(mids))
	out[tmeas] <- y
	out
}


# GenABEL rank normalisation 
"rntransform" <-
		function(formula,data,family=gaussian) {
	if ( is(try(formula,silent=TRUE),"try-error") ) { 
		if ( is(data,"gwaa.data") ) data1 <- phdata(data)
		else if ( is(data,"data.frame") ) data1 <- data
		else stop("'data' must have 'gwaa.data' or 'data.frame' class")
		formula <- data1[[as(match.call()[["formula"]],"character")]] 
	}
	var <- ztransform(formula,data,family)
	out <- rank(var) - 0.5
	out[is.na(var)] <- NA
	mP <- .5/max(out,na.rm=T)
	out <- out/(max(out,na.rm=T)+.5)
	out <- qnorm(out)
	out
}

# --------------------------------------------------------
# Read in meta-data
# --------------------------------------------------------

 # TCR phenotypes 
if (!is.na(opt$subset_file)) {
    tcr_list <- scan(opt$subset_file, what=character()) # Read in TCR names from subset file or; 
} else {
    tcr_list <- unique(fread(opt$input_file, header = TRUE)[[opt$pheno_name]]) # Read in all unique TCRs from input file of TCR calls
}


 # Sample IDs 
 all_samples <- fread(opt$input_file)[, .(Sample)] %>% unique()  


# ------------------------------------------------------------------------
# --- PREDICTION MODEL FUNCTION
# ------------------------------------------------------------------------

# BSLMM (NB: Make sure BLUP/BSLMM weights are being scaled properly based on MAF)
weights.bslmm = function(input, bv_type, snp, out = NA) {
  if (is.na(out)) out = paste(input, ".BSLMM", sep = "")

  out_dir  = dirname(out)
  out_base = basename(out)

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  arg = paste(
    opt$PATH_gemma,
    " -miss 1 -maf 0 -r2 1 -rpace 1000 -wpace 1000",
    " -bfile ", input,
    " -bslmm ", bv_type,
    " -o ", out_base,
    " -outdir ", out_dir,
    sep = ""
  )

  system(arg, ignore.stdout = SYS_PRINT, ignore.stderr = SYS_PRINT)

  param_file = file.path(out_dir, paste0(out_base, ".param.txt"))
  eff = read.table(param_file, head = TRUE, as.is = TRUE)

  eff.wgt = rep(NA, length(snp))
  m = match(snp, eff$rs)
  m.keep = !is.na(m)
  m = m[m.keep]
  eff.wgt[m.keep] = (eff$alpha + eff$beta * eff$gamma)[m]

  return(eff.wgt)
}


# LASSO
weights.lasso = function( genos , pheno , alpha=1 ) {
        eff.wgt = matrix( 0 , ncol=1 , nrow=ncol(genos) )
        # remove monomorphics
        sds = apply( genos  , 2 , sd )
        keep = sds != 0 & !is.na(sds)
        lasso = cv.glmnet( x=genos[,keep] , y=pheno , alpha=alpha , nfold=5 , intercept=T , standardize=F )
        eff.wgt[ keep ] = coef( lasso , s = "lambda.min")[2:(sum(keep)+1)]
        return( eff.wgt )
}


# Marginal Z-scores (used for top1)
weights.marginal = function( genos , pheno , beta=F ) {
	if ( beta ) eff.wgt = t( genos ) %*% (pheno) / ( nrow(pheno) - 1)
	else eff.wgt = t( genos ) %*% (pheno) / sqrt( nrow(pheno) - 1 )
	return( eff.wgt )
}


# Elastic Net
weights.enet = function( genos , pheno , alpha=0.5 ) {
	eff.wgt = matrix( 0 , ncol=1 , nrow=ncol(genos) )
	# remove monomorphics
	sds = apply( genos  , 2 , sd )
	keep = sds != 0 & !is.na(sds)
	enet = cv.glmnet( x=genos[,keep] , y=pheno , alpha=alpha , nfold=5 , intercept=T , standardize=F )
	eff.wgt[ keep ] = coef( enet , s = "lambda.min")[2:(sum(keep)+1)]
	return( eff.wgt )
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


# ------------------------------------------------------------------------
# --- I/O CHECKS
# ------------------------------------------------------------------------

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

if ( sum(models=="bslmm" | models=="blup") != 0 && system( paste(opt$PATH_gemma,"-h") , ignore.stdout=T,ignore.stderr=T ) != 0 ){
	cat( "ERROR: gemma could not be executed, set with --PATH_gemma or remove 'bslmm' and 'blup' from models\n" , sep='', file=stderr() )
	cleanup()
	q()
 }


# ------------------------------------------------------------------------
# --- Compute ImmunoXcan weights over TCR phenotype i to k
# ------------------------------------------------------------------------

# Load fam file
fam = read.table(paste(opt$bfile,".fam",sep=''),as.is=T)


# Loop of TCR phenotypes
 for (tcr_pheno in tcr_list) {


  # Make temporary directory for this run
  tmp_dir <- dirname(opt$tmp)
  arg <- paste("mkdir -p", tmp_dir)
  system(arg)
   

  # Read in rows for phenotype i
  df_sub <- fread(opt$input_file, header = TRUE)[ get(opt$pheno_name) == tcr_pheno, .(Sample, rate) ] 


  # Merge with all samples, fill missing phenotype expression rates with 0
  df_sub <- merge(all_samples, df_sub, by = "Sample", all.x = TRUE)
  df_sub[is.na(rate), rate := 0]


  # Check missingness (rate == 0) and filter on expression prevalence
  tcr.pe <- mean(df_sub$rate > 0, na.rm = TRUE)
  if (tcr.pe < opt$pe_filter) { 
     if (opt$verbose >= 1) {
         message(sprintf("Skipping %s - only %.1f%% of samples have non-zero rate",
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
  tmp_file_res <- file.path(input_dir, paste0(gsub("-", "_", tcr_pheno), ".pheno.resid"))


  # Write temp file
  fwrite(df_sub, file = tmp_file, sep = "\t", quote = FALSE, col.names = FALSE)
  

  # Output path
  out_name <- paste0(opt$weight_prefix,"_",tcr_pheno)
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

if ( opt$rn ) {
	pheno[,3] = rntransform( pheno[,3] )
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
	reg = summary(lm( pheno[,3] ~ as.matrix(covar[,3:ncol(covar)]) ))
	if ( opt$verbose >= 1 ) cat( reg$r.sq , "variance in phenotype explained by covariates\n" )
	pheno[,3] = scale(reg$resid)
	raw.pheno.file = pheno.file
	pheno.file = paste(pheno.file,".resid",sep='')
	write.table(pheno,quote=F,row.names=F,col.names=F,file=pheno.file)
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
arg = paste( opt$PATH_gcta ," --grm ",opt$tmp," --pheno ",raw.pheno.file," --qcovar ",opt$covar," --out ",opt$tmp," --reml --reml-no-constrain --reml-lrt 1",sep='')
} else {
arg = paste( opt$PATH_gcta ," --grm ",opt$tmp," --pheno ",pheno.file," --out ",opt$tmp," --reml --reml-no-constrain --reml-lrt 1",sep='')
}
system(arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)

# 3. Evaluate LRT and V(G)/Vp
if ( !file.exists( paste(opt$tmp,".hsq",sep='') ) ) {
	cat(opt$tmp,"does not exist, likely GCTA could not converge\n",file=stderr())
        hsq <- c(NA, NA)   # hsq and SE
        hsq.pv <- 1
    } else {
    hsq.file = read.table(file=paste(opt$tmp,".hsq",sep=''),as.is=T,fill=T)
    hsq = as.numeric(unlist(hsq.file[hsq.file[,1] == "V(G)/Vp",2:3]))
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

# read in genotypes
genos = read_plink(geno.file,impute="avg")

# OPTIONAL: Call heterodimers 
if (opt$variant == "allele") {
  if (!is.na(opt$heterodimer) && is.null(opt$heterodimer_maf)) {
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
pheno = genos$fam[,c(1,2,6)]
pheno[,3] = scale(pheno[,3])

# check if any genotypes are NA
nasnps = apply( is.na(genos$bed) , 2 , sum )
if ( sum(nasnps) != 0 ) {
	cat( "WARNING :",sum(nasnps != 0),"alleles could not be scaled and were zeroed out, make sure all are polymorphic\n" , file=stderr())
	genos$bed[,nasnps != 0] = 0
}

# regress covariates out of the genotypes as well (this is more accurate but slower)
if ( !is.na(opt$covar) && opt$resid ) {
	if ( opt$verbose >= 1 ) cat("regressing covariates out of the genotypes\n")
	for ( i in 1:ncol(genos$bed) ) {
		genos$bed[,i] = summary(lm( genos$bed[,i] ~ as.matrix(covar[,3:ncol(covar)]) ))$resid
	}
	genos$bed = scale(genos$bed)
}

N.tot = nrow(genos$bed)
if ( opt$verbose >= 1 ) cat(nrow(pheno),"phenotyped samples, ",nrow(genos$bed),"genotyped samples, ",ncol(genos$bed)," markers\n")

if (opt$verbose >= 1) cat("### Computing marginal associations\n")
tcr_name <- tcr_pheno

sum.stats <- data.frame(
  tcr  = character(0),
  geno = character(0),
  beta = numeric(0),
  se   = numeric(0),
  pval = numeric(0),
  N    = integer(0),
  stringsAsFactors = FALSE
)

for (i in 1:ncol(genos$bed)) {
  snp_name <- colnames(genos$bed)[i]

  y <- pheno[[3]]
  x <- as.numeric(genos$bed[, i])

  # Fit the model 
  lm_fit <- lm(y ~ x, na.action = na.omit)
  lm_res <- summary(lm_fit)

  # Number of samples used = rows included in the model
  n_samples <- nobs(lm_fit)  # directly counts the rows used

  # Only extract beta if the predictor exists (non-constant)
  if ("x" %in% rownames(lm_res$coefficients)) {
    beta <- lm_res$coefficients["x", 1]
    se   <- lm_res$coefficients["x", 2]
    pval <- lm_res$coefficients["x", 4]
  } else {
    beta <- NA
    se   <- NA
    pval <- NA
  }

  sum.stats <- rbind(sum.stats, data.frame(
    tcr  = tcr_name,
    geno = snp_name,
    beta = beta,
    se   = se,
    pval = pval,
    N    = n_samples,
    stringsAsFactors = FALSE
  ))
}

# ------------------------------------------------------------------------
# --- REFERENT ALLELE EXCLUSION
# ------------------------------------------------------------------------

# Remove referent alleles
arg = paste( opt$PATH_plink ," --allow-no-sex --bfile ",opt$bfile," --freq --keep-allele-order "," --out ",geno.file,"_maf",sep='')
system(arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)

frq_file <- paste0(geno.file,"_maf.frq")
frq <- read.table(frq_file,header=T)

ref <- frq %>% # NB: assumes plink files are consistently bi-allelic (presence/absence) encoded, without multi-residue variants 
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

ref_filt <- paste0(geno.file,"_ref.filt")
write.table(ref$SNP, ref_filt,
            quote = FALSE, row.names = FALSE, col.names = FALSE)

geno.file.alt <- paste0(geno.file,"_alt")
arg = paste( opt$PATH_plink ," --allow-no-sex --bfile ",opt$tmp," --exclude ", ref_filt," --keep-allele-order --make-bed --out ",geno.file.alt,sep='') 
system(arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)



# Reload genotypes with scaling, standardisation, residualisation
genos_alt = read_plink(geno.file.alt,impute="avg")
sds = apply(genos_alt$bed,2,sd)
genos_alt$bed = scale(genos_alt$bed)
nasnps = apply( is.na(genos_alt$bed) , 2 , sum )
if ( sum(nasnps) != 0 ) {
	cat( "WARNING :",sum(nasnps != 0),"alternate alleles could not be scaled and were zeroed out, make sure all are polymorphic\n" , file=stderr())
	genos_alt$bed[,nasnps != 0] = 0
}
if ( !is.na(opt$covar) && opt$resid ) {
	if ( opt$verbose >= 1 ) cat("regressing covariates out of the genotypes\n")
	for ( i in 1:ncol(genos_alt$bed) ) {
		genos_alt$bed[,i] = summary(lm( genos_alt$bed[,i] ~ as.matrix(covar[,3:ncol(covar)]) ))$resid
	}
	genos_alt$bed = scale(genos_alt$bed)
}
N.tot = nrow(genos_alt$bed)
if ( opt$verbose >= 1 ) cat(nrow(pheno),"phenotyped samples, ",nrow(genos_alt$bed),"genotyped samples, ",ncol(genos_alt$bed)," alternative alleles\n")


# ------------------------------------------------------------------------
# --- OPTIONAL: LOCUS-LEVEL OMNIBUS TEST
# ------------------------------------------------------------------------

# NB: genotypes scaled above
# NB: phenotype (optionally) residualised and rank-normalised above

if (opt$omnibus) {

if ( opt$verbose >= 1 ) cat("### Running locus-level omnibus tests\n")

pheno = genos_alt$fam[,c(1,2,6)]
pheno[,3] = scale(pheno[,3])
y <- pheno[[3]]

prefix <- sub("^([^_]+_[^_]+)_.*$", "\\1", colnames(genos_alt$bed))
groups <- split(seq_along(prefix), prefix)

omnibus <- data.frame(
  tcr  = character(0),
  locus = character(0),
  pval   = numeric(0),
  df     = integer(0),
  delta_deviance = numeric(0),
  N    = integer(0),
  rsq = numeric(0),
  stringsAsFactors = FALSE
)

for (g in names(groups)) {
  idx <- groups[[g]]
  X <- genos_alt$bed[, idx, drop = FALSE]
  
  if (ncol(X) == 0) next
  
  dat <- data.frame(y = y, X)
  
  full_form <- as.formula(
    paste("y ~", paste(colnames(X), collapse = " + ")))
  
  fit_full <- lm(full_form, data = dat, na.action = na.omit)
  fit_null <- lm(y ~ 1, data = dat, na.action = na.omit)
  
  lrt <- anova(fit_null, fit_full)
  pval <- lrt$`Pr(>F)`[2]
  df   <- lrt$Df[2]
  dev = deviance(fit_null) - deviance(fit_full) # F-statistic
  n_samples <- nobs(fit_full)

  reg <- summary(fit_full)
  r2 <- reg$r.sq
  
  omnibus <- rbind(omnibus, data.frame(
    tcr = tcr_name,
    locus = g,
    df     = df,
    rsq = r2,
    delta_deviance = dev,
    pval   = pval,
    N    = n_samples,
    stringsAsFactors = FALSE
  ))
  }
}

# ------------------------------------------------------------------------
# --- CROSS-VALIDATION
# ------------------------------------------------------------------------

set.seed(1)
cv.performance = matrix(NA,nrow=2,ncol=M)
rownames(cv.performance) = c("rsq","pval")
colnames(cv.performance) = models

if ( opt$crossval <= 1 ) {
if ( opt$verbose >= 1 ) cat("### Skipping cross-validation\n")
} else {
if ( opt$verbose >= 1 ) cat("### Performing",opt$crossval,"fold cross-validation\n")
cv.all = pheno
N = nrow(cv.all)
cv.sample = sample(N)
cv.all = cv.all[ cv.sample , ]
folds = cut(seq(1,N),breaks=opt$crossval,labels=FALSE)

cv.calls = matrix(NA,nrow=N,ncol=M)

for ( i in 1:opt$crossval ) {
	if ( opt$verbose >= 1 ) cat("- Crossval fold",i,"\n")
	indx = which(folds==i,arr.ind=TRUE)
	cv.train = cv.all[-indx,]
	# store intercept
	intercept = mean( cv.train[,3] )
	cv.train[,3] = scale(cv.train[,3])
	
	# hide current fold
	cv.file = paste(opt$tmp,".cv",sep='')
	write.table( cv.train , quote=F , row.names=F , col.names=F , file=paste(cv.file,".keep",sep=''))	

	arg = paste( opt$PATH_plink ," --allow-no-sex --bfile ",geno.file.alt," --keep ",cv.file,".keep --out ",cv.file," --keep-allele-order --make-bed",sep='') 
	system(arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)

	for ( mod in 1:M ) {
		if ( models[mod] == "blup" ) {
			pred.wgt = weights.bslmm( cv.file , bv_type=2 , snp=genos_alt$bim[,2] ) 
		}
		else if ( models[mod] == "bslmm" ) {
			pred.wgt = weights.bslmm( cv.file , bv_type=1 , snp=genos_alt$bim[,2] )
		}		
		else if ( models[mod] == "lasso" ) {
			 pred.wgt = weights.lasso( genos_alt$bed[ cv.sample[ -indx ],] , as.matrix(cv.train[,3]) , alpha=1 )
		}
		else if ( models[mod] == "enet" ) {
			pred.wgt = weights.enet( genos_alt$bed[ cv.sample[ -indx ],] , as.matrix(cv.train[,3]) , alpha=0.5 )
		}		
		else if ( models[mod] == "top1" ) { 
			pred.wgt = weights.marginal( genos_alt$bed[ cv.sample[ -indx ],] , as.matrix(cv.train[,3,drop=F]) , beta=T )
			pred.wgt[ - which.max( pred.wgt^2 ) ] = 0
		}

		# predict from weights into sample
		pred.wgt[ is.na(pred.wgt) ] = 0
		
		if (models[mod] == "top1") {
			cv.calls[indx, mod] = genos_alt$bed[ cv.sample[ indx ] , ] %*% pred.wgt
			} else {
				cv.calls[indx, mod] = genos_alt$bed[ cv.sample[ indx ] , ] %*% pred.wgt
				}
	}
}

# compute rsq + P-value for each model
for ( mod in 1:M ) {
	if ( !is.na(sd(cv.calls[,mod])) && sd(cv.calls[,mod]) != 0 ) {
		reg = summary(lm( cv.all[,3] ~ cv.calls[,mod] ))
		cv.performance[ 1, mod ] = reg$adj.r.sq
		cv.performance[ 2, mod ] = reg$coef[2,4]
	} else {
		cv.performance[ 1, mod ] = NA
		cv.performance[ 2, mod ] = NA
	}
}
if ( opt$verbose >= 1 ) write.table(cv.performance,quote=F,sep='\t')
}


# ------------------------------------------------------------------------
# --- FULL ESTIMATION AND OUTPUT WRITING
# ------------------------------------------------------------------------

if ( opt$verbose >= 1 ) cat("Computing full-sample weights\n")

# call models to get weights
wgt.matrix = matrix(0,nrow=nrow(genos_alt$bim),ncol=M)
colnames(wgt.matrix) = models
rownames(wgt.matrix) = genos_alt$bim[,2]
for ( mod in 1:M ) {
	if ( models[mod] == "blup" ) {
		wgt.matrix[,mod] = weights.bslmm( geno.file.alt , bv_type=2 , snp=genos_alt$bim[,2] , out=geno.file.alt ) 
	}
	else if ( models[mod] == "bslmm" ) {
		wgt.matrix[,mod] = weights.bslmm( geno.file.alt , bv_type=1 , snp=genos_alt$bim[,2] , out=geno.file.alt ) 
	}		
	else if ( models[mod] == "lasso" ) {
		 wgt.matrix[,mod] = weights.lasso( genos_alt$bed , as.matrix(pheno[,3]) , alpha=1 )
	}
	else if ( models[mod] == "enet" ) {
		wgt.matrix[,mod] = weights.enet( genos_alt$bed , as.matrix(pheno[,3]) , alpha=0.5 )
	}	
	else if ( models[mod] == "top1" ) { 
		wgt.matrix[,mod] = weights.marginal( genos_alt$bed , as.matrix(pheno[,3]) , beta=F ) 
	}
}

# save weights output
genotypes = genos$bim 
if (opt$omnibus) {
save(sum.stats, omnibus, wgt.matrix, genotypes, cv.performance, hsq, hsq.pv, N.tot, tcr.pe,
     file = paste0(out_path, ".wgt.RDat"))
} else {
save(sum.stats, wgt.matrix, genotypes, cv.performance, hsq, hsq.pv, N.tot, tcr.pe,
     file = paste0(out_path, ".wgt.RDat"))
}
# --- CLEAN-UP
if ( opt$verbose >= 1 ) cat("Cleaning up\n")
cleanup()

# Remove temporary files 
file.remove(tmp_file)
file.remove(tmp_file_res)

}
