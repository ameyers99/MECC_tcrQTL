#####################################################
# Rscript: ImmunoXcan_compute-weights_multivariate.R
# Purpose: Compute HLA-based prediction models \
#          for multivariate TCR expression vector
# Author: Aaron Meyers
# Date: 2nd Oct 2025
#####################################################

# References
# - FUSION TWAS (PMID:26854917)
# - isoTWAS (PMID:38036788)

# --------------------------------------------------------
# Load packages and define options
# --------------------------------------------------------

set.seed(1999)

# Load packages
suppressMessages({
  suppressWarnings({
    library("optparse", quietly=TRUE, warn.conflicts=FALSE)
    library("plink2R", quietly=TRUE, warn.conflicts=FALSE)
    library("methods", quietly=TRUE, warn.conflicts=FALSE)
    library("tidyverse", quietly=TRUE, warn.conflicts=FALSE)
    library("data.table", quietly=TRUE, warn.conflicts=FALSE)
    library("isotwas", quietly=TRUE, warn.conflicts=FALSE)    
  })
})

# Command options
option_list = list(
  make_option("--input_file", action="store", default=NA, type='character', 
              help="Path to full TCR call file (space-separated, long format; columns: Sample, *pheno_name*, rate)"),
  make_option("--sample_file", action="store", default=NA, type='character', 
              help="Path to sample IDs for analysis (single-column, no header"),
  make_option("--output_dir", action="store", default=NA, type='character', 
              help="Path to output directory"), 
  make_option("--tmp", action="store", default=NA, type='character',
              help="Path to temporary files"),              
  make_option("--pheno_name", action="store", default=NA, type='character', 
              help="Name of molecular phenotype in input_file"), 
  make_option("--PATH_plink", action="store", default="plink", type='character',
              help="Path to plink executable [%default]"),              
  make_option("--bfile", action="store", default=NA, type='character',
              help="Path to PLINK binary input file prefix (minus bed/bim/fam)"),
  make_option("--covar", action="store", default=NA, type='character',
              help="Path to quantitative covariates (PLINK format) [optional]"),
  make_option("--resid", action="store_true", default=FALSE,
              help="Also regress the covariates out of the genotypes [default: %default]"),              
  make_option("--pe_filter", action="store", default=0.20, type='double', 
              help="Minimum expression prevalence for which to compute weights [default: %default]"), 
  make_option("--crossval", action="store", default=5, type='double',
              help="How many folds of cross-validation [default: %default]"),
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
  make_option("--models", action="store", default='mrce_lasso,multi_enet,joinet,spls,sgl,mtlasso,stacking', type='character',
              help="Comma-separated list of multivariate-response prediction models;\n
              will use all molecular phenotypes read from input file or specified in subset file [default: %default]\n
					Available models:\n
					mrce_lasso:\t Multivariate regression with covariance estimation)\n
					multi_enet:\t Multivariate elastic net\n
          joinet:\t Joinet stacked elastic net\n
					spls:\t Sparse partial least squares \n
					sgl:\t Sparse group lasso\n
					mtlasso:\t Multi-task lasso (L21 regularization)\n
					stacking:\t Super learner stacking\n")
)

 opt = parse_args(OptionParser(option_list=option_list))

 models = unique( c(unlist(strsplit(opt$models,',')),"multi_enet") ) # multi_enet always included
 M = length(models)
 if ( ! all(models %in% c('mrce_lasso','multi_enet','joinet','spls','sgl','mtlasso','stacking')) ) {
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

files = paste(opt$bfile,c(".bed",".bim",".fam"),sep='')
  if ( !is.na(opt$input_file) ) files = c(files,opt$input_file)
  if ( !is.na(opt$covar) ) files = c(files,opt$covar)
  for ( f in files ) {
        if ( !file.exists(f) ){
                cat( "ERROR: ", f , " input file does not exist\n" , sep='', file=stderr() )
                cleanup()
                q()
        }
}

# ------------------------------------------------------------------------
# MULTIVARIATE-RESPONSE PREDICTION MODEL FUNCTION
# ------------------------------------------------------------------------

weights.mv = function(genos, pheno, method) {

  sds = apply(genos, 2, sd)
  keep = sds != 0 & !is.na(sds)
  
  boot_list = lapply(1:10,function(i){jitter(pheno)})
  pheno_rep = rlist::list.rbind(boot_list)
  
  fit = compute_isotwas(X = genos[,keep], Y = pheno, Y.rep = pheno_rep, R = 10, 
                        id = rownames(pheno_rep), omega_est = 'replicates',
                        omega_nlambda = 10, method = c(method), predict_nlambda = 10,
                        family = "gaussian", scale = FALSE, alpha = 0.5, nfolds = opt$crossval,
                        verbose = FALSE, par = TRUE, n.cores = 5, seed = 1789,
                        run_all = FALSE, return_all = TRUE, tol.in = 0.001, maxit.in = 1000)
  
  if (!is.list(fit$best_models) ||
      is.null(fit$best_models$transcripts)) {
    warning("No valid isotwas output for this model/trait combination")
    return(list())
  }
  return(fit)

}


# ------------------------------------------------------------------------
# PREPARE PHENOTYPES
# ------------------------------------------------------------------------

# Make temporary directory for this run
  tmp_dir <- dirname(opt$tmp)
  arg <- paste("mkdir -p", tmp_dir)
  system(arg)

# Read in tcr file
tcr <- read.delim(opt$input_file,sep=" ")
tcr <- tcr %>% rename(IID = Sample)

# Transpose, with expression value 0 if missing
tcr <- tcr %>%
  tidyr::pivot_wider(names_from = opt$pheno_name, values_from = rate, values_fill = 0)

# Include all samples in the sample file, with expression value 0 if missing
all_samples <- read.delim(opt$sample_file, header=F)
all_samples <- all_samples %>% rename(IID = V1)
tcr <- all_samples %>%
  select(IID) %>%
  left_join(tcr, by = "IID")
tcr <- tcr %>%
  mutate(across(-IID, ~replace_na(., 0)))

# Filter phentoypes on prevalence
tcr <- (function(df) {
  cols <- setdiff(names(df), "IID")
  zp <- sapply(df[cols], function(x) mean(x == 0, na.rm = TRUE))
  df[, c("IID", cols[zp <= opt$pe_filter])]
})(tcr)
tcr.og <- tcr

# Regress out covariates 
if (opt$rn) {
  tcr[ , setdiff(names(tcr), "IID")] <- 
    lapply(tcr[ , setdiff(names(tcr), "IID")], rntransform)
}
if ( !is.na(opt$covar) ) {
	covar = ( read.table(opt$covar,as.is=T,head=T) )
	if ( opt$verbose >= 1 ) cat( "Loaded",ncol(covar)-2,"covariates\n")
  m = match(tcr$IID, covar$IID)
  m.keep = !is.na(m)
  tcr = tcr[m.keep, ]
  covar = covar[m, ]
  covar = covar[m.keep, ]
  for (j in 2:ncol(tcr)) {
  col <- colnames(tcr)[j]
  reg <- summary(lm(tcr[, j] ~ as.matrix(covar[, 3:ncol(covar)])))
  if (opt$verbose >= 1) {
    cat(reg$r.squared, "variance in phenotype", col, "explained by covariates\n")
  }
  tcr[, j] <- scale(reg$residuals)
}
}

# Reformat
rownames(tcr) <- tcr$IID
tcr$IID <- NULL
tcr <- as.matrix(tcr)


# ------------------------------------------------------------------------
# PREPARE GENOTYPES
# ------------------------------------------------------------------------

# Referent allele exclusion
geno.file = opt$tmp
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
arg = paste( opt$PATH_plink ," --allow-no-sex --bfile ",opt$bfile," --exclude ", ref_filt," --keep-allele-order --make-bed --out ",geno.file.alt,sep='') 
system(arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)



# Load genotypes
genos_alt = read_plink(geno.file.alt,impute="avg")
# - Standardise and scale
genos_alt$bed = scale(genos_alt$bed)
# - Check for NAs
nasnps = apply( is.na(genos_alt$bed) , 2 , sum )
if ( sum(nasnps) != 0 ) {
	cat( "WARNING :",sum(nasnps != 0),"alternate alleles could not be scaled and were zeroed out, make sure all are polymorphic\n" , file=stderr())
	genos_alt$bed[,nasnps != 0] = 0
}
# - Regress out covariates (more accurate but slower)
if ( !is.na(opt$covar) && opt$resid ) {
	if ( opt$verbose >= 1 ) cat("regressing covariates out of the genotypes\n")
	for ( i in 1:ncol(genos_alt$bed) ) {
		genos_alt$bed[,i] = summary(lm( genos_alt$bed[,i] ~ as.matrix(covar[,3:ncol(covar)]) ))$resid
	}
	genos_alt$bed = scale(genos_alt$bed)
}
N.tot = nrow(genos_alt$bed)
if ( opt$verbose >= 1 ) cat(nrow(genos_alt$bed),"genotyped samples, ",ncol(genos_alt$bed)," alternative alleles\n")

# Reformat IIDs
rownames(genos_alt$bed) <- sub(".*:", "", rownames(genos_alt$bed))

# ------------------------------------------------------------------------
# COMPUTE EXPRESSION PREDICTION MODELS
# ------------------------------------------------------------------------

if ( opt$verbose >= 1 ) cat("Computing full-sample weights\n")

# Align ID order
geno_mat <- as.matrix(genos_alt$bed)

common_ids <- intersect(rownames(geno_mat), rownames(tcr))
geno_mat <- geno_mat[common_ids, , drop = FALSE]
tcr <- tcr[common_ids, , drop = FALSE]

# Prediction models
for (mod in 1:M) {
  if ( models[mod] == "mrce_lasso" ) {
		wgt.matrix.mrce_lasso = weights.mv(geno_mat, tcr, 'mrce_lasso') 
	}
  if ( models[mod] == "multi_enet" ) {
		wgt.matrix.multi_enet = weights.mv(geno_mat, tcr, 'multi_enet') 
	}
  if ( models[mod] == "joinet" ) {
		wgt.matrix.joinet = weights.mv(geno_mat, tcr, 'joinet') 
	}
  if ( models[mod] == "spls" ) {
		wgt.matrix.spls = weights.mv(geno_mat, tcr, 'spls') 
	}
  if ( models[mod] == "sgl" ) {
		wgt.matrix.sgl = weights.mv(geno_mat, tcr, 'sgl') 
	}
  if ( models[mod] == "mtlasso" ) {
		wgt.matrix.mtlasso = weights.mv(geno_mat, tcr, 'mtlasso') 
	}
  if ( models[mod] == "stacking" ) {
		wgt.matrix.stacking = weights.mv(geno_mat, tcr, 'stacking') 
	}
}


# ------------------------------------------------------------------------
# OUTPUT
# ------------------------------------------------------------------------

# Store
for (trait in colnames(tcr)) {

  # Performance
  cv.performance = matrix(NA,nrow=2,ncol=M)
  rownames(cv.performance) = c("rsq","pval")
  colnames(cv.performance) = models
  for (mod in models) {
    df <- get(paste0("wgt.matrix.",mod))
  if (!is.null(df$best_models$transcripts[[trait]])) {
  cv.performance["rsq", mod] <- df$best_models$transcripts[[trait]]$r2
  cv.performance["pval", mod] <- df$best_models$transcripts[[trait]]$pvalue
} else {
  cv.performance["rsq", mod] <- NA
  cv.performance["pval", mod] <- NA
}
  }

  # Output weight matrix
  wgt.matrix = matrix(0, nrow = nrow(genos_alt$bim), ncol = M)
  colnames(wgt.matrix) = models
  rownames(wgt.matrix) = genos_alt$bim[,2]
  for (mod in models) {
  df <- get(paste0("wgt.matrix.", mod))
  trait_obj <- df$best_models$transcripts[[trait]]
  if (is.null(trait_obj) || is.null(trait_obj$weights)) {
    next
  }
  weights.df <- trait_obj$weights
  if (nrow(weights.df) > 0) {
    idx <- match(weights.df$SNP, rownames(wgt.matrix))
    wgt.matrix[idx, mod] <- weights.df$Weight
  }
}

  # Genotypes
  genotypes = genos_alt$bim

  # tcr.pe
  tcr.pe <- mean(tcr.og[, trait] != 0, na.rm = TRUE)

  # Output path
  out_name <- paste0(opt$weight_prefix,"_",trait)
  out_path <- file.path(opt$output_dir,out_name)

  # Save
  save(wgt.matrix, genotypes, cv.performance, N.tot, tcr.pe,
     file = paste0(out_path,".wgt.RDat"))

}


# CLEAN-UP
if ( opt$verbose >= 1 ) cat("Cleaning up\n")
cleanup()