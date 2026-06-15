library("data.table")
library("parallel")
library("entropy")

cdr3 <- read.table("./data/tcr_phenotype/cdr3_seq.txt.gz",header=T)
setDT(cdr3)
cdr3[, L := nchar(tcr)]

# Sample-level NMI and Shannon entropy (template-weighted)
samps <- as.character(unique(cdr3$Sample))
l=15
all_res_list <- mclapply(samps, function(s) {
  df <- subset(cdr3, Sample==s & L==l)
  if(nrow(df) <= 1) return(NULL) # restrict to repertoires with >1 unique CDR3
  df[, rate := rate / sum(rate)]
  mat <- sapply(as.character(df$tcr), function(x){
    strsplit(x, "")[[1]]
  })
  mat <- t(mat)
  res <- data.frame()
  for(i in 1:l){
    for(k in 1:l){
      d1 <- data.frame(
        X = mat[, i],
        Y = mat[, k],
        rate = df$rate)
      d1 <- d1[d1$X != "NA" & d1$Y != "NA", ]
      xlist <- unique(as.character(d1$X))
      ylist <- unique(as.character(d1$Y))
      freq <- matrix(
        0,
        nrow = length(xlist),
        ncol = length(ylist),
        dimnames = list(xlist, ylist))
      for(x in xlist){
        for(y in ylist){
          both <- d1[d1$X == x & d1$Y == y, ]
          freq[x, y] <- sum(both$rate)
        }
      }
      freq <- freq / sum(freq)
      MI <- mi.empirical(freq, unit = "log2")
      Hx <- entropy.empirical(rowSums(freq), unit = "log2")
      Hy <- entropy.empirical(colSums(freq), unit = "log2")
      NMI <- 2 * MI / (Hx + Hy)
      dump <- data.frame(i=i+103, k=k+103, MI, Hx, Hy, NMI, Sample=s) # IMGT positions for L=15
      res <- rbind(res, dump)
      res$NMI[is.nan(res$NMI)] <- ifelse(
        res$i[is.nan(res$NMI)] == res$k[is.nan(res$NMI)],
        1,
        0)
    }
  }
  return(res)
}, mc.cores = 8)
all_res <- do.call(rbind, all_res_list)
saveRDS(all_res, "./data/tcr_phenotype/cdr3_covariance/cdr3_aa_covariance_weighted.rds")


# Sample-level NMI and Shannon entropy (unweighted)
all_res_list_2 <- mclapply(samps, function(s) {
  d1 <- subset(cdr3, Sample==s & L==l)
  if(nrow(d1) <= 1) return(NULL) # restrict to repertoires with >1 unique CDR3
    mat <- sapply(as.character(d1$tcr), function(x){
    x <- as.character(x);
    x <- strsplit(x, "");
    x <- unlist( x )
  })
  mat <- t(mat) #cdr3 sequence matirx: 15 col matrix
  res <- data.frame()
  for( i in 1:l ){
    for( k in 1:l ){
      d1 <- mat[, c(i,k)]
      d1 <- as.data.frame(d1)
      colnames(d1) <- c("X","Y")
      d1 <- d1[ d1$X!="NA" &  d1$Y!="NA",]
      xlist <- as.character( unique(d1$X) )
      ylist <- as.character( unique(d1$Y) )
      freq <- matrix(0, nrow=length(xlist), ncol=length(ylist))
      row.names(freq) <- xlist
      colnames(freq) <- ylist
      for( x in xlist ){
        for( y in ylist ){
          both <- subset(d1, X==x & Y==y)
          freq[ x, y ] <- nrow( both )
        }
      }
      MI <- mi.empirical(freq, unit=c("log2"))
      Hx <- entropy.empirical( rowSums(freq), unit = "log2" )
      Hy <- entropy.empirical( colSums(freq), unit = "log2" )
      NMI <- 2 * MI / ( Hx + Hy ) #normalized mutual entropy
      dump <- data.frame(i=i+103, k=k+103, MI, Hx, Hy, NMI, Sample=s) # IMGT positions for L=15
      res <- rbind(res, dump)
      res$NMI[is.nan(res$NMI)] <- ifelse(
        res$i[is.nan(res$NMI)] == res$k[is.nan(res$NMI)],
        1,
        0)
    }
  }
  return(res)
}, mc.cores = 8)
all_res_2 <- do.call(rbind, all_res_list_2)
saveRDS(all_res_2, "./data/tcr_phenotype/cdr3_covariance/cdr3_aa_covariance_unweighted.rds")
