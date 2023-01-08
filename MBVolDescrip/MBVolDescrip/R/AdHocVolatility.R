#' Summarizes measures of quantitative, qualitative, or global change in microbiome between time points 
#' 
#' This function calculates several ad hoc measures of microbiome volatility between 
#' two consecutive time points based on changes in relative abundance (from calcMicrobiomeChanges()). 
#'
#' @param mbchanges List of microbiome change matrices, obtained from calcMicrobiomeChanges()
#' @return Data frame containing the following columns: 
#' \item{AvgAbsAdditive}{Average absolute additive change in relative abundance (per taxon per subject).}
#' \item{AvgAbsAdditiveNZ}{Average absolute additive change in relative abundance *between time points with nonzero abundances* (per taxon per subject)}
#' \item{AvgMultNZ}{Average multiplicative change in relative abundance *between time points with nonzero abundances* (per taxon per subject)}
#' \item{AvgNZLogFC}{Average log fold change in relative abundance *between time points with nonzero abundances* (per taxon per subject)}
#' \item{NumPairsNZ}{For each taxon for each subject, number of pairs of consecutive time points for which the taxon had nonzero abundance at both time points.}
#' \item{NumPairs}{Total number of pairs of time points for each subject (identical for all taxa)}  
#' \item{PropQualChange}{Proportion of time points at which a taxon appeared or disappeared (per taxon per subject)}
#' \item{AvgTaxAbund}{Average relative abundance for the taxon (across all subjects and samples included in mbchanges)} 
#' \item{PropTaxNZ}{Proportion of samples for which this taxon has nonzero abundance (across all subjects and samples included in mbchanges)} 
#' @export
#' 
#' @import tidyverse tibble 
#' @importFrom dplyr rename all_of
#' @importFrom tidyr pivot_longer as_tibble
#' @importFrom stats time 
#' 
summMicrobiomeVolatility <- function(mbchanges) {
  
  # Elements of mbchanges: ChangeMeta, AddChangeMat, MultChangeMat, QualChangeMat, TaxCharacteristics
  res <- expand.grid(subjID = unique(mbchanges$ChangeMeta$subjID), 
                     taxID = colnames(mbchanges$AddChangeMat))
  
  # Handy references 
  subjIDs <- unique(mbchanges$ChangeMeta$subjID)
  
  # Nonzero pair matrix 
  nzpair.mat <- !is.na(mbchanges$MultChangeMat)
  
  # Average additive change 
  avgadd.res <- sapply(subjIDs, FUN = function(s) {
    if (length(grep(s, rownames(mbchanges$AddChangeMat))) == 0) {
      NA
    } else if (length(grep(s, rownames(mbchanges$AddChangeMat))) == 1) {
      abs(mbchanges$AddChangeMat[grep(s, rownames(mbchanges$AddChangeMat)), ])
    } else {
      apply(mbchanges$AddChangeMat[grep(s, rownames(mbchanges$AddChangeMat)), ], 2, 
            FUN = function(x) mean(abs(x)))
    }
  })
  colnames(avgadd.res) <- subjIDs
  avgadd <- avgadd.res %>% 
    as_tibble(rownames = "taxID") %>% 
    pivot_longer(all_of(subjIDs), names_to = "subjID", values_to = "AvgAbsAdditive")
  
  # Average additive change among nonzero pairs 
  nzadd.mat <- matrix(NA, 
                      nrow=nrow(mbchanges$AddChangeMat), 
                      ncol=ncol(mbchanges$AddChangeMat))
  nzadd.mat[nzpair.mat] <- mbchanges$AddChangeMat[nzpair.mat]
  rownames(nzadd.mat) <- rownames(mbchanges$AddChangeMat)
  colnames(nzadd.mat) <- colnames(mbchanges$AddChangeMat)
  nzadd.res <- sapply(subjIDs, FUN = function(s) {
    if (length(grep(s, rownames(nzadd.mat))) == 0) {
      NA
    } else if (length(grep(s, rownames(nzadd.mat))) == 1) {
      abs(nzadd.mat[grep(s, rownames(nzadd.mat)), ])
    } else {
      apply(nzadd.mat[grep(s, rownames(nzadd.mat)), ], 2, 
            FUN = function(x) mean(abs(x), na.rm=T))
    }
  })
  nzadd.res[is.nan(nzadd.res)] <- NA
  nzadd <- nzadd.res %>% 
    as_tibble(rownames = "taxID") %>% 
    pivot_longer(all_of(subjIDs), names_to = "subjID", values_to = "AvgAbsAdditiveNZ")
  
  # Average multiplicative change among nonzero pairs 
  nzmult.mat <- matrix(NA, nrow=nrow(mbchanges$MultChangeMat), ncol=ncol(mbchanges$MultChangeMat))
  nzmult.mat[nzpair.mat] <- mbchanges$MultChangeMat[nzpair.mat]
  rownames(nzmult.mat) <- rownames(mbchanges$MultChangeMat)
  colnames(nzmult.mat) <- colnames(mbchanges$MultChangeMat)
  nzmult.res <- sapply(subjIDs, FUN = function(s) {
    if (length(grep(s, rownames(nzmult.mat))) == 0) {
      NA
    } else if (length(grep(s, rownames(nzmult.mat))) == 1) {
      nzmult.mat[grep(s, rownames(nzmult.mat)), ]
    } else {
      apply(nzmult.mat[grep(s, rownames(nzmult.mat)), ], 2, 
            FUN = function(x) mean(x, na.rm=T))
    }
  })
  nzmult.res[is.nan(nzmult.res)] <- NA
  nzmult <- nzmult.res %>% 
    as_tibble(rownames = "taxID") %>% 
    pivot_longer(all_of(subjIDs), names_to = "subjID", values_to = "AvgMultNZ")
  
  # Average log fold change among nonzero pairs 
  logfc.mat <- matrix(NA, nrow=nrow(mbchanges$MultChangeMat), ncol=ncol(mbchanges$MultChangeMat))
  logfc.mat[nzpair.mat] <- log(mbchanges$MultChangeMat[nzpair.mat])
  rownames(logfc.mat) <- rownames(mbchanges$MultChangeMat)
  colnames(logfc.mat) <- colnames(mbchanges$MultChangeMat)
  logfc.res <- sapply(subjIDs, FUN = function(s) {
    if (length(grep(s, rownames(logfc.mat))) == 0) {
      NA
    } else if (length(grep(s, rownames(logfc.mat))) == 1) {
      logfc.mat[grep(s, rownames(logfc.mat)), ]
    } else {
      apply(logfc.mat[grep(s, rownames(logfc.mat)), ], 2, 
            FUN = function(x) mean(x, na.rm=T))
    }
  })
  logfc.res[is.nan(logfc.res)] <- NA
  logfc <- logfc.res %>% 
    as_tibble(rownames = "taxID") %>% 
    pivot_longer(all_of(subjIDs), names_to = "subjID", values_to = "AvgNZLogFC")
  
  # Number of nonzero pairs 
  nzpair.res <- sapply(subjIDs, FUN = function(s) {
    if (length(grep(s, rownames(nzpair.mat))) == 0) {
      NA
    } else if (length(grep(s, rownames(nzpair.mat))) == 1) {
      as.numeric(nzpair.mat[grep(s, rownames(nzpair.mat)), ]) 
    } else {
      apply(nzpair.mat[grep(s, rownames(nzpair.mat)), ], 2, 
            FUN = function(x) sum(x))
    }
  })
  rownames(nzpair.res) <- colnames(nzpair.mat) 
  nzpair <- nzpair.res %>% 
    as_tibble(rownames = "taxID") %>% 
    pivot_longer(all_of(subjIDs), names_to = "subjID", values_to = "NumPairsNZ")
  
  # Number of pairs total 
  pair.res <- sapply(subjIDs, FUN = function(s) {
    if (length(grep(s, rownames(nzpair.mat))) == 0) {
      NA
    } else if (length(grep(s, rownames(nzpair.mat))) == 1) {
      rep(1, ncol(nzpair.mat))
    } else {
      apply(nzpair.mat[grep(s, rownames(nzpair.mat)), ], 2, 
            FUN = function(x) length(x))
    }
  })
  rownames(pair.res) <- colnames(nzpair.mat) 
  pairdat <- pair.res %>% 
    as_tibble(rownames = "taxID") %>% 
    pivot_longer(all_of(subjIDs), names_to = "subjID", values_to = "NumPairs")
  
  # Proportion of pairs with qualitative change 
  qual.res <- sapply(subjIDs, FUN = function(s) {
    if (length(grep(s, rownames(mbchanges$QualChangeMat))) == 0) {
      NA
    } else if (length(grep(s, rownames(mbchanges$QualChangeMat))) == 1) {
      as.numeric(mbchanges$QualChangeMat[grep(s, rownames(mbchanges$QualChangeMat)), ] != 0) 
    } else {
      apply(mbchanges$QualChangeMat[grep(s, rownames(mbchanges$QualChangeMat)), ], 2, 
            FUN = function(x) mean(x != 0))
    }
  })
  rownames(qual.res) <- colnames(mbchanges$QualChangeMat) 
  qualdat <- qual.res %>% 
    as_tibble(rownames = "taxID") %>% 
    pivot_longer(all_of(subjIDs), names_to = "subjID", values_to = "PropQualChange")
  
  # Taxon characteristics 
  taxchars <- mbchanges$TaxCharacteristics 
  taxchars <- taxchars %>% 
    rename(AvgTaxAbund = avgAbund, 
           PropTaxNZ = propNZ) 
  taxchars$taxID <- rownames(taxchars) 
  
  # Prepare final result 
  out <- res %>% 
    merge(., avgadd, by=c("subjID", "taxID")) %>% 
    merge(., nzadd, by=c("subjID", "taxID")) %>% 
    merge(., nzmult, by=c("subjID", "taxID")) %>% 
    merge(., logfc, by=c("subjID", "taxID")) %>% 
    merge(., nzpair, by=c("subjID", "taxID")) %>% 
    merge(., pairdat, by=c("subjID", "taxID")) %>% 
    merge(., qualdat, by=c("subjID", "taxID")) %>% 
    merge(., taxchars, by=c("taxID"))
  
  # Return 
  return(out)
  
}