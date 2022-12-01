#' Calculates amount of quantitative, qualitative, or global change in microbiome between time points 
#' 
#' This function calculates several measures of change in the microbiome between 
#' two consecutive time points based on relative abundance data, including: 
#' - Additive: t2-t1 (difference in taxon proportions. This is defined even when one or both proportions are zero. 
#' - Multiplicative: t2/t1 (ratio of taxon proportions). This will be NA for observations for which t1 = 0 and/or t2 = 0. 
#' - Qualitative: I(t2 > 0) - I(t1 > 0). This will be +1 if the taxon appeared (t2 > 0, t1 = 0); -1 if the taxon disappeared (t2 = 0, t1 > 0); and 0 otherwise. 
#' It also returns a distances between consecutive time points if a distance matrix is provided by the user. 
#' Finally, it calculates a few useful characteristics of the taxa themselves (average abundance across all samples, proportion of non-zero values). 
#' 
#' Note that the function does NOT perform rarefaction prior to computing proportions. 
#' The user will want to consider the possible impact of variable sampling depth 
#' on results and account for it through rarefaction or another approach prior to 
#' using these functions. 
#'
#' @param otus Matrix of taxon abundances (OTUs or ASVs). By default, rows are 
#' samples (taxa are columns); this can be changed using the taxaAreRows argument. 
#' Sample IDs should be appropriate dimension names and should match metadata. If 
#' the matrix does not contain relative abundances, it will be converted to proportions 
#' prior to calculating any volatility measures. 
#' @param changemeta Data frame containing metadata, including columns for 
#' the subject ID ("subjID"), time point 1 for the comparison ("time1"), time point 2 ("time2"), 
#' sample ID corresponding to time point 1 ("sampID1"), and sample ID corresponding to 
#' time point 2 ("sampID2"). The easiest way to generate this metadata file is through 
#' the function temporalSubsampleMeta(). 
#' @param Ds Distance matrix or list of distance matrices, if distances between 
#' comparison time points are desired. One column will be added to changemeta 
#' for each distance matrix. Default NULL. 
#' @param taxaAreRows Indicates whether taxa are rows (versus columns). Default FALSE. 
#' @return List containing the following elements: 
#' \item{ChangeMeta}{New metadata for changes, with unique subject + time pair identifier and (if applicable) columns for distances }
#' \item{AddChangeMat}{Matrix with additive changes. Taxa are columns. Rows are labeled with new identifier.}
#' \item{MultChangeMat}{Matrix with multiplicative changes. Taxa are columns. Rows are labeled with new identifier.}
#' \item{QualChangeMat}{Matrix with qualitative changes. Taxa are columns. Rows are labeled with new identifier.}
#' \item{TaxCharacteristics}{Data frame with taxon IDs, average relative abundance, percentile of average relative abundance, and proportion of non-zero values across all samples. }
#' @export
#' 
#' @import tidyverse
#' 
calcMicrobiomeChanges <- function(otus, changemeta, Ds = NULL, taxaAreRows) {
  
  # check OTU input 
  if (taxaAreRows & is.null(colnames(otus))) {
    stop("Please ensure OTU matrix column names are sample IDs")
  } else if (!taxaAreRows & is.null(rownames(otus))) {
    stop("Please ensure OTU matrix row names are sample IDs")
  }
  
  # standardize OTU format so taxa are columns 
  # and convert to proportions 
  if (taxaAreRows) {
    otus2 <- t(otus) 
  } else {
    otus2 <- otus 
  }
  otuprops <- otus2/rowSums(otus2)
  
  # check metadata 
  if (!is.data.frame(changemeta)) {
    stop("This function expects metadata to be a data frame object")
  }
  if (!all(c("subjID", "sampID1", "sampID2", "time1", "time2") %in% colnames(changemeta))) {
    stop("This function expects metadata to contain the columns subjID, time1, time2, sampID1, and sampID2. Consider using the temporalSubsampleMeta() function to generate this object.")
  }
  if (!inherits(changemeta$sampID1, "character")) {
    stop("Please ensure sampID1 column has class character (not, e.g., factor) so that sampIDs may be used to reference proper OTU profiles")
  }
  if (!inherits(changemeta$sampID2, "character")) {
    stop("Please ensure sampID2 column has class character (not, e.g., factor) so that sampIDs may be used to reference proper OTU profiles")
  }
  if (!all(changemeta$sampID1 %in% rownames(otuprops)) | !all(changemeta$sampID2 %in% rownames(otuprops))) {
    warning("Not all metadata sample IDs match OTU sample IDs. Please check that you are using the same OTU and changemeta files as from the temporalSubsampleMeta() function.")
  }
  
  # Generate time 1 and time 2 OTU proportion matrices 
  changemeta$changeID = paste0(changemeta$subjID, "_", changemeta$time1, "_", changemeta$time2)
  time1mat <- otuprops[changemeta$sampID1, ]
  time2mat <- otuprops[changemeta$sampID2, ]
  rownames(time1mat) = rownames(time2mat) = changemeta$changeID
  
  # remove taxa absent at both time points 
  abstime1 <- which(apply(time1mat, 2, FUN = function(x) sum(x != 0)) == 0)
  abstime2 <- which(apply(time2mat, 2, FUN = function(x) sum(x != 0)) == 0)
  abstax <- intersect(abstime1, abstime2) 
  time1mat <- time1mat[, -abstax]
  time2mat <- time2mat[, -abstax]
  
  # Matrix-based measures of change: ADDITIVE 
  AddChangeMat <- time2mat - time1mat 
  
  # Matrix-based measures of change: MULTIPLICATIVE
  MultChangeMat <- matrix(NA, nrow = nrow(AddChangeMat), ncol = ncol(AddChangeMat))
  MultChangeMat[time2mat > 0 & time1mat > 0] <- time2mat[time2mat > 0 & time1mat > 0] / time1mat[time2mat > 0 & time1mat > 0] 
  rownames(MultChangeMat) <- rownames(AddChangeMat)
  colnames(MultChangeMat) <- colnames(AddChangeMat)
  
  # Matrix-based measures of change: QUALITATIVE ONLY 
  QualChangeMat <- (time2mat > 0) - (time1mat > 0)
  
  # Distance between consecutive time points 
  if (!is.null(Ds)) {
    
    # Check distance input 
    if (is.matrix(Ds)) {
      Ds <- list(Ds)
      names(Ds) = c("DistanceMetric")
    } else if (!is.list(Ds)) {
      stop("Expect Ds to be either a single distance matrix or a list of such matrices.")
    } 
    if (!is.matrix(Ds[[1]])) {
      stop("Expect each element of Ds to be a matrix. Element 1 is not a matrix.")
    } 
    if (!all(changemeta$sampID1 %in% rownames(Ds[[1]]))) {
      stop("Ensure row and column names of distance matrices are sample IDs.")
    }
    
    # Add one column with distances to each row of changemeta
    if (!is.null(names(Ds))) {
      dist.names <- names(Ds)
    } else {
      dist.names <- paste0("Distance", 1:length(Ds))
    }
    changemeta[, dist.names] <- NA 
    
    # Fill in columns 
    for (i in 1:nrow(changemeta)) {
      sid1 <- changemeta$sampID1[i] 
      sid2 <- changemeta$sampID2[i]
      changemeta[i, dist.names] <- lapply(Ds, FUN = function(d) d[sid1, sid2])
    }
    
  }
  
  # Taxon characteristics 
  bothtimes_otus <- otuprops[union(changemeta$sampID1,changemeta$sampID2), -abstax]
  TaxCharacteristics <- data.frame(
    avgAbund = apply(bothtimes_otus, 2, FUN = function(x) mean(x)), 
    propNZ = apply(bothtimes_otus, 2, FUN = function(x) mean(x != 0))
  )
  
  # Return 
  return(list(ChangeMeta = changemeta, AddChangeMat = AddChangeMat, 
              MultChangeMat = MultChangeMat, QualChangeMat = QualChangeMat, 
              TaxCharacteristics = TaxCharacteristics))
  
}

