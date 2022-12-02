#' Identifies time point subsamples
#' 
#' This function cross-checks metadata samples with OTU samples, then identifies 
#' pairs of samples separated by a fixed amount of time (with the ability to 
#' specify an acceptable window around that fixed time difference). It returns a 
#' new metadata file with the beginning and ending time point and sample ID for 
#' each of the non-overlapping fixed-length windows.
#' 
#'
#' @param otus Matrix of taxon abundances (OTUs or ASVs). By default, rows are 
#' samples (taxa are columns); this can be changed using the taxaAreRows argument. 
#' Sample IDs should be appropriate dimension names and should match metadata. 
#' @param metadata Data frame containing metadata, including columns containing 
#' the subject ID ("subjID"), the sample ID ("sampID"), and time point. 
#' @param desired_spacing Desired spacing between consecutive samples (in the 
#' same units as the time column in metadata). For example, if time is measured 
#' in days, then desired_spacing = 7 would give weekly sampling. 
#' @param window_width Acceptable amount of "wiggle room" to each side of the 
#' desired spacing. For example, if a 6-8 day range is acceptable for nominally 
#' weekly sampling, use window_width = 1. 
#' @param taxaAreRows Indicates whether taxa are rows (versus columns). Default FALSE. 
#' @return New metadata file with columns subject ID, time point 1, time point 2, sample ID 1, and sample ID 2 
#' @export
#' 
#' @importFrom dplyr filter arrange
#' 
#' 
temporalSubsampleMeta <- function(otus, metadata, desired_spacing, window_width, taxaAreRows = FALSE) {
  
  # check OTU input 
  if (taxaAreRows & is.null(colnames(otus))) {
    stop("Please ensure OTU matrix column names are sample IDs")
  } else if (!taxaAreRows & is.null(rownames(otus))) {
    stop("Please ensure OTU matrix row names are sample IDs")
  }
  
  # standardize data format so taxa are columns 
  if (taxaAreRows) {
    otus2 <- t(otus) 
  } else {
    otus2 <- otus 
  }
  
  # check metadata 
  if (!is.data.frame(metadata)) {
    stop("This function expects metadata to be a data frame object")
  }
  if (!all(c("subjID", "sampID", "time") %in% colnames(metadata))) {
    stop("This function expects metadata to contain the columns subjID, sampID, and time")
  }
  if (!inherits(metadata$sampID, "character"))  {
    stop("Please ensure sampID column has class character (not, e.g., factor) so that sampIDs may be used to reference proper OTU profiles")
  }
  if (!all(metadata$sampID %in% rownames(otus2)) | !all(rownames(otus2) %in% metadata$sampID)) {
    warning("Not all metadata sample IDs match OTU sample IDs; OTU samples without matching metadata will be excluded")
  }
  if (!inherits(metadata$subjID, "character")) {
    warning("Subject IDs must be character class; converting now")
    metadata$subjID <- as.character(metadata$subjID)
  }
  
  # ensure metadata and OTUs are arranged by subject ID and time 
  # (and identically to each other) 
  metadata <- metadata %>% 
    arrange(subjID, time) 
  otus2 <- otus2[metadata$sampID, ]
  
  # find appropriate time points for each subject 
  newmeta <- data.frame(subjID = metadata$subjID, 
                        time1 = NA, time2 = NA, 
                        sampID1 = NA, sampID2 = NA)
  
  allsubj <- unique(metadata$subjID) 
  allsubj.times <- lapply(allsubj, FUN = function(sid) metadata$time[metadata$subjID == sid])
  names(allsubj.times) <- allsubj
  
  for (i in 1:length(allsubj)) {
    ## which rows of new metadata to fill in 
    which.idx <- which(newmeta$subjID == allsubj[i])
    
    ## fill in first time point (earliest observed for the subject)
    this.times <- allsubj.times[[i]] 
    newmeta[which.idx[1], "time1"] <- this.times[1]
    current.row = 1 
    current.time = this.times[1]
    current.time.j = 1 
    j = 2 
    
    ## iterate through subsequent time points looking for appropriate gaps 
    while (j <= length(this.times)) {
      
      # desired next time point 
      next.time = current.time + desired_spacing 
      
      if (any(this.times == next.time)) {
        
        # if exactly desired spacing is available, then locate it 
        # and update this row's time2 and next row's time1 
        j = which(this.times == next.time) 
        newmeta[which.idx[current.row], "time2"] = this.times[j] 
        current.row = current.row + 1 
        newmeta[which.idx[current.row], "time1"] = this.times[j] 
        current.time = this.times[j] 
        current.time.j = j 
        j <- j+1 
        
      } else if (any(abs(this.times - next.time) <= window_width)) {
        
        # if exact desired spacing is not available, check window 
        # choose time point closest to desired spacing 
        # in case of tie, choose earlier time point (closer to current time) 
        this.diffs <- abs(this.times - next.time) 
        check.times <- which(this.diffs <= window_width) 
        keep.time <- which.min(this.diffs[check.times])
        j <- check.times[keep.time] 
        
        # update 
        newmeta[which.idx[current.row], "time2"] = this.times[j] 
        current.row = current.row + 1 
        newmeta[which.idx[current.row], "time1"] = this.times[j] 
        current.time = this.times[j] 
        current.time.j = j 
        j <- j+1 
        
      } else {
        
        # if no time point within window is available, choose different time1 
        # (j=j+1) and try again 
        j <- current.time.j + 1 
        newmeta[which.idx[current.row], "time1"] = this.times[j] 
        current.time = this.times[j] 
        current.time.j = j 
        j <- j+1 
      }
    }
  }
  
  # keep only pairs for which both time points were found 
  newmeta <- newmeta %>% 
    filter(!is.na(time2))
  
  # fill in sample IDs from metadata
  for (i in 1:nrow(newmeta)) {
    newmeta$sampID1[i] <- metadata$sampID[metadata$subjID == newmeta$subjID[i] & 
                                            metadata$time == newmeta$time1[i]]
    newmeta$sampID2[i] <- metadata$sampID[metadata$subjID == newmeta$subjID[i] & 
                                            metadata$time == newmeta$time2[i]]
  }
  
  return(newmeta)
  
}
