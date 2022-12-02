#' Makes long-format data frame of changes in the microbiome between identified pairs of samples
#' 
#' This function makes a long-format data frame from the changes in relative abundance 
#' (obtained from calcMicrobiomeChanges()). 
#'
#' @param mbchanges List of microbiome change matrices, obtained from calcMicrobiomeChanges()
#' @return Data frame containing the following columns: 
#' \item{taxID}{Taxon identifier}
#' \item{changeID}{Identifier of pair of time points for a subject}
#' \item{subjID}{Subject ID}
#' \item{time1}{Earlier time point in change calculations}
#' \item{time2}{Later time point in change calculations}
#' \item{sampID1}{Sample ID for earlier time point}
#' \item{sampID2}{Sample ID for later time point}
#' \item{AdditiveChg}{(t2-t1) relative abundance}
#' \item{MultNZChg}{(t2/t1) relative abundance, if both are greater than 0}  
#' \item{NZLogFC}{log(t2/t1) relative abundance, if both are greater than 0}
#' \item{NZPairInd}{Indicator that both relative abundances are greater than 0, so fold change and log fold change are defined} 
#' \item{QualChg}{Qualitative change: +1 indicates that the taxon appeared (t2 > 0, t1 = 0), and -1 indicates that the taxon disappeared (t2 = 0, t1 > 0). 0 indicates no qualitative change.}
#' \item{AvgTaxAbund}{Average relative abundance for the taxon (across all subjects & samples)}
#' \item{PropTaxNZ}{Proportion of observations for which taxon is nonzero (across all subjects & samples)}
#' @export
#' 
#' @import tidyverse tibble 
#' @importFrom dplyr rename all_of select mutate
#' @importFrom tidyr pivot_longer as_tibble
#' @importFrom stats time 
#' 
longMicrobiomeChanges <- function(mbchanges) {
  
  # Elements of mbchanges: ChangeMeta, AddChangeMat, MultChangeMat, QualChangeMat, TaxCharacteristics
  
  # Useful reference 
  otu_ids <- colnames(mbchanges$AddChangeMat) 
  
  # Additive changes: wide to long 
  addchange_wide <- mbchanges$AddChangeMat %>% 
    as_tibble() 
  addchange_wide$changeID = rownames(mbchanges$AddChangeMat)
  addchange_long <- addchange_wide %>% 
    pivot_longer(cols = all_of(otu_ids), 
                 names_to = "taxID", values_to = "AdditiveChg")

  # Additive changes: wide to long 
  logfc_wide <- mbchanges$MultChangeMat %>% 
    as_tibble() 
  logfc_wide$changeID = rownames(mbchanges$MultChangeMat)
  logfc_long <- logfc_wide %>% 
    pivot_longer(cols = all_of(otu_ids), 
                 names_to = "taxID", values_to = "MultNZChg") %>% 
    mutate(NZLogFC = log(MultNZChg), 
           NZPairInd = as.numeric(!is.na(MultNZChg))) 
  
  # Qualitative changes: wide to long 
  qualchange_wide <- mbchanges$QualChangeMat %>% 
    as_tibble() 
  qualchange_wide$changeID = rownames(mbchanges$QualChangeMat)
  qualchange_long <- qualchange_wide %>% 
    pivot_longer(cols = all_of(otu_ids), 
                 names_to = "taxID", values_to = "QualChg")
 
  # Taxon characteristics 
  taxchars <- mbchanges$TaxCharacteristics 
  taxchars <- taxchars %>% 
    rename(AvgTaxAbund = avgAbund, 
           PropTaxNZ = propNZ) 
  taxchars$taxID <- rownames(taxchars) 
  
  # Prepare final result 
  out <- mbchanges$ChangeMeta %>% 
    select(changeID, subjID, time1, time2, sampID1, sampID2) %>% 
    merge(., addchange_long, by=c("changeID"), all=T) %>% 
    merge(., logfc_long, by=c("changeID", "taxID"), all=T) %>% 
    merge(., qualchange_long, by=c("changeID", "taxID"), all=T) %>% 
    merge(., taxchars, by=c("taxID"), all=T)
  
  # Return 
  return(out)
  
}
