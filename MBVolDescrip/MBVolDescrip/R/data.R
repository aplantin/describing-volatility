#' Moving Pictures gut timeseries OTU table 
#'
#' Gut microbiome profiles for subjects M3 and F4 in Caporaso et al.'s Moving
#'  Pictures study (all time points). 
#'
#' @format ## `mp_otu`
#' A data frame with 3,962 rows (taxa) and 467 columns (samples). Metadata 
#' associated with sample IDs is included in mp_meta. 
#' 
#' @source MicrobeDS R package
"mp_otu"



#' Moving Pictures gut timeseries metadata
#'
#' Metadata associated with gut microbiome profiles for subjects M3 and F4 in
#'  Caporaso et al.'s Moving Pictures study (all time points). 
#'
#' @format ## `mp_meta`
#' A data frame with 467 rows (samples) and 7 columns: 
#' \describe{
#'   \item{sampID}{Sample ID, corresponding to column names of OTU table}
#'   \item{age}{Subject age (years)}
#'   \item{collectiontime}{Date of sample collection}
#'   \item{subjID}{Subject ID}
#'   \item{day}{Current day} 
#'   \item{firstday}{Earliest sampling date for that subject}
#'   \item{time}{Number of days elapsed since first sampling date for that subject, calculated as day - firstday}
#' }
#' 
#' @source MicrobeDS R package
"mp_meta"


#' Moving Pictures phylogenetic tree
#'
#' Phylogenetic tree associated with gut microbiome profiles for subjects M3 and F4 in
#'  Caporaso et al.'s Moving Pictures study (all time points).
#'
#' @format ## `mp_tree`
#' A phylogenetic tree. 
#' 
#' @source MicrobeDS R package
"mp_tree"


#' Moving Pictures taxonomic assignment for OTUs 
#'
#' Data frame showing taxonomic assignment for each OTU in Caporaso et al.'s 
#' Moving Pictures study (all time points). 
#' 
#' @format ## `mp_tax`
#' A phylogenetic tree. 
#' 
#' @source MicrobeDS R package
"mp_tax"



#' Student Microbiome Project gut timeseries OTU table 
#'
#' Gut microbiome profiles for 58 college-aged adults who did not use antibiotics 
#' during the study period, taken from Flores et al. (2014)'s Student Microbiome 
#' Project dataset. 
#'
#' @format ## `smp_otu`
#' A data frame with 10,509 rows (taxa) and 496 columns (samples). Metadata 
#' associated with sample IDs is included in smp_meta. 
#' 
#' @source <https://github.com/caporaso-lab/student-microbiome-project>
"smp_otu"



#' Student Microbiome Project gut timeseries metadata
#' 
#' Metadata associated with gut microbiome profiles for 58 college-aged adults who
#' did not use antibiotics during the study period, taken from Flores et al. 
#' (2014)'s Student Microbiome Project dataset. 
#'
#' @format ## `smp_meta`
#' A data frame with 496 rows (samples) and 17 columns, including sampID (sample ID, 
#' to match OTU table column names); subjID (subject ID); and time (in days), 
#' calculated as WeeksSinceStart * 7. 
#' 
#' @source <https://github.com/caporaso-lab/student-microbiome-project>
"smp_meta"


#' Student Microbiome Project gut timeseries metadata
#' 
#' Taxonomic assignment for OTUs in the gut microbiome profiles for 58 college-aged
#'  adults who did not use antibiotics during the study period, taken from Flores et al. 
#' (2014)'s Student Microbiome Project dataset. 
#'
#' @format ## `smp_tax_key`
#' A data frame with 1537 rows (genus-level taxa) and 2 columns: OTU ID (to match 
#' dimension names of OTU table) and taxonomic assignment corresponding to that OTU ID. 
#' 
#' @source <https://github.com/caporaso-lab/student-microbiome-project>
"smp_tax_key"


#' Gajer et al., Healthy Vaginal Microbiome OTU table
#'
#' Gajer et al. (2012) sampled the vaginal microbiome of 32 healthy women 
#' twice-weekly for 16 weeks, sequencing the V1-V2 region of the 16S rRNA gene. 
#' Participants reported menstrual bleeding, sexual activity, medications and 
#' contraceptives, and other characteristics in daily diaries. OTU counts were 
#' computed based on the original study’s reported OTU proportions and total 
#' read counts, rounded to the nearest whole number. 
#'
#' @format ## `gaj_otu`
#' A data frame with 937 rows (samples) and 331 columns (taxa). Metadata 
#' associated with sample IDs is included in gaj_meta. Data access and preparation 
#' is described on GitHub at https://github.com/aplantin/describing-volatility. 
#' 
#' @source <https://doi.org/10.1126/scitranslmed.3003605>
"gaj_otu"



#' Gajer et al., Healthy Vaginal Microbiome metadata
#' 
#' Metadata associated with Gajer et al. (2012) OTU data.
#'
#' @format ## `gaj_meta`
#' A data frame with 937 rows (samples) and 10 columns, including sampID (sample ID, 
#' to match OTU table column names); subjID (subject ID); and time (in days). 
#' Data access and preparation is described on GitHub at 
#' https://github.com/aplantin/describing-volatility. 
#' 
#' @source <https://doi.org/10.1126/scitranslmed.3003605>
"gaj_meta"





#' Ravel et al., Vaginal Microbiome (pre-BV) OTU data
#'
#' Ravel et al. (2013) characterized the vaginal microbiome daily for 10 weeks 
#' in 4 women without bacterial vaginosis (BV), 6 women with asymptomatic BV (ABV), 
#' and 15 women with symptomatic BV (SBV), sequencing the V1-V3 regions of the 
#' 16S rRNA gene. For women with episodes of ABV or SBV, we include only 
#' samples prior to the BV episode. Notably, women later diagnosed with ABV or 
#' SBV typically have Lactobacillus-depleted vaginal microbiomes at earlier time 
#' points, despite the lack of an active BV diagnosis based on Nugent score. 
#' We exclude women with fewer than 20 remaining time points.
#'
#' @format ## `rav_otu`
#' A data frame with 186 rows (samples) and 122 columns (taxa). Metadata 
#' associated with sample IDs is included in rav_meta. 
#' 
#' @source <https://microbiomejournal.biomedcentral.com/articles/10.1186/2049-2618-1-29 >
"rav_otu"



#' Ravel et al., Vaginal Microbiome (pre-BV) 
#' 
#' Metadata associated with Ravel et al. (2013) OTU data. Data access and preparation is 
#' described on GitHub at https://github.com/aplantin/describing-volatility. 
#'
#' @format ## `rav_meta`
#' A data frame with 186 rows (samples) and 21 columns, including sampID (sample ID, 
#' to match OTU table column names); subjID (subject ID); and time (in days). 
#' 
#' @source <https://microbiomejournal.biomedcentral.com/articles/10.1186/2049-2618-1-29 >
"rav_meta"