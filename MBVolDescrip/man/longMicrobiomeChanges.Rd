% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LongFormatMbChanges.R
\name{longMicrobiomeChanges}
\alias{longMicrobiomeChanges}
\title{Makes long-format data frame of changes in the microbiome between identified pairs of samples}
\usage{
longMicrobiomeChanges(mbchanges)
}
\arguments{
\item{mbchanges}{List of microbiome change matrices, obtained from calcMicrobiomeChanges()}
}
\value{
Data frame containing the following columns:
\item{taxID}{Taxon identifier}
\item{changeID}{Identifier of pair of time points for a subject}
\item{subjID}{Subject ID}
\item{time1}{Earlier time point in change calculations}
\item{time2}{Later time point in change calculations}
\item{sampID1}{Sample ID for earlier time point}
\item{sampID2}{Sample ID for later time point}
\item{NominalLag}{Nominal time lag between samples}
\item{WindowWidth}{Allowable difference between actual and nominal time lag}
\item{AdditiveChg}{(t2-t1) relative abundance}
\item{MultNZChg}{(t2/t1) relative abundance, if both are greater than 0}
\item{NZLogFC}{log(t2/t1) relative abundance, if both are greater than 0}
\item{NZPairInd}{Indicator that both relative abundances are greater than 0, so fold change and log fold change are defined}
\item{QualChg}{Qualitative change: +1 indicates that the taxon appeared (t2 > 0, t1 = 0), and -1 indicates that the taxon disappeared (t2 = 0, t1 > 0). 0 indicates no qualitative change.}
\item{AvgTaxAbund}{Average relative abundance for the taxon (across all subjects & samples)}
\item{PropTaxNZ}{Proportion of observations for which taxon is nonzero (across all subjects & samples)}
}
\description{
This function makes a long-format data frame from the changes in relative abundance
(obtained from calcMicrobiomeChanges()).
}
