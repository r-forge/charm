\name{dmrFinder}
\alias{dmrFinder}
\title{
Find differentially methylated regions (DMRs)
}
\description{
Find differentially methylated regions (DMRs) from tiling microarray data.
}
\usage{
dmrFinder(eset=NULL, groups, p=NULL, l=NULL, chr=NULL, pos=NULL, pns=NULL,
	sdBins=NULL, controlIndex=NULL, 
	controlProbes=c("CONTROL_PROBES", "CONTROL_REGIONS"), Indexes=NULL, 
	filter=NULL, package=NULL, ws=7, verbose=TRUE, compare="all",  
	withinSampleNorm="loess", betweenSampleNorm="quantile", 
	cutoff=0.995, sortBy="ttarea",...)
}
\arguments{
  \item{eset}{
a TilingFeatureSet
}
  \item{groups}{
a vector of group labels for the samples in eset 
}
  \item{p}{
a matrix of percentage methylation values (scale: 0, 1). One column per sample
}
  \item{l}{
a matrix of methylation values (scale: -Inf, Inf), typically log-ratios.}
  \item{chr}{
vector of chromosome labels for the probes in eset, p or l
}
  \item{pos}{
vector of chromosomal coordinates for the probes in eset, p or l
}
  \item{pns}{
vector of region names for the probes in eset, p or l
}
  \item{sdBins}{
not currently implemented
}
  \item{controlIndex}{
vector of indices of non-CpG control probes
}
  \item{controlProbes}{
not currently used
}
  \item{Indexes}{
not currently used
}
  \item{filter}{
smoothing window weights. See details
}
  \item{package}{
annotation package name
}
  \item{ws}{
smoothing window size parameter. See details.
}
  \item{verbose}{
Verbose progress reporting
}
  \item{compare}{
the groups between which to find DMRs.
}
  \item{withinSampleNorm}{
within-sample normalization method. "loess" or "none"
}
  \item{betweenSampleNorm}{
between-sample normalization method. "quantile", "sqn" or "none"
}
  \item{cutoff}{
t-statistic cutoff used to identify probes as being in a DMR
}
  \item{sortBy}{
sort column for the DMR table. "area" or "ttarea"
}
  \item{\dots}{
further options to be passed to methp
}
}
\details{
This function finds differentially methylated regions (DMRs). The sortby parameter can be used to sort the DMR by area (# probes x length), or t-statistic area (# probes x t-statistic)
}
\value{
A list with
\item{tabs}{A list of DMR tables, one per comparison with columns:
	\describe{
		\item{start}{start of DMR (bp)}
		\item{end}{end of DMR (bp)}
		\item{p1}{average percentage methylation of all probes between start and end for group 1}
		\item{p2}{average percentage methylation of all probes between start and end for group 2}
		\item{regionName}{name of the tiling region in which the DMR is found (These names come from the NDF file)}
		\item{indexStart}{index of first probe in DMR}
		\item{indexEnd}{index of last probe in DMR}
		\item{area}{(indexEnd-indexStart) x (p1-p2), i.e. length x average difference}	
		\item{ttarea}{(indexEnd-indexStart) x (average probe level t-stat for between group difference)}		
	}
}
\item{p}{A matrix of percentage methylation estimates (NOTE: the probe order may differ from that of the input p matrix since probes are sorted into chromosomal order)}
\item{l}{This contains methylation log-ratios if they were passed to the function. Otherwise it contains logit-transformed percentage methylation estimates}
\item{chr }{a vector of chromosomes corresponding to the rows of p and l}
\item{pos }{a vector of positions}
\item{pns }{a vector of probe region names}
\item{controlIndex }{a vector of control probe indices}
\item{gm }{Group medians of the l matrix}
\item{groups }{a vector of group labels}
\item{args }{the DMR finder parameter vector}
\item{comps }{the vector of pairwise group comparisons}
\item{package }{the array annotation package name}
}
\author{
Martin Aryee <aryee@jhu.edu>, Peter Murakami, Rafael Irizarry
}
\seealso{
\code{\link{readCharm}}, \code{\link{methp}}, \code{\link{dmrFdr}}
}
\examples{
# See dmrFdr
}
