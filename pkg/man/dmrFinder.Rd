\name{dmrFinder}
\alias{dmrFinder}
\title{
Find differentially methylated regions (DMRs)
}
\description{
Find differentially methylated regions (DMRs) from tiling microarray data.
}
\usage{
dmrFinder(eset = NULL, groups, p = NULL, l = NULL, chr = NULL, pos = NULL, pns = NULL, sdBins = NULL, controlIndex = NULL, Indexes = NULL, filter = NULL, package = NULL, ws = 7, verbose = TRUE, compare = "all", bgSubtract = TRUE, withinSampleNorm = "loess", binSize = 500, numSegments = 3, betweenSampleNorm = "quantile", minQCScore = NULL, controlProbes = c("CONTROL_PROBES", "CONTROL_REGIONS"), cluster = NULL, cutoff = 0.995, ...)
}
%- maybe also 'usage' for other objects documented here.
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
a matrix of methylation values (scale: -Inf, Inf). One column per sample
}
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
not currently used
}
  \item{controlIndex}{
vector of indices of non-CpG control probes
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
  \item{bgSubtract}{
boolean. Remove background signal?
}
  \item{withinSampleNorm}{
within-sample normalization method. "loess" or "none"
}
  \item{binSize}{
not currently used
}
  \item{numSegments}{
not currently used
}
  \item{betweenSampleNorm}{
between-sample normalization method. "quantile", "sqn" or "none"
}
  \item{minQCScore}{
not currently used
}
  \item{controlProbes}{
not currently used
}
  \item{cluster}{
not currently used
}
  \item{cutoff}{
not currently used
}
  \item{\dots}{
further options to be passed to methp
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a list of the same format as that from \code{\link{dmrFinder}} with the extra columns qval (Q-value) in the tabs slot.

}
\author{
Martin Aryee <aryee@jhu.edu>
}


\seealso{
\code{\link{dmrFinder}}
}
\examples{
## See 	\code{\link{dmrFinder}}
}
