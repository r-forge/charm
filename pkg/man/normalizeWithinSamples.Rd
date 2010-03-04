\name{normalizeWithinSamples}
\alias{normalizeWithinSamples}
\title{
Within-sample normalization for two-color data
}
\description{
Within-sample (between-channel) normalization for two-color DNA methylation microarray data.
}
\usage{
normalizeWithinSamples(dat, method = "loess", scale=c(0.99, 0.99), 
	controlProbes = c("CONTROL_PROBES", "CONTROL_REGIONS"), 
	controlIndex = NULL, approx=TRUE, breaks=1000, verbose=FALSE)
}

\arguments{
  \item{dat}{
a TilingFeatureSet
}
  \item{method}{
normalization method. "loess" or "none"
}
  \item{scale}{
a numeric vector (x,y). The xth percentile of each sample is scaled to represent y\% methylation. The default c(0.99, 0.99) means probes in the 99\% percentile represent 99\% methylation. Set to NA for no scaling.
}
  \item{controlProbes}{
character string of the label assigned to non-CpG control probes in the annotation file (i.e. the container column of the .ndf file).
}
  \item{controlIndex}{
a vector of non-CpG control probe indices
}
  \item{approx}{
Bin probes by signal intensity when loess normalizing. Much faster when TRUE
}
  \item{breaks}{
Number of bins to use when approx=TRUE
}
  \item{verbose}{
boolean: Verbose output?
}
}
\details{
This function is used by \code{\link{methp}} performs within-sample (between-channel) normalization. It is normally not used directly by the user.
}
\value{
a TilingFeatureSet
}

\author{
Martin Aryee <aryee@jhu.edu>, Rafael Irizarry
}

\examples{
# See normalizeBetweenSamples
}
