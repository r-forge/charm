\name{methp}
\alias{methp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate DNA methylation 
}
\description{
Estimate DNA methylation from McrBC/CHARM microarray data in terms of log-ratios or percentages.
}
\usage{
methp(dat, spatial = TRUE, bgSubtract = TRUE, withinSampleNorm = "loess", 
	scale = c(0.99, 0.99), betweenSampleNorm = "quantile", 
	controlProbes = c("CONTROL_PROBES", "CONTROL_REGIONS"), 
	controlIndex = NULL, commonMethPercentParams = NULL, 
	verbose = TRUE, returnM = FALSE, 
	plotDensity = NULL, plotDensityGroups = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
a TilingFeatureSet object
}
  \item{spatial}{
boolean indicating whether to correct spatial artefacts
}
  \item{bgSubtract}{
boolean indicating whether to estimate and remove background signal before computing log-ratios
}
  \item{withinSampleNorm}{
within-sample normalization method. Choices are "loess" and "none"
}
  \item{scale}{
a numeric vector (x,y). The xth percentile of each sample is scaled to represent y\% methylation. The default c(0.99, 0.99) means probes in the 99\% percentile represent 99\% methylation. 
}

  \item{betweenSampleNorm}{
between-sample normalization method. Choices are "quantile", "sqn", and "none". See Details for more fine-grained control.
}
  \item{controlProbes}{
character string of the label assigned to non-CpG control probes in the annotation file (i.e. the container column of the .ndf file).
}
  \item{controlIndex}{
a vector of non-CpG control probe indices
}
  \item{commonMethPercentParams}{
boolean indicating whether a common set of parameters should be used for all samples when converting M-values to percentage methylation.
}
  \item{verbose}{
boolean: Verbose output?
}
  \item{returnM}{
boolean. Return M-values without converting to percentage methylation estimates
}
  \item{plotDensity}{
if specified this is the filename of the pdf diagnostic density plots.
}
  \item{plotDensityGroups}{
numeric vector of group labels used to color lines in the diagnostic density plots (see plotDensity option)
}
}
\details{
This function provides probe-level estimates of percentage DNA methylation from CHARM microarray data.
}
\value{
A matrix of probe-level percentage methylation estimates, one column per sample.
}

\author{
Martin Aryee <aryee@jhu.edu>
}

\seealso{
\code{\link{readCharm}}
}
\examples{
# See dmrFdr
}