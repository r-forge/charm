\name{normalizeBetweenSamples}
\alias{normalizeBetweenSamples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Between-sample normalization
}
\description{
Between-sample normalization for two-color DNA methylation microarray data.
}
\usage{
normalizeBetweenSamples	(dat, m="allQuantiles", untreated="none",
	enriched="none", controlProbes=c("CONTROL_PROBES", "CONTROL_REGIONS"), 
	controlIndex=NULL, verbose=FALSE)
}

\arguments{
  \item{dat}{
a TilingFeatureSet object
}
  \item{m}{
normalization method for log-ratios. "allQuantiles" for full quantile normalization, or "none"
}
  \item{untreated}{
normalization method for the untreated channel. "complete", "allQuantiles" or "none"
}
  \item{enriched}{
normalization method for the untreated channel. "sqn", "allQuantiles" or "none"

}
  \item{controlProbes}{
character string of the label assigned to non-CpG control probes in the annotation file (i.e. the container column of the .ndf file).
}
  \item{controlIndex}{
a vector of non-CpG control probe indices
}
  \item{verbose}{
boolean: Verbose output?
}

}
\details{
This function is used by \code{\link{methp}} performs between-sample normalization. It is normally not used directly by the user.
}
\value{
a TilingFeatureSet
}
\author{
Martin Aryee <aryee@jhu.edu>
}

\seealso{
\code{\link{methp}}
}

\examples{
	if (require(charmData) & require(BSgenome.Hsapiens.UCSC.hg18)) {
		dataDir <- system.file("data", package="charmData")
		setwd(dataDir)
		pd <- read.delim("phenodata.txt")
		pd <- subset(pd, sampleID=="441_liver")
		rawData <- readCharm(files=pd$filename, sampleKey=pd)
		# Correct spatial artifacts
		dat <- spatialAdjust(rawData)
		# Remove background signal
		dat <- bgAdjustBgp(dat)
		# Find non-CpG control probes
		ctrlIdx <- getControlIndex(rawData, subject=Hsapiens)
		# Within-sample normalization
		dat <- normalizeWithinSamples(dat, controlIndex=ctrlIdx)
		# Within-sample normalization
		dat <- normalizeBetweenSamples(dat)
	}
}
