\name{methPercent}
\alias{methPercent}
\title{
Estimate percentage DNA methylation from log-ratios
}
\description{
Estimate percentage DNA methylation from log-ratios
}
\usage{
methPercent(m, ngc, commonParams = TRUE)
}
\arguments{
  \item{m}{
a matrix of M-values (methylation log-ratios). One column per sample.
}
  \item{ngc}{
a vector with GC-content of probes. Same length as nrow(m)
}
  \item{commonParams}{
boolean indicating whether a common set of parameters should be used for all samples when converting M-values to percentage methylation.
}
}
\details{
This function estimates percentage DNA methylation from normalized methylation log-ratios (M-values).
}
\value{
a matrix of percentage methylation estimates. Same dimensions as m
}
\author{
Martin Aryee <aryee@jhu.edu>
}

\examples{
	if (require(charmData) & require(BSgenome.Hsapiens.UCSC.hg18)) {
		dataDir <- system.file("data", package="charmData")
		setwd(dataDir)
		# Read in raw data
		pd <- read.delim("phenodata.txt")
		pd <- subset(pd, sampleID=="441_liver")
		rawData <- readCharm(files=pd$filename, sampleKey=pd)
		# Find non-CpG control probes
		ctrlIdx <- getControlIndex(rawData, subject=Hsapiens)	
		# Get normalized methylation log-ratios
		m <- methp(rawData, controlIndex=ctrlIdx, returnM=TRUE)
		# Estimate percentage methylation
		ngc <- countGC(rawData)
		p <- methPercent(m, ngc)
	}
}