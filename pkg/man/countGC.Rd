\name{countGC}
\alias{countGC}
\title{
Count probe GC content 
}
\description{
Return the GC content for each probe 
}
\usage{
countGC(dat, type = "pm")
}
\arguments{
  \item{dat}{
a TilingFeatureSet object 
}
  \item{type}{
pm or bg probes
}
}
\details{
This function returns the sum of #G + #C in the pm or bg probes.
}
\value{
a numeric vector
}
\author{
Martin Aryee <aryee@jhu.edu>
}

\seealso{
\code{\link{readCharm}}
}
\examples{
	if (require(charmData)) {
		dataDir <- system.file("data", package="charmData")
		setwd(dataDir)
		pd <- read.delim("phenodata.txt")
		pd <- subset(pd, sampleID=="441_liver")
		rawData <- readCharm(files=pd$filename, sampleKey=pd)
		ngc <- countGC(rawData)
		head(ngc)
	}
}
