\name{getControlIndex}
\alias{getControlIndex}
\title{
Get indices of control probes from CpG-free regions
}
\description{
Get indices of control probes from CpG-free regions. 
}
\usage{
getControlIndex(dat, controlProbes = c("CONTROL_PROBES", "CONTROL_REGIONS"), noCpGWindow = NULL, subject, onlyGood = FALSE, matrix = TRUE)
}
\arguments{
  \item{dat}{
TilingFeatureSet
}
  \item{controlProbes}{
vector of names used to denote control probes in the 'container' column of the Nimblegen annotation (ndf) file. Optional
}
  \item{noCpGWindow}{
Size of the window centered on the probe that must be CpG-free
}
  \item{subject}{
A BSgenome object
}
  \item{onlyGood}{
deprecated option
}
  \item{matrix}{
deprecated option
}
}
\details{
The probes can either be identified as control probes in the microarray annotation package, or alternatively the function will search the genome (given an appropriate BSgenome object) for suitable probes. 
}
\value{
a vector
}
\author{
Martin Aryee <aryee@jhu.edu>
}


\examples{
	if (require(charmData) & require(pd.feinberg.hg18.me.hx1)) {
		dataDir <- system.file("data", package="charmData")
		setwd(dataDir)
		pd <- read.delim("sample_description_file.txt")
		pd <- subset(pd, Individual=="441" & Tissue_Type=="liver")	
		rawData <- readCharm(files=pd$Filename, sampleKey=pd)
		idx <- getControlIndex(rawData)
	}
}
