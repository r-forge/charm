\name{plotDensity}
\alias{plotDensity}
\title{
Log-ratio density plot for all probes and control probes
}
\description{
Make density plots of log-ratios for two-color microarray data. Two plots are produced: one for all probes on the array, and a second for the control probes.
}
\usage{
plotDensity(dat, rx = c(-4, 6), controlIndex = NULL, pdfFile = NULL, 
	main = NULL, cols = NULL, lwd = NULL)
}
\arguments{
  \item{dat}{
a TilingFeatureSet
}
  \item{rx}{
x-axis range
}
  \item{controlIndex}{
a vector of non-CpG control probe indices
}
  \item{pdfFile}{
name of output pdf file
}
  \item{main}{
main title 
}
  \item{cols}{
vector of colors, one per sample
}
  \item{lwd}{
vector of line widths, one per sample
}
}
\details{
This function makes density plots for a) all probes and b) control probes. It is typically called from within methp when a file name is specified for its plotDensity option. The plots are useful for identifying problematic outlier samples.
}
\value{
No return value. Called for its side-effect of producing a pdf plot.
}
\author{
Martin Aryee <aryee@jhu.edu>
}

\examples{
	if (require(charmData) & require(BSgenome.Hsapiens.UCSC.hg18)) {
		dataDir <- system.file("data", package="charmData")
		# Read in raw data
		pd <- read.delim(file.path(dataDir, "phenodata.txt"))
		rawData <- readCharm(path=dataDir, files=pd$filename, sampleKey=pd)
		ctrlIdx <- getControlIndex(rawData, subject=Hsapiens)
		plotDensity(rawData, controlIndex=ctrlIdx, pdfFile="density.pdf")	
	}
}