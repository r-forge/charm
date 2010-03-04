\name{qcReport}
\alias{qcReport}
\title{
Microarray quality report
}
\description{
Calculate microarray quality scores and produce an optional pdf report}
\usage{
qcReport(dat, file = NULL, utRange = c(30, 100), enRange = c(8, 12), 
	numProbes = 5e+05, blockSize)
}
\arguments{
  \item{dat}{
a TilingFeatureSet
}
  \item{file}{
name of output pdf file
}
  \item{utRange}{
color-scale range for the untreated channel plots
}
  \item{enRange}{
color-scale range for the methyl-depleted channel plots
}
  \item{numProbes}{
maximum number of probes to use for plots. If smaller than the number of probes on the array numProbes are chosen at random, speeding up calculations for high-density arrays with several million probes.
}
  \item{blockSize}{
The array is divided into a series of blockSize x blockSize rectangular blocks and the average signal level calculated for each. If blockSize is unspecified a size is chosen that gives about 1250 probes per block. 
}
}
\details{
This function calculates microarray quality scores and produces an optional pdf report. Three quality metrics are calculated for each array:
\describe{ 
\item{Average signal strength.}{The average percentile rank of untreated channel signal probes among the background (anti-genomic) probes. Since the untreated channel contains total DNA a successful hybridization would have strong signal for all untreated channel genomic probes.}
\item{Untreated channel signal standard deviation.}{The array is divided into a series of rectangular blocks and the average signal level calculated for each. Since probes are arranged randomly on the array there should be no large differences between blocks. Arrays with spatial artifacts have a larger standard deviation between blocks.}
\item{Methyl-depleted channel signal standard deviation}{}
}
}
\value{
a matrix with a row for each sample. The 3 columns contain array signal strength score, untreated channel standard deviation and methyl-depleted channel standard deviation.}
\author{
Martin Aryee <aryee@jhu.edu>
}

\examples{
	if (require(charmData)) {
		dataDir <- system.file("data", package="charmData")
		setwd(dataDir)
		pd <- read.delim("phenodata.txt")
		rawData <- readCharm(files=pd$filename, sampleKey=pd)
		qcReport(rawData, file="qcReport.pdf")
	}
}