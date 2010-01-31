\name{dmrFdr}
\alias{dmrFdr}
\title{
Calculate FDR q-values for differentially methylated regions (DMRs)
}
\description{
Estimate false discovery rate q-values for a set of differentially methylated regions using a permutation approach.
}
\usage{
dmrFdr(dmr, compare = 1, numPerms = 1000, seed = NULL, verbose = TRUE)
}
\arguments{
  \item{dmr}{
a dmr object as returned by \code{\link{dmrFinder}}
}
  \item{compare}{
The dmr table for which to calculate DMRs. See details.
}
  \item{numPerms}{
Number of permutations
}
  \item{seed}{
Random seed (for reproducibility)
}
  \item{verbose}{
Boolean
}
}
\details{
This function estimates false discovery rate q-values for a dmr object returned by \code{\link{dmrFinder}}. dmrFinder can return a set of DMR tables with one or more pair-wise comparisons between groups. dmrFdr currently only calculated q-values for one of these at a time. The dmr table to use is specified by the compare option.
}
\value{
a list object in the same format as the input, but with extra p-val and q-val columns for the tabs element.
}
\author{
Martin Aryee <aryee@jhu.edu>
}


\seealso{
\code{\link{dmrFinder}}
}
\examples{
	if (require(charmData) & require(pd.feinberg.hg18.me.hx1)) {
		dataDir <- system.file("data", package="charmData")
		setwd(dataDir)
		# Read in raw data
		pd <- read.delim("sample_description_file.txt")
		rawData <- readCharm(files=pd$Filename, sampleKey=pd)
		# Estimate methylation
		p <- methp(rawData)
		# Find differentially methylated regions
		grp <- pData(rawData)$Tissue_Type
		dmr <- dmrFinder(rawData, p=p, groups=grp)
		names(dmr$tabs)
		head(dmr$tabs[[1]])
		# Estimate false discovery rate for DMRs
		dmr <- dmrFdr(dmr) 
		head(dmr$tabs[[1]])
	}
}
