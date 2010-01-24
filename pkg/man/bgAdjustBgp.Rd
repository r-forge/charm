\name{bgAdjustBgp}
\alias{bgAdjustBgp}
\title{
Remove background
}
\description{
Background removal
}
\usage{
bgAdjustBgp(dat)
}
\arguments{
  \item{dat}{
a TilingFeatureSet
}
}
\details{
Background signal removal using a modified version of the RMA convolution model. The background signal level is estimated within GC-strata using anti-genomic background probes.
}
\value{
a TilingFeatureSet

}
\author{
Martin Aryee <aryee@jhu.edu>
}



\examples{
	if (require(charmData) & require(pd.feinberg.hg18.me.hx1)) {
		dataDir <- system.file("data", package="charmData")
		setwd(dataDir)
		pd <- read.delim("sample_description_file.txt")
		rawData <- readCharm(files=pd$Filename, sampleKey=pd)
		dat <- bgAdjustBgp(rawData)
	}

}
