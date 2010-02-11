\name{spatialAdjust}
\alias{spatialAdjust}
\title{
Correct spatial artifacts
}
\description{
Remove spatial artifacts from microarray data stored in TilingFeatureSet objects
}
\usage{
spatialAdjust(dat, blockSize, theta = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
TilingFeatureSet
}
  \item{blockSize}{
The array is divided into a series of blockSize x blockSize rectangular blocks and the average signal level calculated for each. If blockSize is unspecified a size is chosen that gives about 1250 probes per block. 
}
  \item{theta}{
smoothing parameter
}
}
\details{
The array is divided into a set of blockSize x blockSize squares. A kernel smoother is then used to even out spatial artifacts.
}
\value{
a TilingFeatureSet
}
\author{
Martin Aryee <aryee@jhu.edu>
}

\examples{
# See normalizeBetweenSamples
}