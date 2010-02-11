\name{bgAdjustBgp}
\alias{bgAdjustBgp}
\title{
Remove background
}
\description{
Estimate and remove background signal using anti-genomic background probes
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
# See normalizeBetweenSamples
}
