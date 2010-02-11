\name{readCharm}
\alias{readCharm}
\title{
Read in McrBC/CHARM DNA methylation microarray data
}
\description{
Read in DNA methylation microarray data from the McrBC/CHARM platform
}
\usage{
readCharm(files, path = ".", ut = "_532.xys", md = "_635.xys", sampleKey, sampleNames = NULL, pkgname, type = NULL, ...)
}
\arguments{
  \item{files}{
a vector of xys filenames
}
  \item{path}{
the path to the xys files
}
  \item{ut}{
the file ending that designates untreated channel files
}
  \item{md}{
the file ending that designates methyl-depleted channel files
}
  \item{sampleKey}{
a data frame with sample description information. One line per xys file.
}
  \item{sampleNames}{
a vector of names to use for the samples. One line per xys file.
}
  \item{pkgname}{
the annotation package name
}
  \item{type}{
deprecated option 
}
  \item{\dots}{
additional options passed on to read.xysfiles2
}
}
\details{
This function is a convenience wrapper to read.xysfiles2 to simplify reading in DNA methylation data from the Nimblegen McrBC/CHARM microarray platform. It makes guesses about the extensions used for the methyl-depleted (md) and untreated channels (ut).
}
\value{
A TilingFeatureSet object.
}
\references{
www.biostat.jhsph.edu/~maryee/charm
}
\author{
Martin Aryee <aryee@jhu.edu>
}

\seealso{
\code{\link{methp}}, \code{\link{dmrFinder}}
}
\examples{
# See normalizeBetweenSamples
}
