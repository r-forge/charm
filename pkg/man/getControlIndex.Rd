\name{getControlIndex}
\alias{getControlIndex}
\title{
Get indices of control probes from CpG-free regions
}
\description{
Get indices of control probes from CpG-free regions. 
}
\usage{
getControlIndex(dat, controlProbes = c("CONTROL_PROBES", "CONTROL_REGIONS"), noCpGWindow = 1000, subject, onlyGood = FALSE, matrix = TRUE)
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
# See dmrFdr
}
