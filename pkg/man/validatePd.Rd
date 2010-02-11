\name{validatePd}
\alias{validatePd}
\title{
Validate a sample description file for two-color microarray data
}
\description{
Checks a sample description file describing two-color arrays for proper formatting and if requested guesses column numbers for file names, sample labels and group labels.
}
\usage{
validatePd(pd, fileNameColumn, sampleNameColumn, groupColumn, ut = "_532.xys", md = "_635.xys")
}
\arguments{
  \item{pd}{
A data frame containing the sample description table
}
  \item{fileNameColumn}{
Number or name of column containing file names (optional)
}
  \item{sampleNameColumn}{
Number or name of column containing sample names (optional)
}
  \item{groupColumn}{
Number or name of column containing group labels (optional)}
  \item{ut}{
the file ending that designates untreated channel files
}
  \item{md}{
the file ending that designates methyl-depleted channel files
}
}
\details{
This function checks the formatting of a sample description file to make sure it has suitable columns for file names, sample names and (optionally) group labels. The sample description file should have one line per channel, i.e. two lines per sample corresponding to the red and green channel data files. Values in the sample name column are used to pair the two channels together. If fileNameColumn, sampleNameColumn and/or groupColumn are unspecified a guess will be made.
}
\value{
If the input data frame is valid: a list containing the fileNameColumn, sampleNameColumn and groupColumn.
If the input data frame is invalid: FALSE
}
\author{
Martin Aryee <aryee@jhu.edu>
}

\seealso{
\code{\link{readCharm}}
}
\examples{
# See dmrFdr
}