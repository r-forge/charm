\name{cpgdensity}
\alias{cpgdensity}
\title{
Get CpG density for genomic regions
}
\description{
Calculate the CpG density for a set of windows
}
\usage{
cpgdensity(subject, chr, pos, windowSize = 500, sequence = "CG")
}
\arguments{
  \item{subject}{
BSGenome object (e.g. Hsapiens)}
  \item{chr}{
character vector
}
  \item{pos}{
numeric vector
}
  \item{windowSize}{
number value
}
  \item{sequence}{
character string
}
}
\details{
Calculate the CpG density for a set of regions. chr and pos specify the region mid-points and windowSize specifies the size of the window to be centered on these mid-points. i.e. The window will stretch from pos-windowSize/2 to pos+windowSize/2.
}
\value{
a numeric vector 
}
\author{
Martin Aryee <aryee@jhu.edu>
}


\examples{
	if (require(BSgenome.Hsapiens.UCSC.hg18)){
		chr <- c("chr1", "chr1", "chr2")
		pos <- c(100000, 100500, 100000)
		cpgd <- cpgdensity(Hsapiens, chr=chr, pos=pos, windowSize = 500)
		cpgd
	}
}


