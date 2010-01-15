\name{countGC}
\alias{countGC}
\title{
Count probe GC content 
}
\description{
Return the GC content for each probe 
}
\usage{
countGC(dat, type = "pm")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
a TilingFeatureSet object 
}
  \item{type}{
pm or bg probes
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}


\seealso{
\code{\link{readCharm}}
}
\examples{
	if (require(charmExample) & require(pd.feinberg.hg18.me.hx1)) {
		dataDir <- system.file("data", package="charmExample")
		setwd(dataDir)
		pd <- read.delim("sample_description_file.txt")
		rawData <- readCharm(files=pd$Filename, sampleKey=pd)
		ngc <- countGC(rawData)
	}
}
