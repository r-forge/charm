\name{bgAdjustBgp}
\alias{bgAdjustBgp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Remove background
}
\description{
Background removal
}
\usage{
bgAdjustBgp(dat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
a TilingFeatureSet
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
	if (require(charmExample) & require(pd.feinberg.hg18.me.hx1)) {
		dataDir <- system.file("data", package="charmExample")
		setwd(dataDir)
		pd <- read.delim("sample_description_file.txt")
		rawData <- readCharm(files=pd$Filename, sampleKey=pd)
		dat <- bgAdjustBgp(rawData)
	}

}
