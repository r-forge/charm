\name{cpgdensity}
\alias{cpgdensity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get CpG density for genomic regions
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cpgdensity(subject, chr, pos, windowSize = 500, sequence = "CG", showProgress = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{subject}{
%%     ~~Describe \code{subject} here~~
}
  \item{chr}{
%%     ~~Describe \code{chr} here~~
}
  \item{pos}{
%%     ~~Describe \code{pos} here~~
}
  \item{windowSize}{
%%     ~~Describe \code{windowSize} here~~
}
  \item{sequence}{
%%     ~~Describe \code{sequence} here~~
}
  \item{showProgress}{
%%     ~~Describe \code{showProgress} here~~
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (subject, chr, pos, windowSize = 500, sequence = "CG", 
    showProgress = FALSE) 
{
    idx <- split(1:length(chr), chr)
    s <- DNAString(sequence)
    cpgdensity <- rep(NA, length(pos))
    for (curchr in (names(idx))) {
        if (showProgress) 
            cat(".")
        if (curchr \%in\% names(subject)) {
            chrseq <- subject[[curchr]]
            curpos <- pos[idx[[curchr]]]
            v <- suppressWarnings(Views(chrseq, start = curpos - 
                windowSize/2, end = curpos + windowSize/2))
            d <- suppressWarnings(DNAStringSet(v))
            numcpg <- vcountPattern(s, d, fixed = TRUE)
            cpgdensity[idx[[curchr]]] <- numcpg/windowSize
        }
    }
    if (showProgress) 
        cat("\n")
    cpgdensity
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
