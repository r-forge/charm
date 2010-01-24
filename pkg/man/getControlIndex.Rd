\name{getControlIndex}
\alias{getControlIndex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get indices of control probes from CpG-free regions
}
\description{
Get indices of control probes from CpG-free regions. 
}
\usage{
getControlIndex(dat, controlProbes = c("CONTROL_PROBES", "CONTROL_REGIONS"), noCpGWindow = NULL, subject, onlyGood = FALSE, matrix = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
%%     ~~Describe \code{dat} here~~
}
  \item{controlProbes}{
%%     ~~Describe \code{controlProbes} here~~
}
  \item{noCpGWindow}{
%%     ~~Describe \code{noCpGWindow} here~~
}
  \item{subject}{
%%     ~~Describe \code{subject} here~~
}
  \item{onlyGood}{
%%     ~~Describe \code{onlyGood} here~~
}
  \item{matrix}{
%%     ~~Describe \code{matrix} here~~
}
}
\details{
The probes can either be identified as control probes in the microarray annotation package, or alternatively the function will search the genome (given an appropriate BSgenome object) for suitable probes. 
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
Martin Aryee <aryee@jhu.edu>
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
function (dat, controlProbes = c("CONTROL_PROBES", "CONTROL_REGIONS"), 
    noCpGWindow = NULL, subject, onlyGood = FALSE, matrix = TRUE) 
{
    if (is.null(noCpGWindow)) {
        controlIndex <- which(getContainer(dat) \%in\% controlProbes)
    }
    else {
        if (class(subject) != "BSgenome") 
            stop("You must supply a BSgenome object for subject if using the noCpGWindow option. E.g. Hsapiens from the BSgenome.Hsapiens.UCSC.hg18 package.\n")
        chr <- pmChr(dat)
        pos <- pmPosition(dat)
        cpgd <- cpgdensity(subject, chr = chr, pos = pos, windowSize = noCpGWindow, 
            showProgress = FALSE)
        controlIndex <- which(cpgd == 0)
    }
    if (onlyGood == FALSE) {
        ret <- controlIndex
    }
    else {
        datPm <- log2(pm(dat)[controlIndex, , ])
        M <- datPm[, , "channel1"] - datPm[, , "channel2"]
        A <- (datPm[, , "channel1"] + datPm[, , "channel2"])/2
        samps <- sampleNames(dat)
        if (matrix) {
            ret <- sapply(samps, function(samp) {
                m <- M[, samp]
                a <- A[, samp]
                breakpoint <- optimize(checkBreakpoint, range(a), 
                  ms = m, as = a)$minimum
                a > breakpoint
            })
        }
        else {
            ret <- lapply(samps, function(samp) {
                m <- M[, samp]
                a <- A[, samp]
                breakpoint <- optimize(checkBreakpoint, range(a), 
                  ms = m, as = a)$minimum
                keep <- a > breakpoint
                controlIndex[keep]
            })
            names(ret) <- samps
        }
    }
    return(ret)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
