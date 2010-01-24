\name{pmQuality}
\alias{pmQuality}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Probe quality metrics
}
\description{
Calculate probe quality scores 
}
\usage{
pmQuality(dat, channel = "channel1", verbose = FALSE, idx = NULL, cluster = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
%%     ~~Describe \code{dat} here~~
}
  \item{channel}{
%%     ~~Describe \code{channel} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{idx}{
%%     ~~Describe \code{idx} here~~
}
  \item{cluster}{
%%     ~~Describe \code{cluster} here~~
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
function (dat, channel = "channel1", verbose = FALSE, idx = NULL, 
    cluster = NULL) 
{
    if (is.null(idx)) 
        idx <- 1:nrow(pm(dat))
    Ngc <- countGC(dat, "pm")[idx]
    bgNgc <- countGC(dat, "bg")
    if (class(dat) == "OfflineTilingFeatureSet2") {
        pms <- pm(dat)
        bgs <- bg(dat)
        cl <- charmCluster(cluster)
        pmq <- ff(dim = c(length(idx), ncol(pms)), vmode = "double", 
            finalizer = "close")
        parSapply(cl, 1:ncol(pms), function(i, pms, bgs, Ngc, 
            bgNgc) {
            if (verbose) 
                cat(".")
            fn <- tapply(bgs[, i, channel], bgNgc, ecdf)
            ret <- rep(NA, length(Ngc))
            for (ngc in unique(Ngc)) {
                idx2 <- Ngc == ngc
                closestIdx <- order(abs(as.numeric(names(fn)) - 
                  ngc))[1]
                bgngc <- as.character(names(fn)[closestIdx])
                ret[idx2] <- 100 * fn[[bgngc]](pms[idx, i, channel][idx2])
            }
            pmq[, i] <- ret
            return(NULL)
        }, pms, bgs, Ngc, bgNgc)
        if (all(class(cluster) != "cluster")) 
            stopCluster(cl)
    }
    else {
        pms <- pm(dat)[idx, , ]
        bgs <- bg(dat)
        pmq <- sapply(1:ncol(pms), function(i) {
            if (verbose) 
                cat(".")
            fn <- tapply(bgs[, i, channel], bgNgc, ecdf)
            ret <- rep(NA, length(Ngc))
            for (ngc in unique(Ngc)) {
                idx <- Ngc == ngc
                closestIdx <- order(abs(as.numeric(names(fn)) - 
                  ngc))[1]
                bgngc <- as.character(names(fn)[closestIdx])
                ret[idx] <- 100 * fn[[bgngc]](pms[idx, i, channel])
            }
            return(ret)
        })
    }
    colnames(pmq) <- sampleNames(dat)
    return(pmq)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
