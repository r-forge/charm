\name{normalizeBetweenSamples}
\alias{normalizeBetweenSamples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Between-sample normalization
}
\description{
Between-sample normalization for two-color DNA methylation microarray data.
}
\usage{
normalizeBetweenSamples(dat, m = "allQuantiles", untreated = "none", enriched = "none", mapping = "affinity", numSegments = 3, controlProbes = "CONTROL_REGIONS", controlIndex = NULL, cluster = NULL, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
%%     ~~Describe \code{dat} here~~
}
  \item{m}{
%%     ~~Describe \code{m} here~~
}
  \item{untreated}{
%%     ~~Describe \code{untreated} here~~
}
  \item{enriched}{
%%     ~~Describe \code{enriched} here~~
}
  \item{mapping}{
%%     ~~Describe \code{mapping} here~~
}
  \item{numSegments}{
%%     ~~Describe \code{numSegments} here~~
}
  \item{controlProbes}{
%%     ~~Describe \code{controlProbes} here~~
}
  \item{controlIndex}{
%%     ~~Describe \code{controlIndex} here~~
}
  \item{cluster}{
%%     ~~Describe \code{cluster} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
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
function (dat, m = "allQuantiles", untreated = "none", enriched = "none", 
    mapping = "affinity", numSegments = 3, controlProbes = "CONTROL_REGIONS", 
    controlIndex = NULL, cluster = NULL, verbose = FALSE) 
{
    pms <- pm(dat)
    if (class(dat) == "OfflineTilingFeatureSet2") {
        M <- getPmM(dat)
    }
    else {
        M <- getM(dat)[pmindex(dat), , drop = FALSE]
    }
    if (m != "none") {
        M <- normQuantile(M, method = m)
        for (i in 1:ncol(M)) {
            pms[, i, "channel2"] <- 2^(log2(pms[, i, "channel1"]) - 
                M[, i])
        }
    }
    if (untreated == "complete") {
        if (any(class(pms) == "ff")) {
            m <- ffrowapply(rowMedians(log2(pms[i1:i2, , "channel1"]), 
                na.rm = TRUE), X = pms, RETURN = TRUE, RETCOL = NULL, 
                BATCHSIZE = 5000)[]
        }
        else {
            m <- rowMedians(log2(pms[, , "channel1"]))
        }
        for (i in 1:ncol(pms)) {
            pms[, i, "channel1"] <- 2^m
        }
    }
    else {
        pms[, , "channel1"] <- 2^(normQuantile(log2(pms[, , "channel1"]), 
            method = untreated))
    }
    for (i in 1:ncol(M)) {
        pms[, i, "channel2"] <- 2^(log2(pms[, i, "channel1"]) - 
            M[, i])
    }
    pm(dat) <- pms
    if (enriched == "sqn") {
        if (is.null(controlIndex)) 
            controlIndex <- getControlIndex(dat, controlProbes = controlProbes)
        if (any(class(pms) == "ff")) {
            pms <- SQN.ff(pms, channel = 2, ctrl.id = controlIndex, 
                cluster = cluster)
        }
        else {
            if (is.null(cluster)) {
                pms[, , "channel2"] <- 2^(SQN(y = log2(pms[, 
                  , "channel2"]), ctrl.id = controlIndex))
            }
            else {
                pms[, , "channel2"] <- 2^(parSQN(y = log2(pms[, 
                  , "channel2"]), ctrl.id = controlIndex, cluster = cluster))
            }
        }
    }
    else {
        pms[, , "channel2"] <- 2^(normQuantile(log2(pms[, , "channel2"]), 
            method = enriched))
    }
    pm(dat) <- pms
    return(dat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
