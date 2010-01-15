\name{normalizeWithinSamples}
\alias{normalizeWithinSamples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Within-sample (between-channel) normalization for two-color data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
normalizeWithinSamples(dat, method = "loess", controlProbes = "CONTROL_REGIONS", controlIndex = NULL, numSegments = 3, bins = NULL, binSize = 500, approx = TRUE, breaks = 1000, useTot = TRUE, cluster = NULL, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
%%     ~~Describe \code{dat} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{controlProbes}{
%%     ~~Describe \code{controlProbes} here~~
}
  \item{controlIndex}{
%%     ~~Describe \code{controlIndex} here~~
}
  \item{numSegments}{
%%     ~~Describe \code{numSegments} here~~
}
  \item{bins}{
%%     ~~Describe \code{bins} here~~
}
  \item{binSize}{
%%     ~~Describe \code{binSize} here~~
}
  \item{approx}{
%%     ~~Describe \code{approx} here~~
}
  \item{breaks}{
%%     ~~Describe \code{breaks} here~~
}
  \item{useTot}{
%%     ~~Describe \code{useTot} here~~
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
function (dat, method = "loess", controlProbes = "CONTROL_REGIONS", 
    controlIndex = NULL, numSegments = 3, bins = NULL, binSize = 500, 
    approx = TRUE, breaks = 1000, useTot = TRUE, cluster = NULL, 
    verbose = FALSE) 
{
    if (grepl("gc", method)) {
        if (class(dat) == "OfflineTilingFeatureSet2") {
            stop("normalizeWithinSamples does not yet support the method='gc' option for OfflineTilingFeatureSet2 objects. Please use the method='affinity' option.")
        }
        mAdj <- diffAmpEstGC(dat, method = method, controlProbes = controlProbes)
        dat <- diffAmpAdjustGC(dat, mAdj)
    }
    if (grepl("loess", method)) {
        dat <- normalizeLoess(dat, controlIndex = controlIndex, 
            controlProbes = controlProbes, approx = approx, breaks = breaks, 
            cluster = cluster)
    }
    if (grepl("affinity", method)) {
        if (is.null(controlIndex)) {
            controlIndex <- getControlIndex(dat, controlProbes = controlProbes)
        }
        if (is.null(bins)) {
            bins <- affinityBin(dat, controlIndex = controlIndex, 
                binSize = binSize, useTot = useTot, verbose = verbose)
        }
        pms <- pm(dat)
        for (i in 1:ncol(pms)) {
            mAdj <- bins$binMed[bins$bin, i]
            pms[, i, "channel2"] <- 2^(log2(pms[, i, "channel2"]) + 
                mAdj)
        }
        pm(dat) <- pms
    }
    if (grepl("median", method)) {
        if (class(dat) == "OfflineTilingFeatureSet2") {
            stop("normalizeWithinSamples does not yet support the method='median' option for OfflineTilingFeatureSet2 objects. Please use the method='affinity' option.")
        }
        if (is.null(controlIndex)) {
            controlIndex <- getControlIndex(dat, controlProbes = controlProbes)
        }
        datPm <- log2(pm(dat))
        M <- datPm[, , "channel1"] - datPm[, , "channel2"]
        mCtrl <- M[controlIndex, ]
        mAdj <- apply(mCtrl, 2, median, na.rm = TRUE)
        datPm[, , "channel2"] <- sweep(datPm[, , "channel2"], 
            2, mAdj, FUN = "+")
        pm(dat) <- 2^datPm
    }
    return(dat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
