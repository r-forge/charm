\name{methp}
\alias{methp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate DNA methylation 
}
\description{
Estimate DNA methylation from McrBC/CHARM microarray data in terms of log-ratios or percentages.
}
\usage{
methp(dat, spatial = TRUE, spatialMethod = "kernel", spatial1d = NULL, spatial2d = NULL, bgSubtract = TRUE, withinSampleNorm = "loess", binSize = 500, numSegments = 3, useTot = TRUE, betweenSampleNorm = "quantile", scale = c(0.99, 0.99), minQCScore = NULL, controlProbes = c("CONTROL_PROBES", "CONTROL_REGIONS"), controlIndex = NULL, ctrlNoCpGWindow = NULL, subject = NULL, commonMethPercentParams = NULL, verbose = TRUE, cluster = NULL, returnM = FALSE, plotDensity = NULL, plotDensityCols = NULL, duplicate = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
a TilingFeatureSet object
}
  \item{spatial}{
boolean indicating whether to correct spatial artefacts
}
  \item{spatialMethod}{
method to use for spatial artefact correction
}
  \item{spatial1d}{
deprecated
}
  \item{spatial2d}{
deprecated
}
  \item{bgSubtract}{
boolean indicating whether to estimate and remove background signal before computing log-ratios
}
  \item{withinSampleNorm}{
within-sample normalization method. Choices are "loess" and "none"
}
  \item{binSize}{
%%     ~~Describe \code{binSize} here~~
}
  \item{numSegments}{
%%     ~~Describe \code{numSegments} here~~
}
  \item{useTot}{
%%     ~~Describe \code{useTot} here~~
}
  \item{betweenSampleNorm}{
%%     ~~Describe \code{betweenSampleNorm} here~~
}
  \item{scale}{
%%     ~~Describe \code{scale} here~~
}
  \item{minQCScore}{
%%     ~~Describe \code{minQCScore} here~~
}
  \item{controlProbes}{
%%     ~~Describe \code{controlProbes} here~~
}
  \item{controlIndex}{
%%     ~~Describe \code{controlIndex} here~~
}
  \item{ctrlNoCpGWindow}{
%%     ~~Describe \code{ctrlNoCpGWindow} here~~
}
  \item{subject}{
%%     ~~Describe \code{subject} here~~
}
  \item{commonMethPercentParams}{
%%     ~~Describe \code{commonMethPercentParams} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{cluster}{
%%     ~~Describe \code{cluster} here~~
}
  \item{returnM}{
%%     ~~Describe \code{returnM} here~~
}
  \item{plotDensity}{
%%     ~~Describe \code{plotDensity} here~~
}
  \item{plotDensityCols}{
%%     ~~Describe \code{plotDensityCols} here~~
}
  \item{duplicate}{
%%     ~~Describe \code{duplicate} here~~
}
}
\details{
This function provides probe-level estimates of percentage DNA methylation from CHARM microarray data.

}
\value{
A matrix of probe-level percentage methylation estimates, one column per sample.
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
function (dat, spatial = TRUE, spatialMethod = "kernel", spatial1d = NULL, 
    spatial2d = NULL, bgSubtract = TRUE, withinSampleNorm = "loess", 
    binSize = 500, numSegments = 3, useTot = TRUE, betweenSampleNorm = "quantile", 
    scale = c(0.99, 0.99), minQCScore = NULL, controlProbes = c("CONTROL_PROBES", 
        "CONTROL_REGIONS"), controlIndex = NULL, ctrlNoCpGWindow = NULL, 
    subject = NULL, commonMethPercentParams = NULL, verbose = TRUE, 
    cluster = NULL, returnM = FALSE, plotDensity = NULL, plotDensityCols = NULL, 
    duplicate = TRUE) 
{
    if (!is.null(spatial1d) | !is.null(spatial2d)) {
        warning("The spatial1d and spatial2d options are deprecated. Please use spatial instead\n")
    }
    if (!is.null(minQCScore)) {
        warning("methp no longer filters out probes based on quality score. Please use the pmQuality function to calculate probe quality.\n")
    }
    if (class(dat) == "OfflineTilingFeatureSet2" & duplicate) {
        dat <- clone(dat)
    }
    if (is.null(cluster)) {
        cl <- NULL
    }
    else if (is.numeric(cluster)) {
        cl <- charmCluster(cluster)
    }
    else if (any(class(cluster) == "cluster")) {
        cl <- cluster
    }
    if (!is.null(plotDensity)) {
        pdf(file = plotDensity, height = 11, width = 8)
        par(mfrow = c(5, 2), mar = c(2, 2, 4, 2))
        lwd <- rep(1, ncol(dat))
        if (is.null(plotDensityCols)) {
            if (!is.null(pData(rawData)$type)) {
                cols <- as.numeric(factor(pData(rawData)$type))
            }
            else {
                cols <- rep(1, ncol(dat))
            }
        }
        else {
            cols <- plotDensityCols
        }
        plotDensity(dat, main = "1. Raw", cols = cols, lwd = lwd)
    }
    if (is.list(betweenSampleNorm)) {
        bs <- betweenSampleNorm
    }
    else {
        if (betweenSampleNorm == "quantile") {
            bs <- list(m = "allQuantiles", untreated = "none", 
                enriched = "none")
            if (is.null(commonMethPercentParams)) 
                commonMethPercentParams <- TRUE
        }
        else if (betweenSampleNorm == "sqn") {
            bs <- list(m = "none", untreated = "complete", enriched = "sqn")
            if (is.null(commonMethPercentParams)) 
                commonMethPercentParams <- FALSE
        }
        else if (betweenSampleNorm == "sqn99") {
            bs <- list(m = "99", untreated = "complete", enriched = "sqn")
            if (is.null(commonMethPercentParams)) 
                commonMethPercentParams <- FALSE
        }
        else if (betweenSampleNorm == "none") {
            bs <- list(m = "none", untreated = "none", enriched = "none")
            if (is.null(commonMethPercentParams)) 
                commonMethPercentParams <- FALSE
        }
    }
    if (spatial) {
        if (verbose) 
            cat("Spatial normalization ")
        if (spatialMethod == "kernel") {
            cat("\n")
            dat <- spatialAdjust(dat, cluster = cl)
        }
        else if (spatialMethod == "poly") {
            cat("(polynomial surface)\n")
            dat <- spatialAdjust.poly(dat, cluster = cl)
        }
        else {
            stop("Invalid value for spatialMethod.\n")
        }
    }
    if (bgSubtract) {
        if (verbose) 
            cat("Background removal\n")
        dat <- bgAdjustBgp(dat, cluster = cl)
    }
    if (!is.null(plotDensity)) {
        plotDensity(dat, main = "2. After spatial & bg", cols = cols, 
            lwd = lwd)
    }
    if (is.null(controlIndex)) {
        controlIndex <- getControlIndex(dat, noCpGWindow = ctrlNoCpGWindow, 
            subject = subject, controlProbes = controlProbes)
    }
    if (verbose) 
        cat("Within sample normalization: ", withinSampleNorm, 
            "\n", sep = "")
    dat <- normalizeWithinSamples(dat, method = withinSampleNorm, 
        useTot = useTot, binSize = binSize, numSegments = numSegments, 
        cluster = cl, controlIndex = controlIndex, verbose = verbose)
    if (!is.null(plotDensity)) {
        plotDensity(dat, main = "3. After within-sample norm", 
            cols = cols, lwd = lwd)
    }
    if (verbose) {
        cat("Between sample normalization")
        if (is.list(betweenSampleNorm)) {
            cat(". M: ", bs$m, ", Untreated channel: ", bs$untreated, 
                ", Methyl-depleted channel: ", bs$enriched, "\n", 
                sep = "")
        }
        else {
            cat(": ", betweenSampleNorm, "\n", sep = "")
        }
    }
    dat <- normalizeBetweenSamples(dat, m = bs$m, untreated = bs$untreated, 
        enriched = bs$enriched, controlProbes = controlProbes, 
        controlIndex = controlIndex, cluster = cl, verbose = verbose)
    if (any(class(dat) == "OfflineTilingFeatureSet2")) {
        M <- getPmM(dat)
    }
    else {
        M <- getM(dat)[pmindex(dat), ]
    }
    if (msqn) {
        cat("The msqn option is not implemented\n")
    }
    if (!is.null(plotDensity)) {
        plotDensity(dat, main = "4. After between-sample norm", 
            cols = cols, lwd = lwd)
    }
    if (returnM == "TRUE" | returnM == "+") {
        retval <- M
    }
    else if (returnM == "-") {
        if (any(class(dat) == "OfflineTilingFeatureSet2")) {
            for (i in 1:ncol(M)) {
                M[, i] <- -M[, i]
                retval <- M
            }
        }
        else {
            retval <- -M
        }
    }
    else {
        if (verbose) 
            cat("Estimating percentage methylation\n")
        retval <- methPercent(m = M, commonParams = commonMethPercentParams, 
            ngc = countGC(dat), cluster = cl)
    }
    if (!is.null(plotDensity)) {
        if (returnM == "FALSE") 
            plotDensity(retval, main = "5. Percentage methylation", 
                controlIndex = getControlIndex(dat), rx = c(0, 
                  1), cols = cols, lwd = lwd)
        dev.off()
    }
    if (all(class(cluster) != "cluster")) 
        stopCluster(cl)
    return(retval)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
