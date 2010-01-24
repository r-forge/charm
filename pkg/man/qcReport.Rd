\name{qcReport}
\alias{qcReport}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Microarray quality report
}
\description{
Calculate microarray quality scores and produce an optional pdf report}
\usage{
qcReport(dat, file = NULL, utRange = c(30, 100), enRange = c(8, 12), numProbes = 5e+05, cluster = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
%%     ~~Describe \code{dat} here~~
}
  \item{file}{
%%     ~~Describe \code{file} here~~
}
  \item{utRange}{
%%     ~~Describe \code{utRange} here~~
}
  \item{enRange}{
%%     ~~Describe \code{enRange} here~~
}
  \item{numProbes}{
%%     ~~Describe \code{numProbes} here~~
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
function (dat, file = NULL, utRange = c(30, 100), enRange = c(8, 
    12), numProbes = 5e+05, cluster = NULL) 
{
    n <- nrow(pm(dat))
    if (numProbes == -1) 
        numProbes <- n
    if (n < numProbes) 
        numProbes <- n
    idx <- as.integer(seq(1, n, length.out = numProbes))
    pmQual <- pmQuality(dat, idx = idx, cluster = cluster)
    X <- getX(dat, "pm")[idx]
    Y <- getY(dat, "pm")[idx]
    imgs1 <- arrayImage(X, Y, pmQual, cluster = cluster)
    if (class(dat) == "OfflineTilingFeatureSet2") {
        pms <- pm(dat)
        c1 <- ff(as.vector(log2(pms[idx, , "channel1"])), dim = c(length(idx), 
            ncol(pms)), finalizer = "close")
        colnames(c1) <- colnames(pms)
        tmp <- arrayImage(X, Y, c1, cluster = cluster)
        sd1 <- unlist(lapply(tmp, function(x) sd(as.vector(x$z), 
            na.rm = TRUE)))
        c2 <- ff(as.vector(log2(pms[idx, , "channel2"])), dim = c(length(idx), 
            ncol(pms)), finalizer = "close")
        colnames(c2) <- colnames(pms)
        imgs2 <- arrayImage(X, Y, c2, cluster = cluster)
        sd2 <- unlist(lapply(imgs2, function(x) sd(as.vector(x$z), 
            na.rm = TRUE)))
    }
    else {
        tmp <- arrayImage(X, Y, log2(pm(dat)[idx, , "channel1"]))
        sd1 <- unlist(lapply(tmp, function(x) sd(as.vector(x$z), 
            na.rm = TRUE)))
        imgs2 <- arrayImage(X, Y, log2(pm(dat)[idx, , "channel2"]))
        sd2 <- unlist(lapply(imgs2, function(x) sd(as.vector(x$z), 
            na.rm = TRUE)))
    }
    sdRange <- c(0, max(sd1, sd2) * 1.1)
    if (any(class(pmQual) == "ff")) {
        pmSignal <- sapply(1:ncol(pmQual), function(i) mean(pmQual[, 
            i]))
        names(pmSignal) <- colnames(pmQual)
    }
    else {
        pmSignal <- colMeans(pmQual, na.rm = TRUE)
    }
    if (!is.null(file)) {
        pdf(file = file, width = 8, height = 10.5)
        lay <- rbind(cbind(rep(1, 3), matrix(rep(2:4, each = 6), 
            ncol = 6)), rep(5, 7))
        layout(lay)
        n <- length(pmSignal)
        if (is.null(names(pmSignal))) {
            names(pmSignal) <- paste("Sample", 1:n)
        }
        l <- length(pmSignal)
        longestName <- max(sapply(names(pmSignal), nchar))
        cexh <- ifelse(n > 50, max(0.25, 1 - (n - 50) * 0.01), 
            1)
        cexw <- ifelse(longestName > 17, max(0.25, 1 - (longestName - 
            15) * 0.05), 1)
        cex <- min(cexh, cexw)
        par(oma = c(2, 2, 4, 2))
        o <- order(names(pmSignal))
        par(mar = c(5, 0, 3, 0))
        plot.new()
        plot.window(xlim = c(0, 10), ylim = c(0, l))
        text(1, l:1 - 0.5, names(pmSignal)[o], adj = 0.1, cex = cex)
        abline(h = 1:(l + 1) - 1, lwd = 0.5, lty = 3)
        xmin <- max(0, floor((min(pmSignal) - 5)/5) * 5)
        xmax <- min(100, ceiling((max(pmSignal) + 5)/5) * 5)
        par(mar = c(5, 1, 3, 1))
        plot(pmSignal[o], l:1, las = 1, ylab = "", xlab = "", 
            type = "p", pch = 19, xlim = c(xmin, xmax), ylim = c(0.5, 
                l + 0.5), yaxt = "n", frame.plot = TRUE, main = "Signal strength")
        abline(h = 1:(l + 1) - 0.5, lwd = 0.5, lty = 3)
        par(mar = c(5, 1, 3, 1))
        plot(sd1[o], l:1, las = 1, ylab = "", xlab = "", type = "p", 
            pch = 19, xlim = sdRange, ylim = c(0.5, l + 0.5), 
            yaxt = "n", frame.plot = TRUE, main = "Channel 1\nstandard deviation")
        abline(h = 1:(l + 1) - 0.5, lwd = 0.5, lty = 3)
        par(mar = c(5, 1, 3, 1))
        plot(sd2[o], l:1, las = 1, ylab = "", xlab = "", type = "p", 
            pch = 19, xlim = sdRange, ylim = c(0.5, l + 0.5), 
            yaxt = "n", frame.plot = TRUE, main = "Channel 2\nstandard deviation")
        abline(h = 1:(l + 1) - 0.5, lwd = 0.5, lty = 3)
        hist(pmSignal, main = "Signal strength histogram", xlab = "Signal strength")
        l <- layout(matrix(c(1, 1, 1, 1, 2:21), 6, 4, byrow = TRUE))
        plot.new()
        text(0.5, 0.2, "Untreated Channel: PM probe quality", 
            cex = 2)
        par(mar = c(4.5, 8, 3.5, 3))
        arrayPlot(imgs1[o], r = utRange)
        l <- layout(matrix(c(1, 1, 1, 1, 2:21), 6, 4, byrow = TRUE))
        plot.new()
        text(0.5, 0.2, "Enriched Channel: PM signal intensity", 
            cex = 2)
        par(mar = c(4.5, 8, 3.5, 3))
        arrayPlot(imgs2[o], r = enRange)
        dev.off()
    }
    return(as.data.frame(cbind(pmSignal, sd1, sd2)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
