\name{plotDensity}
\alias{plotDensity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Log-ratio density plot for signal and control probes
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotDensity(dat, rx = c(-4, 6), controlIndex = NULL, pdfFile = NULL, main = NULL, cols = NULL, lwd = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
%%     ~~Describe \code{dat} here~~
}
  \item{rx}{
%%     ~~Describe \code{rx} here~~
}
  \item{controlIndex}{
%%     ~~Describe \code{controlIndex} here~~
}
  \item{pdfFile}{
%%     ~~Describe \code{pdfFile} here~~
}
  \item{main}{
%%     ~~Describe \code{main} here~~
}
  \item{cols}{
%%     ~~Describe \code{cols} here~~
}
  \item{lwd}{
%%     ~~Describe \code{lwd} here~~
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
function (dat, rx = c(-4, 6), controlIndex = NULL, pdfFile = NULL, 
    main = NULL, cols = NULL, lwd = NULL) 
{
    if (!is.null(pdfFile)) {
        pdf(pdfFile)
        par(mfcol = c(2, 1))
    }
    if (any(class(dat) == "OfflineTilingFeatureSet2")) {
        M <- getPmM(dat)
        lab <- sampleNames(dat)
    }
    else if (any(class(dat) == "TilingFeatureSet2")) {
        M <- getM(dat)[pmindex(dat), , drop = FALSE]
        lab <- sampleNames(dat)
    }
    else {
        M <- dat[, ]
        lab <- colnames(dat)
    }
    if (is.null(cols)) 
        cols <- 1:ncol(M)
    if (is.null(lwd)) 
        lwd <- rep(1, ncol(M))
    if (is.null(controlIndex)) 
        controlIndex <- getControlIndex(dat)
    plotDensityMat(M[, ], xlab = "M", lab = lab, main = paste(main, 
        "\nAll probes"), rx = rx, cols = cols, lwd = lwd)
    plotDensityMat(M[controlIndex, ], xlab = "M", lab = lab, 
        main = paste(main, "\nControl probes"), rx = rx, cols = cols, 
        lwd = lwd)
    if (!is.null(pdfFile)) 
        dev.off()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
