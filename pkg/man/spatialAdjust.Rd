\name{spatialAdjust}
\alias{spatialAdjust}
\title{
Correct spatial artifacts
}
\description{
Remove spatial artifacts from microarray data stored in TilingFeatureSet objects
}
\usage{
spatialAdjust(dat, cluster = NULL, blockSize = 50, theta = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
TilingFeatureSet
}
  \item{cluster}{
%%     ~~Describe \code{cluster} here~~
}
  \item{blockSize}{
The dimension of blocks to use when averaging pixel values
}
  \item{theta}{
smoothing parameter
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
function (dat, cluster = NULL, blockSize = 50, theta = 1) 
{
    x <- getX(dat, "pm")
    y <- getY(dat, "pm")
    bgX <- getX(dat, "bg")
    bgY <- getY(dat, "bg")
    NCOL = dim(dat)["Samples"]
    pms <- pm(dat)
    bgs <- bg(dat)
    nx <- max(x)/blockSize
    ny <- max(y)/blockSize
    d <- discretize.image(cbind(y, x), m = ny, n = nx)
    dBg <- discretize.image(cbind(bgY, bgX), grid = d$grid)
    if (class(dat) == "OfflineTilingFeatureSet2") {
        cl <- charmCluster(cluster)
        tmp <- parLapply(cl, 1:NCOL, function(i, pms, bgs, x, 
            y, bgX, bgY, nx, ny, d, dBg) {
            tot <- log2(pms[, i, "channel1"])
            bgTot <- log2(bgs[, i, "channel1"])
            sm <- spatialAdjustVec(tot, d)
            pms[, i, "channel1"] <- 2^sm$zAdj
            smBg <- spatialAdjustVec(bgTot, dBg, sm$ims)
            bgs[, i, "channel1"] <- 2^smBg$zAdj
            en <- log2(pms[, i, "channel2"])
            bgEn <- log2(bgs[, i, "channel2"])
            sm <- spatialAdjustVec(en, d)
            pms[, i, "channel2"] <- 2^sm$zAdj
            smBg <- spatialAdjustVec(bgEn, dBg, sm$ims)
            bgs[, i, "channel2"] <- 2^smBg$zAdj
            return(NULL)
        }, pms, bgs, x, y, bgX, bgY, nx, ny, d, dBg)
        if (all(class(cluster) != "cluster")) 
            stopCluster(cl)
    }
    else {
        for (i in 1:ncol(pms)) {
            for (channel in 1:2) {
                tot <- log2(pms[, i, channel])
                bgTot <- log2(bgs[, i, channel])
                sm <- spatialAdjustVec(tot, d)
                pms[, i, channel] <- 2^sm$zAdj
                smBg <- spatialAdjustVec(bgTot, dBg, sm$ims)
                bgs[, i, channel] <- 2^smBg$zAdj
            }
        }
        pm(dat) <- pms
        bg(dat) <- bgs
    }
    return(dat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
