\name{methPercent}
\alias{methPercent}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate percentage DNA methylation from log-ratios
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
methPercent(m, ngc, commonParams = TRUE, cluster = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{
%%     ~~Describe \code{m} here~~
}
  \item{ngc}{
%%     ~~Describe \code{ngc} here~~
}
  \item{commonParams}{
%%     ~~Describe \code{commonParams} here~~
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
function (m, ngc, commonParams = TRUE, cluster = NULL) 
{
    m <- as.matrix(m)
    param <- t(sapply(1:ncol(m), function(i) logmethParameters(m[, 
        i], ngc)))
    alpha <- unlist(param[, "alpha"])
    sigma <- unlist(param[, "sigma"])
    if (commonParams) {
        alpha[] <- median(alpha)
        sigma[] <- median(sigma)
    }
    if (any(class(m) == "ff")) {
        cl <- charmCluster(cluster)
        tmp <- ff(dim = dim(m), vmode = "double", finalizer = "close")
        parSapply(cl, 1:ncol(m), function(i, m, ngc, alpha, sigma, 
            tmp) {
            x <- m[, i]
            a <- alpha[i]
            b <- sigma[i]
            mu <- 0
            f0 <- dnorm(x, mean = 0, sd = b)
            f1 <- a * exp(0.5 * a^2 * b^2 - a * x) * (pnorm((x - 
                a * b^2)/b))
            p0.prior <- sum(x < 0, na.rm = TRUE)/sum(!is.na(x))
            p0 <- (p0.prior * f0)/((p0.prior * f0) + ((1 - p0.prior) * 
                f1))
            p1 <- 1 - p0
            x <- x - mu - a * b^2
            postM <- p1 * (x + b * ((1/sqrt(2 * pi)) * exp((-1/2) * 
                ((x/b)^2)))/pnorm(x/b))
            tmp[, i] <- 1 - 2^(-postM)
            return(NULL)
        }, m, ngc, alpha, sigma, tmp)
        ret <- tmp
        colnames(ret) <- colnames(m)
        if (all(class(cluster) != "cluster")) 
            stopCluster(cl)
    }
    else {
        ret <- sapply(1:ncol(m), function(i) {
            x <- m[, i]
            a <- alpha[i]
            b <- sigma[i]
            mu <- 0
            f0 <- dnorm(x, mean = 0, sd = b)
            f1 <- a * exp(0.5 * a^2 * b^2 - a * x) * (pnorm((x - 
                a * b^2)/b))
            p0.prior <- sum(x < 0, na.rm = TRUE)/sum(!is.na(x))
            p0 <- (p0.prior * f0)/((p0.prior * f0) + ((1 - p0.prior) * 
                f1))
            p1 <- 1 - p0
            x <- x - mu - a * b^2
            postM <- p1 * (x + b * ((1/sqrt(2 * pi)) * exp((-1/2) * 
                ((x/b)^2)))/pnorm(x/b))
            1 - 2^(-postM)
        })
        colnames(ret) <- colnames(m)
    }
    return(ret)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
