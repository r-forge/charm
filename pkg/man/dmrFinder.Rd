\name{dmrFinder}
\alias{dmrFinder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find differentially methylated regions (DMRs)
}
\description{
Find differentially methylated regions (DMRs) from tiling microarray data.
}
\usage{
dmrFinder(eset = NULL, groups, p = NULL, l = NULL, chr = NULL, pos = NULL, pns = NULL, sdBins = NULL, controlIndex = NULL, Indexes = NULL, filter = NULL, package = NULL, ws = 7, verbose = TRUE, compare = "all", bgSubtract = TRUE, withinSampleNorm = "loess", binSize = 500, numSegments = 3, betweenSampleNorm = "quantile", minQCScore = NULL, controlProbes = c("CONTROL_PROBES", "CONTROL_REGIONS"), cluster = NULL, cutoff = 0.995, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{
%%     ~~Describe \code{eset} here~~
}
  \item{groups}{
%%     ~~Describe \code{groups} here~~
}
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{l}{
%%     ~~Describe \code{l} here~~
}
  \item{chr}{
%%     ~~Describe \code{chr} here~~
}
  \item{pos}{
%%     ~~Describe \code{pos} here~~
}
  \item{pns}{
%%     ~~Describe \code{pns} here~~
}
  \item{sdBins}{
%%     ~~Describe \code{sdBins} here~~
}
  \item{controlIndex}{
%%     ~~Describe \code{controlIndex} here~~
}
  \item{Indexes}{
%%     ~~Describe \code{Indexes} here~~
}
  \item{filter}{
%%     ~~Describe \code{filter} here~~
}
  \item{package}{
%%     ~~Describe \code{package} here~~
}
  \item{ws}{
%%     ~~Describe \code{ws} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{compare}{
%%     ~~Describe \code{compare} here~~
}
  \item{bgSubtract}{
%%     ~~Describe \code{bgSubtract} here~~
}
  \item{withinSampleNorm}{
%%     ~~Describe \code{withinSampleNorm} here~~
}
  \item{binSize}{
%%     ~~Describe \code{binSize} here~~
}
  \item{numSegments}{
%%     ~~Describe \code{numSegments} here~~
}
  \item{betweenSampleNorm}{
%%     ~~Describe \code{betweenSampleNorm} here~~
}
  \item{minQCScore}{
%%     ~~Describe \code{minQCScore} here~~
}
  \item{controlProbes}{
%%     ~~Describe \code{controlProbes} here~~
}
  \item{cluster}{
%%     ~~Describe \code{cluster} here~~
}
  \item{cutoff}{
%%     ~~Describe \code{cutoff} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (eset = NULL, groups, p = NULL, l = NULL, chr = NULL, 
    pos = NULL, pns = NULL, sdBins = NULL, controlIndex = NULL, 
    Indexes = NULL, filter = NULL, package = NULL, ws = 7, verbose = TRUE, 
    compare = "all", bgSubtract = TRUE, withinSampleNorm = "loess", 
    binSize = 500, numSegments = 3, betweenSampleNorm = "quantile", 
    minQCScore = NULL, controlProbes = c("CONTROL_PROBES", "CONTROL_REGIONS"), 
    cluster = NULL, cutoff = 0.995, ...) 
{
    groups = as.character(groups)
    if (identical(compare, "all")) 
        compare = comp(groups)
    if (length(compare)\%\%2 != 0) 
        stop("compare must have an even number of elements.")
    args = list(filter = filter, ws = ws, betweenSampleNorm = betweenSampleNorm, 
        withinSampleNorm = withinSampleNorm, minQCScore = minQCScore, 
        controlProbes = controlProbes, cutoff = cutoff)
    if (is.null(eset)) {
        if (is.null("p") & is.null("l")) 
            stop("p or l must be supplied.")
        args = c("chr", "pos", "pns", "controlIndex")
        nulls = sapply(args, function(x) is.null(get(x)))
        if (any(nulls)) 
            stop(paste("The following arguments are missing:", 
                paste(args[nulls], collapse = ", ")))
        lens = c(nrow(p), nrow(l), length(chr), length(pos), 
            length(pns))
        if (length(unique(lens)) != 1) 
            stop("p, l, chr, pos, and/or pns are incompatible.")
        stopifnot(length(groups) == max(ncol(p), ncol(l)))
        index <- which(!is.na(chr) & !is.na(pos) & !is.na(pns))
        index = index[order(chr[index], pos[index])]
        chr = chr[index]
        pos = pos[index]
        pns = pns[index]
        controlIndex = which(index \%in\% controlIndex)
        if (!is.null(sdBins)) 
            sdBins <- sdBins[index]
        if (!is.null(p)) 
            p = p[index, ]
        if (!is.null(l)) 
            l = l[index, ]
    }
    else if (is.character(eset)) {
        pdInfo = get(eset)
        class(pdInfo) = "TilingFeatureSet"
        chr = pmChr(pdInfo)
        pos = pmPosition(pdInfo)
        index = which(rowSums(is.na(p)) == 0)
        index = index[order(chr[index], pos[index])]
        chr = chr[index]
        pos = pos[index]
        pns = probeNames(pdInfo)[index]
        p = p[index, ]
        controlIndex <- which(getContainer(pdInfo) \%in\% controlProbes)
        controlIndex = which(index \%in\% controlIndex)
        if (!is.null(sdBins)) 
            sdBins <- sdBins[index]
        package = eset
        if (package == "pd.feinberg.mm8.me.hx1") {
        }
        if (package == "pd.feinberg.hg18.me.hx1") {
        }
    }
    else {
        stopifnot(length(groups) == length(eset))
        if (is.null(p) & is.null(l)) {
            p <- methp(eset, bgSubtract = bgSubtract, withinSampleNorm = withinSampleNorm, 
                binSize = binSize, numSegments = numSegments, 
                betweenSampleNorm = betweenSampleNorm, controlProbes = controlProbes, 
                cluster = cluster, verbose = TRUE)[, ]
        }
        chr = pmChr(eset)
        pos = pmPosition(eset)
        index = which(rowSums(is.na(p)) == 0)
        index = index[order(chr[index], pos[index])]
        chr = chr[index]
        pos = pos[index]
        pns = probeNames(eset)[index]
        p = p[index, ]
        controlIndex = getControlIndex(eset, controlProbes = controlProbes)
        controlIndex = which(index \%in\% controlIndex)
        if (!is.null(sdBins)) 
            sdBins <- sdBins[index]
        package = annotation(eset)
        if (package == "pd.feinberg.mm8.me.hx1") {
        }
        if (package == "pd.feinberg.hg18.me.hx1") {
        }
    }
    if (is.null(l)) {
        l = log(p) - log(1 - p)
    }
    Indexes = split(seq(along = pns), pns)
    tog = get.tog(l = l, groups = groups, compare = compare, 
        verbose = verbose)
    lm = tog$lm
    ls = tog$ls
    ns = tog$ns
    COMPS = tog$COMPS
    tt = get.tt(lm = lm, ls = ls, ns = ns, COMPS = COMPS, Indexes = Indexes, 
        filter = filter, ws = ws, verbose = verbose)
    res = vector("list", ncol(tt))
    names(res) = colnames(tt)
    if (verbose) 
        cat("Finding DMRs for each pairwise comparison.")
    for (r in 1:nrow(COMPS)) {
        j = COMPS[r, 1]
        k = COMPS[r, 2]
        if (verbose) 
            cat("\n", colnames(tt)[r])
        DF = ifelse(ns[j] == 1 & ns[k] == 1, 1, ns[j] + ns[k] - 
            2)
        if (length(sdBins) == 0) {
            K = mad(tt[, r], na.rm = TRUE) * qt(cutoff, DF)
        }
        else {
            s <- tapply(tt[, r], sdBins, mad, na.rm = TRUE)
            K = s[sdBins] * qt(cutoff, DF)
        }
        LAST = 0
        segmentation = vector("numeric", nrow(tt))
        type = vector("numeric", nrow(tt))
        for (i in seq(along = Indexes)) {
            if (verbose) 
                if (i\%\%1000 == 0) 
                  cat(".")
            Index = Indexes[[i]]
            y = tt[Index, r]
            if (length(sdBins) == 0) {
                tmp = sign(y) * as.numeric(abs(y) > K)
            }
            else {
                Ki <- K[Index]
                tmp = sign(y) * as.numeric(abs(y) > Ki)
            }
            tmp2 = cumsum(c(1, diff(tmp) != 0)) + LAST
            segmentation[Index] = tmp2
            type[Index] = tmp
            LAST = max(tmp2)
        }
        Index = which(type != 0)
        res[[r]] = data.frame(chr = tapply(chr[Index], segmentation[Index], 
            function(x) x[1]), start = tapply(pos[Index], segmentation[Index], 
            min), end = tapply(pos[Index], segmentation[Index], 
            max), p1 = NA, p2 = NA, regionName = tapply(pns[Index], 
            segmentation[Index], function(x) x[1]), indexStart = tapply(Index, 
            segmentation[Index], min), indexEnd = tapply(Index, 
            segmentation[Index], max))
        if (is.null(p)) {
            colnames(res[[r]]) <- sub("p1", "m1", colnames(res[[r]]))
            colnames(res[[r]]) <- sub("p2", "m2", colnames(res[[r]]))
            res[[r]]$m1 = tapply(lm[Index, j], segmentation[Index], 
                mean)
            res[[r]]$m2 = tapply(lm[Index, k], segmentation[Index], 
                mean)
            length = res[[r]]$indexEnd - res[[r]]$indexStart + 
                1
            area = abs(res[[r]]$m2 - res[[r]]$m1) * length
            res[[r]]$area = area
        }
        else {
            res[[r]]$p1 = tapply(1/(1 + exp(-lm[Index, j])), 
                segmentation[Index], mean)
            res[[r]]$p2 = tapply(1/(1 + exp(-lm[Index, k])), 
                segmentation[Index], mean)
            length = res[[r]]$indexEnd - res[[r]]$indexStart + 
                1
            area = abs(res[[r]]$p2 - res[[r]]$p1) * length
            res[[r]]$area = area
        }
        res[[r]] = res[[r]][order(-area), ]
    }
    if (verbose) 
        cat("\nDone\n")
    return(list(tabs = res, p = p, m = l, chr = chr, pos = pos, 
        pns = pns, index = index, controlIndex = controlIndex, 
        gm = lm, groups = groups, args = args, cutoff = cutoff, 
        filter = filter, ws = ws, comps = COMPS, package = package))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
