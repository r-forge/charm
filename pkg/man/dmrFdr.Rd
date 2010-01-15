\name{dmrFdr}
\alias{dmrFdr}
\title{
Calculate FDR q-values for differentially methylated regions (DMRs)
}
\description{
Estimate false discovery rate q-values for a set of differentially methylated regions using a permutation approach.
}
\usage{
dmrFdr(dmr, compare = 1, numPerms = 1000, seed = NULL, verbose = TRUE)
}
\arguments{
  \item{dmr}{
a dmr object as returned by \code{\link{dmrFinder}}
}
  \item{compare}{
The dmr table for which to calculate DMRs. See details.
}
  \item{numPerms}{
Number of permutations
}
  \item{seed}{
Random seed (for reproducibility)
}
  \item{verbose}{
Boolean
}
}
\details{
This function estimates false discovery rate q-values for a dmr object returned by \code{\link{dmrFinder}}. dmrFinder can return a set of DMR tables with one or more pair-wise comparisons between groups. dmrFdr currently only calculated q-values for one of these at a time. The dmr table to use is specified by the compare option.
}
\value{
a list object in the same format as the input, but with extra p-val and q-val columns for the tabs element.
}
\references{

}
\author{
Martin Aryee <aryee@jhu.edu>
}


\seealso{
\code{\link{dmrFinder}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (dmr, compare = 1, numPerms = 1000, seed = NULL, verbose = TRUE) 
{
    if (length(compare) != 1) 
        stop("You must choose one comparison at a time when calculating FDRs. Please set dmr to be one of: ", 
            paste(names(dmr$tabs), collapse = ", "), "\n")
    if (is.numeric(compare)) 
        compare <- names(dmr$tabs)[compare]
    cat("Calculating q-values for DMRs between", compare, "\n")
    pdInfo = get(dmr$package)
    class(pdInfo) = "TilingFeatureSet"
    chr = pmChr(pdInfo)
    pos = pmPosition(pdInfo)
    o <- order(chr, pos)
    p <- matrix(NA, nrow = nrow(dmr$p), ncol = ncol(dmr$p))
    p[o, ] <- dmr$p
    keep <- dmr$groups \%in\% unlist(strsplit(compare, "-"))
    p <- p[, keep]
    n <- sum(keep)
    n1 <- sum(dmr$groups == unlist(strsplit(compare, "-"))[1])
    maxPerms <- choose(n, n1)
    if (numPerms == "all") 
        numPerms <- maxPerms
    if (numPerms > maxPerms) {
        cat("Given the sample sizes in the two groups the maximum number of permutations is ", 
            maxPerms, ".\n")
        numPerms <- maxPerms
    }
    if (!is.null(seed)) 
        set.seed(seed)
    s <- sample(1:maxPerms, numPerms)
    grp1 <- combinations(n, n1)[s, ]
    if (verbose) 
        cat("Finding permuted data DMRs. Estimating time remaining\n")
    areas <- lapply(1:numPerms, function(i) {
        groups <- rep("grp2", n)
        groups[grp1[i, ]] <- "grp1"
        st <- system.time(dmrPerm <- dmrFinder(dmr$package, p = p, 
            groups = groups, cutoff = dmr$cutoff, filter = dmr$filter, 
            ws = dmr$ws, verbose = FALSE))[3]
        if (verbose & (i \%in\% round(seq(1, numPerms, length.out = 10)))) {
            cat(i, "/", numPerms, " (", prettyTime((numPerms - 
                i) * st), " remaining)\n", sep = "")
        }
        dmrPerm$tabs[[1]]$area
    })
    nullDist <- unlist(areas)
    fn <- ecdf(nullDist)
    pval <- 1 - fn(dmr$tabs[[compare]]$area)
    pi0 <- pi0.est(pval)$p0
    qval <- qvalue.cal(pval, pi0)
    dmr$tabs[[compare]] <- cbind(dmr$tabs[[compare]], pval, qval)
    return(dmr)
  }
}
