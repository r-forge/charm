
oligoVer <- as.numeric(sub("(.*\\..*)\\..*", "\\1", package.version("oligo")))
if (oligoVer>=1.11) {
	setClass("OfflineTilingFeatureSet2", contains="TilingFeatureSet")
} else {
	setClass("OfflineTilingFeatureSet2", contains="TilingFeatureSet2")
}

setMethod("initialize", "OfflineTilingFeatureSet2", function(.Object, ...) {
		callNextMethod(.Object, ...)
	})

setMethod("pm", "OfflineTilingFeatureSet2", function (object, subset=NULL, ...) {
	object@pmIntensities
})

setGeneric("clone", function(object) standardGeneric("clone"))

setMethod("clone", "OfflineTilingFeatureSet2", function (object) {
	object2 <- object
	slot(object2, "pmIntensities", check=FALSE) <- ff::clone(object@pmIntensities)
	slot(object2, "bgIntensities", check=FALSE) <- ff::clone(object@bgIntensities)
	object2
})


setGeneric("getPmM", function(object, rows, cols) standardGeneric("getPmM"))

getPmM.OfflineTilingFeatureSet2 <- function (object, rows, cols) {
	d <- dim(object@pmIntensities)
	if (missing(rows)) rows <- 1:d[1]
	if (missing(cols)) cols <- 1:d[2]
	M <- ff(dim=c(length(rows), length(cols)), vmode="double", finalizer="close")
	for (i in cols) {
		M[,i] <- log2(object@pmIntensities[rows,i,1])-log2(object@pmIntensities[rows,i,2])		
	}
	colnames(M) <-  sampleNames(object)
	M
}

setMethod("getPmM", "OfflineTilingFeatureSet2", getPmM.OfflineTilingFeatureSet2 )

setMethod("bg", "OfflineTilingFeatureSet2", function (object, subset=NULL) {
	object@bgIntensities
})

setMethod("pm<-", "OfflineTilingFeatureSet2", function (object, value) {
	object
})

setMethod("bg<-", "OfflineTilingFeatureSet2", function (object, value) {
	object
})

rowMedians.ff_matrix <- function(imat, na.rm=TRUE) {
	ret <- ffrowapply(rowMedians(imat[i1:i2,,drop=FALSE], na.rm=na.rm), X=imat, 
		RETURN=TRUE, RETCOL=NULL, BATCHSIZE=5000)	
	return(ret[])
}
setMethod("rowMedians", "ff_matrix", rowMedians.ff_matrix)


readXysMatrix2Offline <- function (filenames1, filenames2, pkgname) {
	pkg <- get(pkgname)
	pmIndex <- pmindex(pkg)
	bgIndex <- bgindex(pkg)
    pmIntensities <- NULL
    bgIntensities <- NULL
	#M <- NULL
    for (i in seq(along = filenames1)) {
        tmp <- oligo:::readonexysfile(filenames1[i])
		if (is.null(pmIntensities)) {
			pmIntensities <- ff(dim=c(length(pmIndex), length(filenames1), 2),
			 	vmode="double", finalizer="close")
			bgIntensities <- ff(dim=c(length(bgIndex), length(filenames1), 2),
			 	vmode="double", finalizer="close")
			#M <- ff(dim=c(length(pmIndex), length(filenames1)), vmode="double")
			idx <- tmp[["X"]] + (tmp[["Y"]] - 1) * geometry(pkg)[2]
			idx <- order(idx)
		} 
        pmIntensities[,i,1] <- tmp$SIGNAL[idx][pmIndex]
        bgIntensities[,i,1] <- tmp$SIGNAL[idx][bgIndex]
		rm(tmp)
        tmp <- oligo:::readonexysfile(filenames2[i])
        pmIntensities[,i,2] <- tmp$SIGNAL[idx][pmIndex]
        bgIntensities[,i,2] <- tmp$SIGNAL[idx][bgIndex]
		#M[,i] <- log2(tmp1$SIGNAL[idx][pmIndex]) - log2(tmp2$SIGNAL[idx][pmIndex])
    }
    #list(pm=pmIntensities, bg=bgIntensities, M=M)
	list(pm=pmIntensities, bg=bgIntensities)
}

read.xysfiles2.offline <- function (channel1, channel2, pkgname, phenoData, featureData, 
    experimentData, notes, verbose = TRUE, sampleNames, checkType = TRUE) 
{
    filenames <- c(channel1, channel2)
    oligo:::checkValidFilenames(filenames)
    if (checkType) 
        stopifnot(oligo:::checkChipTypes(filenames, verbose, "nimblegen"))
    firstline <- oligo:::readxysHeader(filenames[1])
    designname <- unlist(strsplit(firstline[grep("designname", 
        firstline, fixed = TRUE, useBytes = TRUE)], "="))[2]
    if (missing(pkgname)) 
        pkgname <- cleanPlatformName(designname)
    if (require(pkgname, character.only = TRUE)) {
        if (verbose) 
            message("Platform design info loaded.")
    }
    else {
        stop("Must install the ", pkgname, " package.")
    }
    arrayType <- kind(get(pkgname))
    #channel1Intensities <- readXysMatrixOffline(channel1, pkgname)
    #channel2Intensities <- readXysMatrixOffline(channel2, pkgname)
	intensities <- readXysMatrix2Offline(channel1, channel2, pkgname)
	emptyIntensities <- matrix(NA, nrow=0, ncol=length(sampleNames))
	colnames(emptyIntensities) <- sampleNames
    metadata <- oligo:::getMetadata(emptyIntensities, channel1, phenoData, 
        featureData, experimentData, notes, sampleNames)
    out <- new("OfflineTilingFeatureSet2", channel1 = emptyIntensities, 
        channel2 = emptyIntensities, manufacturer = "NimbleGen", 
        annotation = pkgname, phenoData = metadata[["phenoData"]], 
        experimentData = metadata[["experimentData"]])
	#slot(out, "channel1ff", check=FALSE) <- channel1Intensities
	#slot(out, "channel2ff", check=FALSE) <- channel2Intensities
	nam <- list(NULL, as.character(sampleNames), c("channel1", "channel2"))
	dimnames(intensities$pm) <- nam
	dimnames(intensities$bg) <- nam
	slot(out, "pmIntensities", check=FALSE) <- intensities$pm
	slot(out, "bgIntensities", check=FALSE) <- intensities$bg
	#slot(out, "M", check=FALSE) <- intensities$M
    return(out)
}

read.celfiles2.offline <- function (channel1, channel2, pkgname, 
	phenoData, featureData, 
    experimentData, notes, verbose = TRUE, sampleNames, rm.mask = FALSE, 
    rm.outliers = FALSE, rm.extra = FALSE, sd = FALSE, checkType = TRUE) 
{
    filenames <- c(channel1, channel2)
    oligo:::checkValidFilenames(filenames)
    if (checkType) 
        stopifnot(oligo:::checkChipTypes(filenames, verbose, "affymetrix", 
            useAffyio=FALSE))
    chiptype <- oligo:::getCelChipType(filenames[1], useAffyio=FALSE)
    if (missing(pkgname)) 
        pkgname <- cleanPlatformName(chiptype)
    if (require(pkgname, character.only = TRUE)) {
        if (verbose) 
            message("Platform design info loaded.")
    }
    else {
        stop("Must install the ", pkgname, " package.")
    }
    arrayType <- kind(get(pkgname))
	
	intensities <- readCelIntensities2.offline(channel1, channel2, pkgname)
	emptyIntensities <- matrix(NA, nrow=0, ncol=length(sampleNames))
	colnames(emptyIntensities) <- sampleNames
    metadata <- oligo:::getMetadata(emptyIntensities, channel1, phenoData, 
        featureData, experimentData, notes, sampleNames)
    out <- new("OfflineTilingFeatureSet2", channel1 = emptyIntensities, 
        channel2 = emptyIntensities, manufacturer = "Affymetrix", 
        annotation = pkgname, phenoData = metadata[["phenoData"]])
	nam <- list(NULL, as.character(sampleNames), c("channel1", "channel2"))
	dimnames(intensities$pm) <- nam
	dimnames(intensities$bg) <- nam
	slot(out, "pmIntensities", check=FALSE) <- intensities$pm
	slot(out, "bgIntensities", check=FALSE) <- intensities$bg
    return(out)
}

readCelIntensities2.offline <- function (filenames1, filenames2,  pkgname) {   
	pkg <- get(pkgname)
	pmIndex <- pmindex(pkg)
	bgIndex <- bgindex(pkg)
		
    filenames1 <- file.path(dirname(filenames1), basename(filenames1))
    filenames2 <- file.path(dirname(filenames2), basename(filenames2))
    all.headers <- lapply(as.list(c(filenames1, filenames2)), readCelHeader)
    chiptype <- unique(sapply(all.headers, function(x) x$chiptype))
    if (length(chiptype) != 1) {
        warning("The CEL files do not have the same chiptype.")
    }
    nfiles <- length(filenames1)

	pmIntensities <- ff(dim=c(length(pmIndex), nfiles, 2),
	 	vmode="double", finalizer="close")
	bgIntensities <- ff(dim=c(length(bgIndex), nfiles, 2),
	 	vmode="double", finalizer="close")

    colnames(pmIntensities) <- basename(filenames1)
    colnames(bgIntensities) <- basename(filenames1)
    for (i in 1:nfiles) {
        tmp <- readCel(filename = filenames1[i], indices=NULL,
            readIntensities = TRUE, readHeader = FALSE, readStdvs = FALSE,
            readPixels = FALSE, readXY = FALSE, readOutliers = FALSE,
            readMasked = FALSE)
        pmIntensities[, i, 1] <- tmp$intensities[pmIndex]
        bgIntensities[, i, 1] <- tmp$intensities[bgIndex]
        tmp <- readCel(filename = filenames2[i], indices=NULL,
            readIntensities = TRUE, readHeader = FALSE, readStdvs = FALSE,
            readPixels = FALSE, readXY = FALSE, readOutliers = FALSE,
            readMasked = FALSE)
        pmIntensities[, i, 2] <- tmp$intensities[pmIndex]
        bgIntensities[, i, 2] <- tmp$intensities[bgIndex]
    }
	list(pm=pmIntensities, bg=bgIntensities)
}


if (FALSE) {
	rowMedians.ff <- function(x, BATCHSIZE=NULL, na.rm=TRUE, returnFF=FALSE) {
		if(is.null(BATCHSIZE)) BATCHSIZE=5000
		ret <- ffrowapply(rowMedians(x[i1:i2,,drop=FALSE], na.rm=na.rm), X=x, 
			RETURN=TRUE, RETCOL=NULL, BATCHSIZE=BATCHSIZE)	
		if (returnFF) {
			return(ret)
		} else {
			return(ret[])
		}
	}

	rowMedians.ff_array <- function(imat, na.rm=TRUE) { # Operate on 1st level of dim 3
		ret <- ffrowapply(rowMedians(imat[i1:i2,,1], na.rm=na.rm), X=imat, 
			RETURN=TRUE, RETCOL=NULL, BATCHSIZE=5000)	
		return(ret[])
	}


	setMethod("rowMedians", "ff_array", rowMedians.ff_array)
}

