methp <- function(dat, spatial=TRUE, bgSubtract=TRUE,
		withinSampleNorm="loess",
		scale=c(0.99, 0.99),
		betweenSampleNorm="quantile", 
		controlProbes=c("CONTROL_PROBES", "CONTROL_REGIONS"),
		controlIndex=NULL, 
		commonMethPercentParams=NULL,
		verbose=TRUE, returnM=FALSE, 
		plotDensity=NULL, plotDensityGroups=NULL) {
		
	if(!is.null(plotDensity)) {
		pdf(file=plotDensity, height=11, width=8)
		par(mfrow=c(5,2), mar=c(2,2,4,2))
		lwd <- rep(1, ncol(dat))
		plotDensity(dat, main="1. Raw", lab=plotDensityGroups, 
		 	controlIndex=controlIndex)
	}
	
	if (is.list(betweenSampleNorm)) {
		bs <- betweenSampleNorm
	} else {
		if (betweenSampleNorm=="quantile") {
			bs <- list(m="allQuantiles", untreated="none", enriched="none")
			if(is.null(commonMethPercentParams)) commonMethPercentParams <- TRUE
		} else if (betweenSampleNorm=="sqn") {
			bs <- list(m="none", untreated="complete", enriched="sqn")
			if(is.null(commonMethPercentParams)) commonMethPercentParams <- FALSE
		} else if (betweenSampleNorm=="sqn99") {
			stop("Option sqn99 is deprecated. The same behaviour is now obtained by using sqn.")
		} else if (betweenSampleNorm=="none") {
			bs <- list(m="none", untreated="none", enriched="none")
			if(is.null(commonMethPercentParams)) commonMethPercentParams <- FALSE
		}
	}	
    # Spatial bias correction
    if (spatial) {
        if (verbose) message("Spatial normalization")
       	dat <- spatialAdjust(dat)
    }
    # Background removal
	if (bgSubtract) {
    	if (verbose) message("Background removal")
		dat <- bgAdjustBgp(dat)
	}
	if(!is.null(plotDensity)) {
		plotDensity(dat, main="2. After spatial & bg", lab=plotDensityGroups, controlIndex=controlIndex)
	}	
	# Within sample normalization
	if (is.null(controlIndex)) {
		controlIndex <- getControlIndex(dat, controlProbes=controlProbes)
	}
	if (verbose) message("Within sample normalization: ", withinSampleNorm) 
	dat <- normalizeWithinSamples(dat, method=withinSampleNorm,
		scale=scale, controlIndex=controlIndex, verbose=verbose)
	if(!is.null(plotDensity)) {
		plotDensity(dat, main="3. After within-sample norm", 
			lab=plotDensityGroups, controlIndex=controlIndex)
	}
    # Between sample normalization    
	if (verbose) {
		message("Between sample normalization", appendLF=FALSE)
		if (is.list(betweenSampleNorm)) {
			message(". M: ", bs$m, ", Untreated channel: ", bs$untreated, 
		        ", Methyl-depleted channel: ", bs$enriched)
		} else {
			message (": ", betweenSampleNorm)
		}
	}	
    dat <- normalizeBetweenSamples(dat, m=bs$m, 
		untreated=bs$untreated, 
        enriched=bs$enriched, controlProbes=controlProbes, 
		controlIndex=controlIndex, verbose=verbose)
	M <- getM(dat)[pmindex(dat),,drop=FALSE]

	if(!is.null(plotDensity)) {
		plotDensity(dat, main="4. After between-sample norm",
		 lab=plotDensityGroups, controlIndex=controlIndex)
	}		
	
	if (returnM=="TRUE" | returnM=="+") {
		retval <- M
	} else if (returnM=="-") {
		retval <- -M
	} else {
	    if (verbose) message("Estimating percentage methylation")
     	retval <- methPercent(m=M, commonParams=commonMethPercentParams,
	 		ngc=countGC(dat))
	}
	if(!is.null(plotDensity)) {
		if (is.null(controlIndex)) controlIndex <- getControlIndex(dat)
		if(returnM=="FALSE") plotDensity(retval, main="5. Percentage methylation", rx=c(0,1), lab=plotDensityGroups, controlIndex=controlIndex)
		dev.off()		
	}
    return(retval)
}


scaleSamples <- function(dat, scale=c(0.99, 0.99)) {
	if (length(scale)==2) {
		pms <- pm(dat)
		c1 <- log2(pms[,,1,drop=FALSE])
		c2 <- log2(pms[,,2,drop=FALSE])
		M <- c1-c2		
		x <- apply(M, 2, quantile, scale[1], na.rm=TRUE)
		adj <- x / -log(1-scale[2])
		M <- sweep(M, 2, adj, FUN="/")
		c2 <- c1-M
		pms[,,2] <- 2^c2
		pm(dat) <- pms
	}
	return(dat)
}		


readCharm <- function(files, path=".", ut="_532.xys", md="_635.xys", 
		sampleKey, sampleNames=NULL, pkgname,
		type=NULL, ...) {
    files <- as.character(files)
    o <- order(files)
    files <- files[o]
    if (!is.null(sampleNames)) sampleNames <- as.character(sampleNames[o])
    utIdx <- grep(ut, files)
	if (length(utIdx)==0) {
		stop("No files match the untreated extension ", ut, "\nPlease use the ut option to set the correct extension\n")
	}
    mdIdx <- grep(md, files)
	if (length(mdIdx)==0) {
		stop("No files match the methyl-depleted extension ", md, "\nPlease use the md option to set the correct extension\n")
	}
    filesUt <- files[utIdx]
    filesMd <- files[mdIdx]
    if (!all(sub(ut, "", filesUt) == sub(md, "", filesMd))) 
        stop(("The untreated (ut) and methyl-depleted (md) file names don't match up\n"))
    if (!is.null(sampleNames)) {
        sampleCheck <- sampleNames[utIdx] == sampleNames[mdIdx]
        if (!all(sampleCheck)) 
            stop("The untreated (ut) and methyl-depleted (md) sample names don't match up\n Check:", 
                sampleNames[utIdx][!sampleCheck])
        sampleNames <- sampleNames[utIdx]
    } else {
        sampleNames <- sub(ut, "", filesUt)
    }

	if (!missing(sampleKey)){
		sampleKey <- sampleKey[o,]
		keep <- apply(sampleKey[utIdx,]==sampleKey[mdIdx,], 2, all)
		extraCols <- sampleKey[utIdx, keep]
	} else {
		extraCols <- matrix(nrow=length(utIdx), ncol=0)
	}
	if (!is.null(type)) {
		type <- as.character(type[o])
		pd <- data.frame(extraCols, type=type[utIdx],
			 arrayUT=filesUt, arrayMD=filesMd, stringsAsFactors=FALSE)
	} else {
		pd <- data.frame(extraCols, arrayUT=filesUt, arrayMD=filesMd, stringsAsFactors=FALSE)
	}
	vpd <- data.frame(
		labelDescription=c(rep(NA, ncol(pd)-2), 
			"Untreated channel file name", 
			"Methyl-depleted channel file name"),
		channel=factor(c(rep("_ALL_", ncol(pd)-2), 
			"channel1", "channel2")))
			
	## TEMPORARY FIX TO REMAIN BACKWARD COMPATABILITY WITH ##
	## oligo 1.10.x (Bioconductor 2.5). This fix avoids    ##
	## a warning 										   ## 		
	version <- as.numeric(sub(".*\\.(.*)\\..*", "\\1",
	 	packageDescription("oligo", field="Version")))
	if (version<11) {	
		vpd <- data.frame(
			labelDescription=c(rep(NA, ncol(pd)-2), 
				"Untreated channel file name", 
				"Methyl-depleted channel file name"))
	}
	## END TEMPORARY FIX ##
	 		
    pdd <- new("AnnotatedDataFrame", data=pd, varMetadata=vpd)
    sampleNames(pdd) <- sampleNames  
   	dat <- read.xysfiles2(channel1=file.path(path, filesUt), 
			channel2=file.path(path, filesMd), pkgname=pkgname,
        	phenoData=pdd, sampleNames=sampleNames, ...)     
	return(dat)
}


plotDensity <- function(dat, rx=c(-4,6), controlIndex=NULL, 
		pdfFile=NULL, main=NULL, lab=NULL) {
	if (!is.null(pdfFile)) {
		pdf(pdfFile)
		par(mfcol=c(2,1))
	}
	if (any(class(dat)=="TilingFeatureSet2") |
	 		any(class(dat)=="TilingFeatureSet")) {
		M <- getM(dat)[pmindex(dat),,drop=FALSE]
		if (is.null(lab)) lab <- sampleNames(dat)
	} else {
		M <- dat
		if (is.null(lab)) lab <- colnames(dat)		
	}
	if (is.null(controlIndex)) controlIndex <- getControlIndex(dat)
	plotDensityMat(M, xlab="M", lab=lab, 
		main=paste(main,"\nAll probes"), rx=rx)
	plotDensityMat(M[controlIndex,,drop=FALSE], xlab="M", lab=lab, 
		main=paste(main, "\nControl probes"), rx=rx)
	if (!is.null(pdfFile)) dev.off()
}

normalizeLoess <- function(dat, controlIndex=NULL, 
		controlProbes=c("CONTROL_PROBES", "CONTROL_REGIONS"), span=0.3,
		by="A", approx=TRUE, breaks=1000) {
	if (is.null(controlIndex)) {
    	controlIndex <- getControlIndex(dat, controlProbes=controlProbes)
	}
	pms <- pm(dat)
	for (i in 1:dim(dat)["Samples"]) {
		c1 <- log2(pms[,i,"channel1"])
		c2 <- log2(pms[,i,"channel2"])
		y <- c1-c2
		if (by=="tot") {
			x <- c1
		} else if (by=="A") {
			x <- (c1+c2)/2
		}
 		fit <- loess(y ~ x, subset = controlIndex,
               na.action = na.exclude, degree=1, surface = "direct",
 			   span=span)
		adj <- predictLoess(fit, newdata=x, approx=approx, breaks=breaks)
		pms[, i, "channel2"] <- 2^(c2 + adj)
	}
	pm(dat) <- pms
	return(dat)
}


predictLoess <- function(fit, newdata, approx=TRUE, breaks=1000) {
	if (approx) {
		newdataBin <- cut(newdata, breaks=breaks, ordered.result=TRUE)
		binMean <- tapply(newdata, newdataBin, mean)
		isNa <- which(is.na(binMean))
		adj <- rep(NA, length(binMean))
		if (length(isNa)==0) {
			adj <- predict(fit, newdata = binMean)
		} else {
			adj[-isNa] <- predict(fit, newdata = binMean[-isNa])
		}
		ret <- adj[newdataBin]
	} else {
		ret <- predict(fit, newdata=newdata)
	}
	return(ret)
}
	

regionFilter <- function(x, region, f) {
	x<-unlist(tapply(x, region, myfilter, f))
	names(x) <- NULL
	return(x)
}

normalizeWithinSamples <- function (dat, method = "loess", 
	scale=c(0.99, 0.99), 
	controlProbes = c("CONTROL_PROBES", "CONTROL_REGIONS"), 
	controlIndex = NULL, approx=TRUE, breaks=1000, 
	verbose=FALSE) {
    #if (method=="gc") {
	if (grepl("gc", method)) {
        mAdj <- diffAmpEstGC(dat, method = method, controlProbes = controlProbes)
        dat <- diffAmpAdjustGC(dat, mAdj)
    }
	if (grepl("loess", method)) {
			dat <- normalizeLoess(dat, controlIndex=controlIndex, 
				controlProbes=controlProbes, 
				approx=approx, breaks=breaks)
	}
	if (grepl("median", method)) {
        if (is.null(controlIndex)) {
            controlIndex <- getControlIndex(dat, controlProbes = controlProbes)
		}
        datPm <- log2(pm(dat))
        M <- datPm[, , "channel1"] - datPm[, , "channel2"]
        mCtrl <- M[controlIndex, ]
        mAdj <- apply(mCtrl, 2, median, na.rm=TRUE)
        datPm[, , "channel2"] <- sweep(datPm[, , "channel2"], 
            2, mAdj, FUN = "+")
        pm(dat) <- 2^datPm
    }
	dat <- scaleSamples(dat, scale) 
    return(dat)
}

normalizeBetweenSamples <- function(dat, m="allQuantiles",
 		untreated="none", enriched="none", 
		controlProbes=c("CONTROL_PROBES", "CONTROL_REGIONS"), 
		controlIndex=NULL, verbose=FALSE) {    

	pms <- pm(dat)
	M <- getM(dat)[pmindex(dat),,drop=FALSE]	
	if (m!="none"){
		M <- normQuantile(M, method=m)
		for (i in 1:ncol(M)) {
			pms[,i,"channel2"] <- 2^(log2(pms[,i,"channel1"]) - M[,i])	
		}		
	}
    ## Normalize untreated (leaving M unchanged)
    if (untreated=="complete") {
       	m <- rowMedians(log2(pms[,,"channel1"]))
		for (i in 1:ncol(pms)) {
			pms[,i,"channel1"] <- 2^m
		}
    } else {
        pms[,,"channel1"] <- 2^(normQuantile(log2(pms[,,"channel1"]), method=untreated))
    }
    for (i in 1:ncol(M)) {
		pms[,i,"channel2"] <- 2^(log2(pms[,i,"channel1"]) - M[,i])
	}
    pm(dat) <- pms
    ## Normalize enriched 
	if (enriched=="sqn") {
		if(is.null(controlIndex)) controlIndex <- getControlIndex(dat, controlProbes=controlProbes) 
		pms[,,"channel2"] <- 2^(SQN(y=log2(pms[,,"channel2"]), ctrl.id=controlIndex))
	} else {
        pms[,,"channel2"] <- 2^(normQuantile(log2(pms[,,"channel2"]), method=enriched))
    }
    pm(dat) <- pms
    return(dat)
} 


normQuantile <- function(x, method="allQuantiles") {
    if (method=="99thQuantile") method="99" # Legacy option
    if (method=="none") {
            # Pass through option
    } else if (method=="allQuantiles") {
		if (any(class(x)=="ff")) {
			x[,] <- normalize.quantiles(x[,], copy=FALSE) 			
		} else {
        	x <- normalize.quantiles(x,copy=FALSE) 
		}
    } else if (!is.na(as.numeric(method)) & (as.numeric(method)>=0) & (as.numeric(method)<=100)) {
		if (any(class(x)=="ff")) {
			qnt <- rep(NA, ncol(x))
			for (i in 1:ncol(x)) {
				qnt[i] <- quantile(x[,i], as.numeric(method)/100, na.rm=TRUE)
			}
			adj <- qnt/median(qnt)
			for (i in 1:ncol(x)) {
				x[,i] <- x[,i] / adj[i]
			}
		} else {
	        qnt <- apply(x, 2, quantile, as.numeric(method)/100, na.rm=TRUE)
	        adj <- qnt / median(qnt)
	        x <- sweep(x, 2, adj, FUN="/")
		}
    } else {
        stop("Invalid option '", method, "' to normQuantile\n")
    }
    return(x)
} 

## Uses all control and bg probes passed to it
normControlBg <- function(pms, bgs=NULL, controlIndex=NULL, affinity=NULL) {
    mBg <- NULL
    ctrl <- NULL
    mCtrl <- NULL
    if (!is.null(bgs)) {
        mBg <- rowMedians(bgs)
    }
    if (!is.null(controlIndex)) {
        ctrl <- pms[controlIndex,]
        mCtrl <- rowMedians(ctrl)
    }
    message(" Normalizing with", length(c(mBg, mCtrl)), "control probes ", appendLF=FALSE)
    if (!is.null(affinity)) {
        message("(affinity-based mapping)")
    } else {
        message("(intensity-based mapping)")
    }
  
    for (i in 1:ncol(pms)) {
        if (is.null(affinity)) {
            map <- getMap(x=c(bgs[,i], ctrl[,i]), m=c(mBg, mCtrl))
            pms[,i] <- pms[,i]-map(pms[,i])   
        } else {
            map <- getMap(x=ctrl[,i], m=mCtrl, affinity=affinity[controlIndex])
            pms[,i] <- pms[,i]-map(affinity)   
        }
    }
    return(pms)
}

getMap <- function (x, m, affinity=NULL, df = 5) 
{
    if (is.null(affinity)) affinity <- x
    nas <- is.na(x) | is.na(m)
    x <- x[!nas]
    m <- m[!nas]
    affinity <- affinity[!nas]
    d <- x - m
    lmFit = lm(d ~ ns(affinity, df = df))
    return(function(x1) predict(lmFit, data.frame(affinity = x1)))
}




getControlIndex <- function(dat, 
		controlProbes=c("CONTROL_PROBES", "CONTROL_REGIONS"), 
		noCpGWindow=1000, subject, onlyGood=FALSE, matrix=TRUE) {
    if (missing(subject)) {
		file <- file.path(system.file(package = annotation(dat)), 
		    "data", "controlIndex.rda")
		if (file.exists(file)) {
		    load(file)
		} else {
        	controlIndex <- which(getContainer(dat) %in% controlProbes)
		}
    } else {
        if (class(subject)!="BSgenome")
            stop("You must supply a BSgenome object for subject if using the noCpGWindow option. E.g. Hsapiens from the BSgenome.Hsapiens.UCSC.hg18 package.\n")
		pv <- providerVersion(subject)	
		file <- file.path(system.file(package = annotation(dat)), 
			    "data",	paste(pv, "_noCpG", noCpGWindow, ".rda", sep=""))
		if (file.exists(file)) {
		    load(file)
		} else {
	    	chr <- pmChr(dat)
	        pos <- pmPosition(dat)
	        cpgd <- cpgdensity(subject, chr=chr, pos=pos, 
				windowSize=noCpGWindow)
	        controlIndex <- which(cpgd==0)
		}
    }
	return(controlIndex)
}



diffAmpEstGC <- function(dat, method="median", controlProbes="CONTROL_REGIONS", smooth=2) {
    #kern <- function(u, type="epanechnikov", smooth=1) {
    #    ifelse(abs(u)<smooth, (1/smooth) * (3/4) * (1-(u/smooth)^2), 0)
    #}
    Ngc <- countGC(dat, "pm")
    controlIndex <- which(getContainer(dat)%in%controlProbes)
    ngc <- Ngc[controlIndex]
    datPm <- log2(pm(dat)) 
    M <- datPm[,,"channel1"] - datPm[,,"channel2"] 
    adj <- sapply (1:ncol(M), function(samp) {
        m <- M[controlIndex, samp]
        a <- (datPm[,,"channel1"][controlIndex,samp] + 
              datPm[,,"channel2"][controlIndex,samp])/2
        nas <- is.na(m) | is.na(a)
        m <- m[!nas]       
        a <- a[!nas]
        ngc <- ngc[!nas]
        breakpoint <- optimize(checkBreakpoint, range(a), ms=m, as=a)$minimum
        keep <- a>breakpoint
        if (sum(keep)/length(controlIndex) < (1/10)) { # Make sure we use at least 1/10 the control probes
            o <- order(a, decreasing=TRUE)
            keep[o[1:(length(controlIndex)/10)]] <- TRUE
        }        
        med <- tapply(m[keep], ngc[keep], method)
        weights <- tapply(m[keep], ngc[keep], length)
        ngc_adj <- est <- w <- rep(NA, max(Ngc))
        est[as.numeric(names(med))] <- med
        w[as.numeric(names(weights))] <- weights
        for (i in which(!is.na(est))) {
            k <- dnorm(abs(i-1:length(ngc_adj))/smooth)       
            #k <- kern(abs(i-1:length(ngc_adj)), smooth=smooth)
            ngc_adj[i] <- sum(k*w*est, na.rm=TRUE) / sum(k*w, na.rm=TRUE)
        }
        for (i in 2:length(ngc_adj)) {
            if (is.na(ngc_adj[i]) & !is.na(ngc_adj[i-1]))
                ngc_adj[i] <- ngc_adj[i-1]
        }
        for (i in (length(ngc_adj)-1):1) {
            if (is.na(ngc_adj[i]) & !is.na(ngc_adj[i+1]))
                ngc_adj[i] <- ngc_adj[i+1]
        }
        ngc_adj
    })   
    colnames(adj) <- sampleNames(dat)
    rownames(adj) <- paste("GC=", 1:nrow(adj), sep="")
    adj
}


diffAmpAdjustGC <- function (dat, adj){
    datPm <- log2(pm(dat)) 
    Ngc <- countGC(dat)
    for (samp in 1:dim(dat)["Samples"]) {
        datPm[,samp,"channel2"] <- datPm[,samp,"channel2"] + adj[Ngc, samp]
    }
    pm(dat) <- 2^datPm
    return(dat)
}


    
countGC <- function(dat, type="pm") {   
    if (type=="pm") {
        seqs=as.character(pmSequence(dat))
    } else if (type=="bg") {
        seqs=as.character(bgSequence(dat))            
    } else {
        stop("Invalid type. Choose 'pm' or 'bg'\n")
    }
    bc <- oligo:::basecontent(seqs)
    bc[,"C"] + bc[,"G"]
}

cpgdensity <-function(subject, chr, pos, windowSize=500, sequence="CG") {
    idx <- split(1:length(chr), chr)
    s <- DNAString(sequence) 
    cpgdensity <- rep(NA, length(pos))
    for (curchr in (names(idx))) {
            if (curchr %in% names(subject)) {
                chrseq <- subject[[curchr]]
                curpos <- pos[idx[[curchr]]]
                v <- suppressWarnings(Views(chrseq, start=curpos-windowSize/2, end=curpos+windowSize/2))
                d <- suppressWarnings(DNAStringSet(v))
                numcpg <- vcountPattern(s, d, fixed=TRUE)
                cpgdensity[idx[[curchr]]] <- numcpg/windowSize       
            }
    }
    cpgdensity
}

demedian <- function(dat) {
    X <- getX(dat, "pm")
    Y <- getY(dat, "pm")
    pms <- pm(dat)
    for (chan in c("channel1", "channel2")) {
	    for (samp in sampleNames(dat)) {
	        pm <- log2(pms[,samp,chan])
	        m <- median(pm)
	        tmp <- matrix(nrow=max(Y), ncol=max(X))
	        for (i in 1:length(pm)) {
	            tmp[Y[i], X[i]] <- pm[i]
	        }
	        rm <- rowMedians(tmp, na.rm=TRUE)
	        tmp <- sweep(tmp, 1, rm-m)
	        tmp <- t(tmp)
	        cm <- rowMedians(tmp, na.rm=TRUE)
	        tmp <- sweep(tmp, 1, cm-m)
	        pms[,samp,chan] <- 2^(sapply(1:length(pm), function(i) tmp[X[i], Y[i]]))
	    }
    }
    pm(dat) <- pms
    return(dat)
}

qcReport <- function(dat, file=NULL, utRange=c(30,100), enRange=c(8,12), numProbes=5e+5, blockSize) {
    # Calculate summary quality scores
	n <- nrow(pm(dat))
	if (numProbes==-1) numProbes <- n
	if (n < numProbes) numProbes <- n
	idx <- as.integer(seq(1, n, length.out=numProbes))
    pmQual <- pmQuality(dat, idx=idx) 
    X <- getX(dat, "pm")[idx]
    Y <- getY(dat, "pm")[idx]
	if(missing(blockSize)) blockSize <- getBlockSize(dat)
    imgs1 <- arrayImage(X,Y, pmQual, blockSize=blockSize)
 	#sd1 <- unlist(lapply(imgs1, function(x) sd(as.vector(x$z), na.rm=TRUE)))

	tmp <- arrayImage(X,Y, log2(pm(dat)[idx,,"channel1"]),
		blockSize=blockSize)
	sd1 <- unlist(lapply(tmp, function(x) sd(as.vector(x$z), na.rm=TRUE)))
	imgs2 <- arrayImage(X,Y, log2(pm(dat)[idx,,"channel2"]),
		blockSize=blockSize)
 	sd2 <- unlist(lapply(imgs2, function(x) sd(as.vector(x$z), na.rm=TRUE)))

	sdRange <- c(0, max(sd1, sd2)*1.1)
	if (any(class(pmQual)=="ff")) {
		pmSignal <- sapply(1:ncol(pmQual), function(i) mean(pmQual[,i]))
		names(pmSignal) <- colnames(pmQual)
	} else {
		pmSignal <- colMeans(pmQual, na.rm=TRUE) 
	}
	
    if (!is.null(file)) {
        pdf(file=file, width=8, height=10.5)    
        lay <- rbind(cbind(rep(1,3), matrix(rep(2:4, each=6), ncol=6)), rep(5,7))
		layout(lay)
		#layout(matrix(c(1,1,1,4, 2,2,2,4, 3,3,3,4 ), ncol=3))

	    n <- length(pmSignal)
		if (is.null(names(pmSignal))) {
			names(pmSignal) <- paste("Sample", 1:n)
		}
		l <- length(pmSignal)

		longestName <- max(sapply(names(pmSignal), nchar))
	    cexh <- ifelse (n>50, max(0.25, 1-(n-50)*.01), 1)
	    cexw <- ifelse (longestName>17, max(0.25, 1-(longestName-15)*0.05), 1)
		cex <- min(cexh, cexw)
		par(oma=c(2,2,4,2))

		o <- order(names(pmSignal))
		#o <- order(names(pmSignal), decreasing=TRUE)
		par(mar=c(5,0,3,0))
		plot.new()
		plot.window(xlim=c(0,10), ylim=c(0,l))
		text(1,l:1-0.5, names(pmSignal)[o], adj=0.1, cex=cex)
	    abline(h=1:(l+1)-1, lwd=0.5, lty=3)

	    xmin <- max(0, floor((min(pmSignal)-5)/5)*5)
	    xmax <- min(100, ceiling((max(pmSignal)+5)/5)*5)
	    par(mar=c(5,1,3,1))
	    plot(pmSignal[o], l:1, las=1, ylab="", xlab="", type="p", 
	        pch=19, xlim=c(xmin, xmax), ylim=c(0.5, l+0.5),
	        yaxt="n", frame.plot=TRUE, main="Signal strength")
	    abline(h=1:(l+1)-0.5, lwd=0.5, lty=3)
	    #axis(2, at=1:l, labels=names(pmSignal), las=1, tick=FALSE, cex.axis=cex)

	    #xmin <- min(sd1)-1
	    #xmax <- max(sd1)+1
	    par(mar=c(5,1,3,1))
	    plot(sd1[o], l:1, las=1, ylab="", xlab="", type="p", 
	        pch=19, xlim=sdRange, ylim=c(0.5, l+0.5),
	        yaxt="n", frame.plot=TRUE, main="Channel 1\nstandard deviation")
	    abline(h=1:(l+1)-0.5, lwd=0.5, lty=3)

	    #xmin <- min(sd2)-0.01
	    #xmax <- max(sd2)+0.01
	    par(mar=c(5,1,3,1))
	    plot(sd2[o], l:1, las=1, ylab="", xlab="", type="p", 
	        pch=19, xlim=sdRange, ylim=c(0.5, l+0.5),
	        yaxt="n", frame.plot=TRUE, main="Channel 2\nstandard deviation")
	    abline(h=1:(l+1)-0.5, lwd=0.5, lty=3)

	    #axis(1, tick=FALSE)
	    hist(pmSignal, main="Signal strength histogram", xlab="Signal strength")

        # Plot UT channel probe quality (PM vs BG)
        l <- layout(matrix(c(1,1,1,1,2:21), 6, 4, byrow = TRUE)) # 
        plot.new()
        text(0.5, 0.2, "Untreated Channel: PM probe quality", cex=2)
        par(mar=c(4.5,8,3.5,3))
        arrayPlot(imgs1[o], r=utRange)
        # Plot MD channel
        l <- layout(matrix(c(1,1,1,1,2:21), 6, 4, byrow = TRUE)) # 
        plot.new()
        text(0.5, 0.2, "Enriched Channel: PM signal intensity", cex=2)
        par(mar=c(4.5,8,3.5,3))
        arrayPlot(imgs2[o], r=enRange)
        dev.off()
    }
    return(as.data.frame(cbind(pmSignal, sd1, sd2)))
}

arrayImage <- function(x,y,z, view="2d", blockSize=50) {
    if (is.null(dim(z))) z <- as.matrix(z)
    if (view=="col") {
		tmp <- vec2array(x,y,z)
        ret <- apply(tmp, 3, rowMeans, na.rm=TRUE)
        colnames(ret) <- colnames(z)
    }
    if (view=="row") {
		tmp <- vec2array(x,y,z)
        ret <- apply(aperm(tmp, c(2,1,3)), 3, rowMeans, na.rm=TRUE)  
        colnames(ret) <- colnames(z)
    }
    if (view=="2d") {
		nx <- max(x)/blockSize
		ny <- max(y)/blockSize
		d <- discretize.image(cbind(y,x), m=ny, n=nx)
		ret <- apply(z, 2, function(vec) {
			as.image(vec, ind=d$index, nx=d$m, ny=d$n, na.rm=TRUE)		
        })		
        names(ret) <- colnames(z)
    }
	return(ret)
}

arrayPlot <- function(imgs, xlab="NULL", r=NULL) {
    if (is.list(imgs)) {
        if(is.null(r)) r <- range(unlist(lapply(imgs, function(x) x$z)), na.rm=TRUE)
        for (i in 1:length(imgs)) {
			curImg <- imgs[[i]]$z
			curImg <- apply(curImg, 1, pmax, r[1])
			curImg <- apply(curImg, 1, pmin, r[2])
            image.plot(curImg, main=names(imgs)[i],  zlim=r, xaxt="n", yaxt="n", horizontal=TRUE)    
        }
    } else {
		# No longer needed
        #r <- range(imgs$z, na.rm=TRUE)
        #for (i in 1:ncol(imgs$z)) {
        #    plot(imgs[,i], main=colnames(imgs)[i], type="l", ylim=r, ylab="pm quality", xlab=xlab)
        #}
    }
}

plotDensityMat <- function (x, main = NULL, 
    ylab = "Density", xlab = "M", lab = NULL, rx = NULL, ry = NULL, 
    legendPos = "topright", cex=0.8) {
	lab <- as.factor(lab)
    d <- apply(x, 2, density, na.rm = TRUE)
    if (is.null(rx)) rx <- range(sapply(d, function(i) i$x))
    if (is.null(ry)) ry <- range(sapply(d, function(i) i$y))
    plot(rx, ry, type = "n", xlab = xlab, ylab = ylab, main = main)
    sapply(1:length(d), function(i) lines(d[[i]], col = lab[i]))
    legend(legendPos, legend = levels(lab), 
		text.col = 1:length(levels(lab)), cex = cex)
}

## fn(pm) where fn is the ECDF of bg
pmQuality <- function(dat, channel="channel1", verbose=FALSE, idx=NULL) {
	if (is.null(idx)) idx <- 1:nrow(pm(dat))
    Ngc <- countGC(dat, "pm")[idx]
    bgNgc <- countGC(dat, "bg")  

	pms <- pm(dat)[idx,,,drop=FALSE]
	bgs <- bg(dat)
	pmq <- sapply(1:ncol(pms), function(i) {
	    fn <- tapply(bgs[,i,channel], bgNgc, ecdf)
	    ret <- rep(NA, length(Ngc))
		for (ngc in unique(Ngc)) {
	        idx <- Ngc==ngc
	        closestIdx <- order(abs(as.numeric(names(fn))-ngc))[1]
	        bgngc <- as.character(names(fn)[closestIdx])
			ret[idx] <- 100 * fn[[bgngc]](pms[idx,i,channel])   
	    }
		return(ret)
    })
	colnames(pmq) <- sampleNames(dat)
	return(pmq)
}

pmvsbg <- function(...) {
	message("The pmvsbg function has been renamed pmQuality. Both names work for now but please update your code soon")
	pmQuality(...)
}

dynamicRange <- function(dat, prob=0.8) {
    Ngc <- countGC(dat, "pm")
    bgNgc <- countGC(dat, "bg")
    c1 <- log2(oligo::pm(dat)[,,"channel1"])
    c2 <- log2(oligo::bg(dat)[,,"channel2"])    
    sapply(sampleNames(dat), function(i) {
            tmp <- tapply(c2[,i], bgNgc, quantile, prob)
            c2med <- rep(NA, max(as.numeric(names(tmp))))
            c2med[as.numeric(names(tmp))] <- tmp
            c1[,i] - c2med[Ngc]
    })    
}

vec2array <- function(X,Y,Z) {
    if (ncol(Z)==1) Z <- as.matrix(Z)
    maxY <- max(Y)
    maxX <- max(X)
    o <- order(X,Y)
    Y <- Y[o]
	X <- X[o]
    tmp <- tapply(Y, X, function(y) as.numeric(y))
    tmp <- unlist(lapply(names(tmp), function(i) tmp[[i]]+((as.numeric(i)-1)*maxY)))
	Z <- as.matrix(Z[o,])
    d <- matrix(NA, nrow=maxX*maxY, ncol=ncol(Z))
    d[tmp,] <- Z
    array(d, dim=c(maxY, maxX, ncol(Z)))
}

countSeq <- function(subject, chr, start, end, seq) {
    idx <- split(1:length(chr), chr)
    retval <- rep(NA, length(chr))
    for (curChr in (names(idx))) {
        if (curChr %in% names(subject)) {
            chrSeq <- unmasked(subject[[curChr]])
            curStart <- start[idx[[curChr]]]
            curEnd <- end[idx[[curChr]]]
            v <- suppressWarnings(Views(chrSeq, start=curStart, end=curEnd))
            d <- suppressWarnings(DNAStringSet(v))
            s <- DNAString(seq)
            numSeq <- vcountPattern(s, d, fixed=TRUE)
            retval[idx[[curChr]]] <- numSeq
        } else {
            #cat("\ncountSeq: Skipping chr", curChr, "\n")
        }
    }    
    return(retval)
}



spatialAdjustVec <- function(z, d, ims=NULL, theta=1) {
	if (is.null(ims)) {
		im <- as.image(z, ind=d$index, nx=d$m, ny=d$n, na.rm=TRUE)
		ims <- image.smooth(im, theta=theta)
	}
	adj <- ims$z - median(ims$z)
	adjV <- as.vector(t(adj))
	ind <- d$index
	idx <- (ind[,1]-1)*d$n + ind[,2]
	zAdj <- z - adjV[idx]
	return(list(zAdj=zAdj, ims=ims))
}

spatialAdjust <- function(dat, blockSize, theta=1) {
	if(missing(blockSize)) blockSize <- getBlockSize(dat)
    x <- getX(dat, "pm")
    y <- getY(dat, "pm")
    bgX <- getX(dat, "bg")  
    bgY <- getY(dat, "bg")     
    NCOL=dim(dat)["Samples"]
	pms <- pm(dat)
	bgs <- bg(dat)
	nx <- max(x)/blockSize
	ny <- max(y)/blockSize
	d <- discretize.image(cbind(y, x), m=ny, n=nx)
	dBg <- discretize.image(cbind(bgY, bgX), grid=d$grid)
	for (i in 1:ncol(pms)) {
		for (channel in 1:2) {
			tot <- log2(pms[,i, channel])
		    bgTot <- log2(bgs[,i, channel])
			sm <- spatialAdjustVec(tot, d)
			pms[,i, channel] <- 2^sm$zAdj
			smBg <- spatialAdjustVec(bgTot, dBg, sm$ims)
			bgs[,i, channel] <- 2^smBg$zAdj
			}
		}
	pm(dat) <- pms
	bg(dat) <- bgs    		
    return(dat)
}

getBlockSize <- function(dat, probesPerBlock=1250) {
	x<-getX(dat, "pm")
	y<-getY(dat, "pm")
	area <-  max(x) * max(y)
	numProbes <- nrow(pm(dat))
	probeDensity <- numProbes/area
	blockSize <- round(sqrt(probesPerBlock/probeDensity))
	return(blockSize)
}


max.density <- function(x, n.pts = 2^14, min.points=30) { # From affy
        if (length(x) >= min.points) {
            aux <- density(x, kernel = "epanechnikov", n = n.pts, 
                na.rm = TRUE)
            aux$x[order(-aux$y)[1]]
        } else {
            return(NA)
        }
}

# Calculate RMA bg parameters using the GC-stratitifed RMA model 
# with background probes
bgParametersBgp <- function (pm, bgpm, Ngc, bgNgc, n.pts = 2^14) 
{
    pmbg <- max.density(bgpm, n.pts)
    pmbg.gc <- tapply(bgpm, bgNgc, max.density, n.pts) # Get mode for each GC bin
    l<-loess(pmbg.gc~as.numeric(names(pmbg.gc)), control = loess.control(surface = "direct")) 
    x <- 1:max(Ngc)
    mubg.gc <- predict(l,(x))
    bg.data <- bgpm[bgpm < pmbg]
    pmbg <- max.density(bg.data, n.pts)
    bg.data <- bgpm[bgpm < pmbg]
    bg.data <- bg.data - pmbg
    bgsd <- sqrt(sum(bg.data^2)/(length(bg.data) - 1)) * sqrt(2)
    sig.data <- pm[pm > pmbg]
    sig.data <- sig.data - pmbg
    expmean <- max.density(sig.data, n.pts)
    alpha <- 1/expmean
    mubg <- pmbg
    list(alpha = alpha, mu = mubg, mu.gc = mubg.gc, sigma = bgsd)
}


# Adjust background using the RMA model with background probes
bgAdjustBgp <- function (dat) {
	pms <- pm(dat)
	bgs <- bg(dat)
    Ngc <- countGC(dat, "pm")
    bgNgc <- countGC(dat, "bg")
	for (samp in 1:ncol(pms)) {
    	for (chan in 1:2) {
            param <- bgParametersBgp(pms[,samp,chan], bgs[,samp,chan], Ngc, bgNgc, n.pts=2^14)
            b <- param$sigma
            pms[,samp,chan] <- pms[,samp,chan] - param$mu.gc[Ngc] - param$alpha * b^2
            pms[,samp,chan] <- pms[,samp,chan] + b * ((1/sqrt(2 * pi)) * exp((-1/2) * ((pms[,samp,chan]/b)^2)))/pnorm(pms[,samp,chan]/b)
            bgs[,samp,chan] <- bgs[,samp,chan] - param$mu.gc[bgNgc] - param$alpha * b^2
            bgs[,samp,chan] <- bgs[,samp,chan] + b * ((1/sqrt(2 * pi)) * exp((-1/2) * ((bgs[,samp,chan]/b)^2)))/pnorm(bgs[,samp,chan]/b)
		}
	}
	isInf <- which(is.infinite(pms))
	if (length(isInf)>0) {
		biggest <- max(pms[-isInf], na.rm=TRUE)		
		pms[isInf] <- biggest
	}
	isInf <- which(is.infinite(bgs))
	if (length(isInf)>0) {
		biggest <- max(bgs[-isInf], na.rm=TRUE)		
		bgs[isInf] <- biggest
	}
   	pm(dat) <- pms
    bg(dat) <- bgs
    return(dat)
}



##Assume the signal is S=X+Y where X~Normal(0, s2) and Y~Exp(alpha)
## Calculate E(Y|S)
logmethParameters <- function (pm, ngc, n.pts = 2^14) 
{
    #max.density <- function(x, n.pts) {
    #    aux <- density(x, kernel = "epanechnikov", n = n.pts, 
    #        na.rm = TRUE)
    #    aux$x[order(-aux$y)[1]]
    #}
    pmbg <- 0
    bg.data <- pm[pm < pmbg]
    #bgsd <- mad(bg.data, center=0, constant=1)
    #bgsd <- mad(bg.data, center=0)
    idx <- pm < pmbg
    bgsd <- median(tapply(pm[idx], ngc[idx], mad, center=0, na.rm=TRUE))
    sig.data <- pm[pm > pmbg]
    sig.data <- sig.data - pmbg
    #cat("Debug logmethParameters(): sig.data =", 
    #    100*round(length(sig.data)/length(pm), 2), "%\n")
    #expmean <- max.density(sig.data, n.pts)
	expmean <- mean(sig.data, na.rm=TRUE)
    alpha <- 1/expmean
    mubg <- pmbg
    list(alpha = alpha, mu = mubg, sigma = bgsd)
}


methPercent <- function(m, ngc, commonParams=TRUE) {
	m <- as.matrix(m)
	param <- t(sapply(1:ncol(m), 
		function(i) logmethParameters(m[,i], ngc)))
	alpha <- unlist(param[,"alpha"])
	sigma <- unlist(param[,"sigma"])
	if(commonParams) {
		alpha[] <- median(alpha)
		sigma[] <- median(sigma)
	}

	ret <- sapply(1:ncol(m), function(i) {
		x <- m[,i]
		a <- alpha[i]
		b <- sigma[i]
		mu <- 0
        f0 <- dnorm(x, mean=0, sd=b)
        f1 <- a * exp(0.5*a^2*b^2-a*x) * (pnorm((x-a*b^2)/b))
        p0.prior <- sum(x<0, na.rm=TRUE) / sum(!is.na(x))
        p0 <- (p0.prior*f0) / ((p0.prior*f0) + ((1-p0.prior)*f1))
        p1 <- 1-p0
        x <- x - mu - a * b^2
        postM <- p1 * (x + b * ((1/sqrt(2 * pi)) * exp((-1/2) * ((x/b)^2)))/pnorm(x/b))
		1-2^(-postM)
   	})
	colnames(ret) <- colnames(m)
	return(ret)
}




checkBreakpoint <- function(breakpoint, ms, as) {
    #plot(a, data, type="p", pch=".")
    tmp <- sapply(breakpoint, function(b) {
        #browser()
        left <- as<b
        right <- as>=b
        l <- length(as)
        #print(sum(left))
        if (sum(left)==l | sum(right)==l) {
            return(1000)
        }else {   
            #print(as[left])
            lmLeft <- lm(ms[left] ~ as[left])
            lmRight <- lm(ms[right] ~ as[right])
            #abline(reg=lmLeft, col="red")
            #abline(reg=lmRight, col="red")
            res <- c(lmLeft$residuals, lmRight$residuals)
            #return(mad(res))
            return(sum(res^2))
        }
    })
    names(tmp) <- breakpoint
    tmp
}

reshapePd <- function(data, sample, channel, fileName, direction="wide", ...) {
    reshape(data, idvar=sample, timevar=channel, v.names=fileName, direction=direction, ...)   
}



myfilter <- function(x,filter,...){
  L=length(x)
  if(L>length(filter)){
    res=filter(x,filter,...)
    ##now fill out the NAs
    M=length(filter)
    N=(M- 1)/2
    for(i in 1:N){
      w=filter[(N-i+2):M]
      y=x[1:(M-N+i-1)]
      res[i] = sum(w*y)/sum(w)
      
      w=rev(w)
      ii=(L-(i-1))
      y=x[(ii-N):L]
      res[ii]<-sum(w*y)/sum(w)
    }
  } else{ res=rep(mean(x),L)}
  return(res)
}

myfilterse <- function(x,filter,...){
  ##filter must add up to 1!!!
  ##x are the SD of the variable being smoothed
  L=length(x)
  if(L>length(filter)){
    res=filter(x,filter^2,...)
    ##now fill out the NAs
    M=length(filter)
    N=(M- 1)/2
    for(i in 1:N){
      w=filter[(N-i+2):M]
      y=x[1:(M-N+i-1)]
      res[i] = sum(w^2*y)/sum(w)^2
      
      w=rev(w)
      ii=(L-(i-1))
      y=x[(ii-N):L]
      res[ii]<-sum(w^2*y)/sum(w)^2
    }
  } else { res=rep(mean(x)/L,L)}
  return(sqrt(res))
}

rowMads <- function(x,center=NULL,constant=1.4826){
  if(is.null(center)) center=rowMedians(x)
  constant*rowMedians(abs(x-center))
}

smoothTile <- function(mat,pns,filter=NULL,ws=3,verbose=TRUE,...){
  ##... goes to filter
  if(is.null(dim(mat))) mat=as.matrix(mat)
  if(is.null(filter)){
    Tukey = function(x) pmax(1 - x^2,0)^2
    filter= Tukey(seq(-ws,ws)/(ws+1));filter=filter/sum(filter)
  }
  ##this function assumes genome position of mat[Indexes[[1]] are ordered
  Indexes=split(seq(along=pns),pns)
  for(i in seq(along=Indexes)){
    #if(verbose) if(i%%1000==0) cat(i,",")
    Index=Indexes[[i]]
    for(j in 1:ncol(mat)){
      mat[Index,j] <- myfilter(mat[Index,j],filter)
    }
  }
  return(mat)
}

comp = function(g){
    g = sort(unique(g))
    ng = length(g)
    v = vector("character",length=ng*(ng-1)/2)
    count=1
    for(j in 1:(ng-1)){
        for(k in (j+1):ng){
            v[count]   = g[j]
            v[count+1] = g[k]
            count = count+2
        }
    }
    return(v)
}
 
get.tog <- function(l,groups,compare,verbose){
  #require(genefilter)
  
  gIndex=split(seq(along=groups),groups)
  gIndex=gIndex[which(names(gIndex)%in%compare)]
  ng=length(gIndex)
  
  lm=matrix(0,nrow(l),ng)
  ls=matrix(0,nrow(l),ng)
  ns=sapply(gIndex,length)
  
  if(verbose) message("Computing group medians and SDs for",ng,"groups:")
  for(i in seq(along=gIndex)){
    if(verbose) message("\n",i)
    Index=gIndex[[i]]
    if(length(Index)>1){
        lm[,i]=rowMedians(l[,Index,drop=FALSE])
        ls[,i]=rowMads(l[,Index,drop=FALSE],center=lm[,i])
    } else{
        message(paste(" ",names(gIndex)[i],"has only 1 array!"))
        lm[,i]=l[,Index,drop=FALSE]
        ls[,i]=NA
    }
  }
  colnames(lm)=names(gIndex)
  colnames(ls)=names(gIndex)

  nums  <- match(compare,colnames(lm))
  COMPS <- matrix(nums,ncol=2,byrow=TRUE)
  
  if(verbose) message("\nDone.\n")
  return(list(lm=lm,ls=ls,ns=ns,COMPS=COMPS))
}

get.tt <- function(lm,ls,ns,filter,Indexes,COMPS,ws,verbose){
##this function assumes genome position of mat[Indexes[[1]] are ordered:
  if(is.null(filter)){
    Tukey = function(x) pmax(1 - x^2,0)^2
    filter= Tukey(seq(-ws,ws)/(ws+1));filter=filter/sum(filter)
  }
  if(!isTRUE(all.equal(sum(filter),1))) stop("filter must sum to 1.")
  
  if(verbose) message("Smoothing:\n")
  dm=matrix(0,nrow(lm),nrow(COMPS))
  cnames = vector("character",nrow(COMPS))
  for(r in 1:ncol(dm)) cnames[r]=paste(colnames(lm)[COMPS[r,]],collapse="-")
  colnames(dm) = cnames
  vr = dm
  sdm = dm
  svr = dm
  tt = dm
  for(r in 1:nrow(COMPS)){
      j=COMPS[r,1]
      k=COMPS[r,2]
      dm[,r]=lm[,j]-lm[,k]
      vr[,r]=(ls[,j]^2)/ns[j]+(ls[,k]^2)/ns[k]
  }
  if(verbose) pb = txtProgressBar(min=1, max=length(Indexes), initial=0, style=3)
  for(i in seq(along=Indexes)){
    #if(verbose) if(i%%1000==0) cat(".")
    Index=Indexes[[i]]
    for(r in 1:nrow(COMPS)){
        j=COMPS[r,1]
        k=COMPS[r,2]
        sdm[Index,r]=myfilter(dm[Index,r],filter)
        if(ns[j]==1|ns[k]==1){ svr[Index,r] = 1} else{ svr[Index,r] = myfilterse(vr[Index,r],filter) }
    }
    if(verbose) setTxtProgressBar(pb, i)
  }
  if(verbose) close(pb)
  for(r in 1:nrow(COMPS)) tt[,r] = sdm[,r]/svr[,r]
  if(verbose) message("Done.")
  return(tt)
}


dmrFinder <- function(eset=NULL, groups, p=NULL, l=NULL, chr=NULL, pos=NULL, pns=NULL, sdBins=NULL, controlIndex=NULL, controlProbes=c("CONTROL_PROBES", "CONTROL_REGIONS"), Indexes=NULL, filter=NULL, package=NULL, ws=7, verbose=TRUE, compare="all",  withinSampleNorm="loess", betweenSampleNorm="quantile", cutoff=0.995, sortBy="ttarea",...){

  groups = as.character(groups)
  if(identical(compare,"all")) compare=comp(groups)
  if(length(compare)%%2!=0) stop("compare must have an even number of elements.")
  if(length(cutoff)==1) cutoff <- rep(cutoff, length(compare)/2)
  if(length(compare)/2!=length(cutoff)) stop(length(compare)/2," comparisons requested but ", length(cutoff)," cutoff(s) provided.")

  args=list(filter=filter, ws=ws, betweenSampleNorm=betweenSampleNorm, 
	    withinSampleNorm=withinSampleNorm, sdBins=sdBins,
            controlProbes=controlProbes, cutoff=cutoff, sortBy=sortBy)

  # dmrFinder must be given either eset or p/l,chr,pos,pns, and controlIndex
  # If eset is supplied then all the latter will be taken from it (with any
  # that were given as arguments being ignored) (except p).
  # l=logit(p). Indexes=split(seq(along=pns),pns).
  if(is.null(eset)){
      if (is.null("p") & is.null("l")) stop("p or l must be supplied.")
      Args = c("chr","pos","pns","controlIndex","package")
      nulls = sapply(Args,function(x) is.null(get(x)))
      if(any(nulls))
        stop(paste("The following arguments are missing:", paste(Args[nulls], collapse=", ")))
      lens = c( nrow(p), nrow(l), length(chr), length(pos), length(pns) ) 
      if(length(unique(lens))!=1)
        stop("p, l, chr, pos, and/or pns are incompatible.")
      stopifnot(length(groups)==max(ncol(p), ncol(l)))
	  index <- which(!is.na(chr) & !is.na(pos) & !is.na(pns))
      index=index[order(chr[index],pos[index])]
      chr=chr[index]
      pos=pos[index]
      pns=pns[index]
      controlIndex=which(index%in%controlIndex)
      if(!is.null(sdBins)) sdBins<-sdBins[index]
      if(!is.null(p)) p=p[index,]
      if(!is.null(l)) l=l[index,]
  } else if (is.character(eset)) {
          if (is.null("p") & is.null("l")) stop("p or l must be supplied.")
	  pdInfo=get(eset)
	  class(pdInfo)="TilingFeatureSet" # Trick oligo so that pmChr, pmPosition, probeNames work
	  chr=pmChr(pdInfo)
	  pos=pmPosition(pdInfo)
	  if (!is.null(l)) {
		index=which(rowSums(is.na(l))==0)	
        index=index[order(chr[index],pos[index])]
		l=l[index,]
	  } else {
		index=which(rowSums(is.na(p))==0)			
        index=index[order(chr[index],pos[index])]
        p=p[index,]
	  }
      chr=chr[index]
      pos=pos[index]
      pns=probeNames(pdInfo)[index]
	  if (is.null(controlIndex)) {
	  	controlIndex <- which(getContainer(pdInfo) %in% controlProbes)
	  }
      controlIndex=which(index%in%controlIndex)
      if(!is.null(sdBins)) sdBins<-sdBins[index]
      package = eset
      if(package=="pd.feinberg.mm8.me.hx1"){
      #add some code here to break up regions with gaps of >300 bp
      }
      if(package=="pd.feinberg.hg18.me.hx1"){
      #add some code here to break up regions with gaps of >300 bp
      }
  } else {
      stopifnot(length(groups)==ncol(eset))
      if(is.null(p) & is.null(l)){
		p <- methp(eset, 
					withinSampleNorm=withinSampleNorm, 
					betweenSampleNorm=betweenSampleNorm,
					controlIndex=controlIndex,
					controlProbes=controlProbes, 
					verbose=TRUE)
	  }
      chr=pmChr(eset)
      pos=pmPosition(eset)
	  if (!is.null(l)) {
		index=which(rowSums(is.na(l))==0)	
        index=index[order(chr[index],pos[index])]
		l=l[index,]
	  } else {
		index=which(rowSums(is.na(p))==0)			
        index=index[order(chr[index],pos[index])]
        p=p[index,]
	  }
      chr=chr[index]
      pos=pos[index]
      pns=probeNames(eset)[index]
	  if (is.null(controlIndex)) {
        controlIndex=getControlIndex(eset, controlProbes=controlProbes)
      }
      controlIndex=which(index%in%controlIndex)
      if(!is.null(sdBins)) sdBins<-sdBins[index]
      package = annotation(eset)
      if(package=="pd.feinberg.mm8.me.hx1"){
      #add some code here to break up regions with gaps of >300 bp
      }
      if(package=="pd.feinberg.hg18.me.hx1"){
      #add some code here to break up regions with gaps of >300 bp
      }
  }
  if(is.null(l)) {
	  l=log(p)-log(1-p)	
  }
  if(is.null(Indexes)) {
  	  Indexes=split(seq(along=pns),pns)
  }
  tog = get.tog(l=l,groups=groups,compare=compare,verbose=verbose)
  lm=tog$lm
  ls=tog$ls
  ns=tog$ns
  COMPS=tog$COMPS
  tt = get.tt(lm=lm,ls=ls,ns=ns,COMPS=COMPS,Indexes=Indexes,filter=filter,ws=ws,verbose=verbose)

  res=vector("list",ncol(tt))
  names(res)=colnames(tt)
  if(verbose) message("Finding DMRs for each pairwise comparison.")
  for(r in 1:nrow(COMPS)){
      j = COMPS[r,1]
      k = COMPS[r,2]
      if(verbose) message("\n",colnames(tt)[r])
      DF=ifelse(ns[j]==1 & ns[k]==1, 1, ns[j]+ns[k]-2)
	
	  if (length(sdBins)==0) {
	      K=mad(tt[,r], na.rm=TRUE)*qt(cutoff[r],DF)	
	  }	else {
		  s <- tapply(tt[,r], sdBins, mad, na.rm=TRUE)
		  K=s[sdBins]*qt(cutoff[r],DF)	
	  }
      LAST=0
      segmentation=vector("numeric",nrow(tt))
      type=vector("numeric",nrow(tt))
      for(i in seq(along=Indexes)){
        if(verbose) if(i%%1000==0) message(".", appendLF=FALSE)
        Index=Indexes[[i]]
        y=tt[Index,r]
		if(length(sdBins)==0) {
			tmp=sign(y)*as.numeric(abs(y)>K)	
		} else {
			Ki <- K[Index]
			tmp=sign(y)*as.numeric(abs(y)>Ki)	
		}
        tmp2=cumsum(c(1,diff(tmp)!=0))+LAST
        segmentation[Index]=tmp2
        type[Index]=tmp
        LAST=max(tmp2)
      }
      
      Index=which(type!=0)
	  rows <- length(unique(segmentation[Index]))
      res[[r]]=data.frame(
           chr=tapply(chr[Index],segmentation[Index],function(x) x[1]),
           start=tapply(pos[Index],segmentation[Index],min),
           end=tapply(pos[Index],segmentation[Index],max),
           p1=rep(NA, rows),
           p2=rep(NA, rows),
           regionName=tapply(pns[Index],segmentation[Index],function(x) x[1]),
           indexStart=tapply(Index,segmentation[Index],min),
           indexEnd=tapply(Index,segmentation[Index],max))
      length=res[[r]]$indexEnd-res[[r]]$indexStart+1
	  if (is.null(p)) { #  We return log-ratios
	      colnames(res[[r]]) <- sub("p1", "m1", colnames(res[[r]]))
	      colnames(res[[r]]) <- sub("p2", "m2", colnames(res[[r]]))
	      res[[r]]$m1=tapply(lm[Index,j],segmentation[Index],mean)
              res[[r]]$m2=tapply(lm[Index,k],segmentation[Index],mean)
	      area=abs(res[[r]]$m2-res[[r]]$m1)*length
	      res[[r]]$area=area
	  } else { # We return percentages
	      res[[r]]$p1=tapply(1/(1+exp(-lm[Index,j])),segmentation[Index],mean)
              res[[r]]$p2=tapply(1/(1+exp(-lm[Index,k])),segmentation[Index],mean)
	      area=abs(res[[r]]$p2-res[[r]]$p1)*length
	      res[[r]]$area=area
	  }
      ttarea = abs(tapply(tt[Index,r],segmentation[Index],mean)) *length
      res[[r]]$ttarea = ttarea
      if(sortBy=="area")   res[[r]]=res[[r]][order(-area),]
      if(sortBy=="ttarea") res[[r]]=res[[r]][order(-ttarea),]

  }
  if(verbose) message("\nDone")
  return(list(tabs=res, p=p, l=l, chr=chr, pos=pos, pns=pns, 
              index=index, controlIndex=controlIndex, gm=lm,
              groups=groups, args=args, comps=COMPS, package=package))
}

dmrFdr <- function(dmr, compare=1, numPerms=1000, seed=NULL, verbose=TRUE) {
	if (length(compare)!=1) stop("You must choose one comparison at a time when calculating FDRs. Please set dmr to be one of: ", 
	paste(names(dmr$tabs), collapse=", "), "\n")
	if (is.numeric(compare)) compare <- names(dmr$tabs)[compare]
	message("Calculating q-values for DMRs between", compare, "\n")
	# Get probe order from TilingFeatureSet object
	pdInfo=get(dmr$package)
	class(pdInfo)="TilingFeatureSet" # Trick oligo so that pmChr, pmPosition work
	chr=pmChr(pdInfo)
	pos=pmPosition(pdInfo)
	chrpos <- paste(chr, pos)
	dchrpos <- paste(dmr$chr, dmr$pos)
	#idx <- which(!(chrpos %in% dchrpos))
	mis <- setdiff(chrpos, dchrpos)
	o <- order(chr, pos)
	chrpos <- chrpos[o]
	keepProbes <- !(chrpos %in% mis)
	o <- o[keepProbes]
	keep <- dmr$groups %in% unlist(strsplit(compare, "-"))
	# Recreate p or l with same sort order as in annotation package 
	if (is.null(dmr$p)) {
		l <- matrix(NA, nrow=length(pos), ncol=ncol(dmr$l))
		l[o,] <- dmr$l
		l <- l[,keep]
	} else {
		p <- matrix(NA, nrow=length(pos), ncol=ncol(dmr$p))
		p[o,] <- dmr$p
		p <- p[,keep]
	}
	n <- sum(keep)
	n1 <- sum(dmr$groups==unlist(strsplit(compare, "-"))[1])
	maxPerms <- choose(n, n1)
	if (numPerms=="all") numPerms <- maxPerms
	if (numPerms>maxPerms) {
		message("Given the sample sizes in the two groups the maximum number of permutations is ", maxPerms, sep="")
		numPerms <- maxPerms
	} 

	## Reshuffled group label DMRs
	if (!is.null(seed)) set.seed(seed)
	s <- sample(1:maxPerms, numPerms)
	grp1 <- combinations(n,n1)[s,]

	if (verbose) message("Finding permuted data DMRs. Estimating time remaining")
	areas <- lapply(1:numPerms, function(i) {
		groups <- rep("grp2", n)
		groups[grp1[i,]] <- "grp1"
		if (is.null(dmr$p)) {
			st <- system.time(dmrPerm <- dmrFinder(dmr$package, l=l, 
				groups=groups, cutoff=dmr$args$cutoff, 
				filter=dmr$args$filter, ws=dmr$args$ws, 
				verbose=FALSE))[3]
		} else {
			st <- system.time(dmrPerm <- dmrFinder(dmr$package, p=p, 
				groups=groups, cutoff=dmr$args$cutoff, 
				filter=dmr$args$filter, ws=dmr$args$ws,
				verbose=FALSE))[3]
		}
		if (verbose & (i %in% round(seq(1, numPerms, length.out=10)))) {
			message(i, "/", numPerms, " (", prettyTime((numPerms-i)*st), 
				" remaining)", sep="")
		}
		dmrPerm$tabs[[1]][,dmr$args$sortBy]
	})

	nullDist <- unlist(areas)
	fn <- ecdf(nullDist)
	pval <- 1-fn(dmr$tabs[[compare]][,dmr$args$sortBy])
	pi0<-pi0.est(pval)$p0
	qval<-qvalue.cal(pval, pi0)
	dmr$tabs[[compare]] <- cbind(dmr$tabs[[compare]], pval, qval)
	if (!("numPerms" %in% names(dmr))) {
		dmr$numPerms <- rep(NA, length(dmr$tabs))
		names(dmr$numPerms) <- names(dmr$tabs)
	}
	dmr$numPerms[compare] <- numPerms
	return(dmr)
}


prettyTime <- function(seconds) {
	hours <- floor(seconds / 60 / 60)
	minutes <- floor((seconds/60) - (hours*60))
	ret <- ""
	if (hours>0) ret <- paste(ret, hours, 
		ifelse(hours==1, " hour", " hours"), ", ", sep="")
	ret <- paste(ret, minutes, ifelse(minutes==1, " minute", " minutes"),
		sep="")
	ret
}


validatePd <- function(pd, fileNameColumn, sampleNameColumn, 
		groupColumn, ut = "_532.xys", md = "_635.xys") {
	msg <- ""
	message("Filenames:")
	if (missing(fileNameColumn)) {
		fileNameColumn <- which.max(apply(pd, 2, function(col) {
			u <- grep(ut, col, ignore.case=TRUE)
			m <- grep(md, col, ignore.case=TRUE)
			length(c(u,m))
		}))
		message("  fileNameColumn not specified. Trying to guess.\n")	
	}
	files <- pd[,fileNameColumn]
	o <- order(files)
	files <- files[o]
	pd <- pd[o,]
	utIdx <- grep(ut, files)
	if (length(utIdx)==0) {
		msg <- paste(msg, "  No files in column ", fileNameColumn, " match the untreated extension ", ut, ". Please use the ut option to set the correct extension\n", sep="")
	}
	mdIdx <- grep(md, files)
	if (length(mdIdx)==0) {
		msg <- paste(msg, "  No files in column ", fileNameColumn, " match the methyl-depleted extension ", md, ". Please use the md option to set the correct extension\n", sep="")
	}
	filesUt <- sub(ut, "", files[utIdx])
	filesMd <- sub(md, "", files[mdIdx])
	missing <- c(filesUt[!filesUt%in%filesMd], filesMd[!filesMd%in%filesUt])
	if (length(missing)>0) {
		msg <- paste(msg, "  The untreated (ut) and methyl-depleted (md) file names in column ", fileNameColumn, " don't match up. Check ", paste(missing, collapse=", "), sep="")
	}
	if (length(missing)>0 | length(utIdx)==0 | length(mdIdx)==0) {
		message("  ERROR - ", msg)
		return(list(status="ERROR", message=msg))	
	} else {
		message("  OK - Found in column", fileNameColumn)
	}

	channelMatch <- apply(pd[utIdx,]==pd[mdIdx,], 2, all)
	numUnique <- apply(pd[utIdx,], 2, function(x) length(unique(x)))
	onePerSample <- numUnique==length(utIdx)

	message("Sample names:")
	if (missing(sampleNameColumn)) {
		message("  sampleNameColumn column not specified. Trying to guess.")	
		sampleNameColumn <- which(channelMatch & onePerSample)	
		if (length(sampleNameColumn)==1) {
			message("  OK - Found in column", sampleNameColumn)
		} else if (length(sampleNameColumn)==0) {
			msg <- paste(msg, "No suitable sample ID column found. This column should have the same entry for each pair of untreated and methyl-depleted filenames.\n")
			message("  ERROR - ", msg, appendLF=FALSE )
			return(list(status="ERROR", message=msg))
		} else if (length(sampleNameColumn)>1) {
			message("  WARNING - Multiple columns (", paste(sampleNameColumn, collapse=", "), 
				") are candidates for sample ID.", sep="")
				sampleNameColumn <- sampleNameColumn[1]
		}
	} else { # User specified sampleNameColumn
		if (channelMatch[sampleNameColumn] & onePerSample[sampleNameColumn]) {
			message("  OK - Found in column", sampleNameColumn)
		} else {
			if (!channelMatch[sampleNameColumn]) {
				msg <- paste(msg, "  There are samples with different values in column ", sampleNameColumn, " for their untreated and methyl-depleted rows\n", sep="")
			} else if (!onePerSample[sampleNameColumn]) {
				msg <- paste(msg, "  The sample IDs in column ", sampleNameColumn, " are not unique to each sample\n", sep="")
			}
			message("  ERROR -", msg)
			return(list(status="ERROR", message=msg))	
		}
	}
	message("Group labels:")
	if (missing(groupColumn)) {
		message("  groupColumn column not specified. Trying to guess.")	
		groupColumn <- which(channelMatch & !onePerSample & numUnique>1)
		if (length(groupColumn)==1) {
			message("  OK - Found in column", groupColumn)
		} else if (length(groupColumn)==0) {
			msg <- paste(msg,  "No suitable group label column found. If there are no replicates you can set groupColumn to be the same as sampleNameColumn\n")
			message("  ERROR -", msg, appendLF=FALSE)
			return(list(status="ERROR", message=msg))	
		} else if (length(groupColumn)>1) {
			warning("  WARNING - Multiple columns (", 
				paste(groupColumn, collapse=", "), 
				") are candidates for group labels.\n", sep="")
				groupColumn <- groupColumn[1]
		}
	} else { # User specified groupColumn
		if (channelMatch[groupColumn] & !onePerSample[groupColumn]) {
			message("  OK - Found in column", groupColumn)
		} else {
			if (!channelMatch[groupColumn]) {
				msg <- paste(msg,  "There are samples with different values in column", groupColumn, "for their untreated and methyl-depleted rows\n")
				message("  ERROR -", msg, appendLF=FALSE)
				return(list(status="ERROR", message=msg))	
			} else if (onePerSample[groupColumn]) {
				message("  WARNING - Each group has only 1 sample")
			}
		}
	}
	return(list(status="OK", fileNameColumn=fileNameColumn, 
			sampleNameColumn=sampleNameColumn,
			groupColumn=groupColumn))
}


.onAttach <- function(libname, pkgname) {
 message("Welcome to charm version ", packageDescription("charm", field="Version"))
}
