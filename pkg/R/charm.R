## Temporary fix for sampleNames
#sampleNames <- function(dat) {
#    rownames(pData(dat))
#}


methp <- function(dat, spatial=TRUE, spatialMethod="kernel", 
		spatial1d=NULL, spatial2d=NULL, 
		bgSubtract=TRUE,
		withinSampleNorm="loess", binSize=500, numSegments=3, useTot=TRUE,
		betweenSampleNorm="quantile", msqn=FALSE,
        minQCScore=NULL, 
		controlProbes=c("CONTROL_PROBES", "CONTROL_REGIONS"),
		controlIndex=NULL, 
		ctrlNoCpGWindow=NULL, subject=NULL,
		commonMethPercentParams=NULL,
		verbose=TRUE, cluster=NULL, returnM=FALSE, 
		plotDensity=NULL, plotDensityCols=NULL,
		duplicate=TRUE) {
	
	if (!is.null(spatial1d) | !is.null(spatial2d)) {
		warning("The spatial1d and spatial2d options are deprecated. Please use spatial instead\n")
	}
    if (!is.null(minQCScore)) {
        warning("methp no longer filters out probes based on quality score. Please use the pmQuality function to calculate probe quality.\n")
    }	
			
	if (class(dat)=="OfflineTilingFeatureSet2" & duplicate) {
		dat <- clone(dat)
	}		
			
	# Start parallel processes
	if (is.null(cluster)){
		#cl <- charmCluster()
		cl <- NULL
	} else if (is.numeric(cluster)) {
		cl <- charmCluster(cluster)
	} else if (any(class(cluster)=="cluster")){
		cl <- cluster
	}
	
	if(!is.null(plotDensity)) {
		pdf(file=plotDensity, height=11, width=8)
		par(mfrow=c(5,2), mar=c(2,2,4,2))
		lwd <- rep(1, ncol(dat))
		if (is.null(plotDensityCols)) {
			if (!is.null(pData(rawData)$type)){
				cols <- as.numeric(factor(pData(rawData)$type))
			} else {
				cols <- rep(1, ncol(dat))
			}
		} else {
			cols <- plotDensityCols
		}
		plotDensity(dat, main="1. Raw", cols=cols, lwd=lwd)
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
			bs <- list(m="99", untreated="complete", enriched="sqn")
			if(is.null(commonMethPercentParams)) commonMethPercentParams <- FALSE
		} else if (betweenSampleNorm=="none") {
			bs <- list(m="none", untreated="none", enriched="none")
			if(is.null(commonMethPercentParams)) commonMethPercentParams <- FALSE
		}
	}	
	
    # Spatial bias correction
    if (spatial) {
        if (verbose) cat("Spatial normalization ")
		if (spatialMethod=="kernel") {
			cat ("\n")
        	dat <- spatialAdjust(dat, cluster=cl)
		} else if (spatialMethod=="poly") {
			cat ("(polynomial surface)\n")
			dat <- spatialAdjust.poly(dat, cluster=cl)	
		} else {
			stop("Invalid value for spatialMethod.\n")
		}
    }

    # Background removal
	if (bgSubtract) {
    	if (verbose) cat("Background removal\n")
		dat <- bgAdjustBgp(dat, cluster=cl)
	}
	if(!is.null(plotDensity)) {
		plotDensity(dat, main="2. After spatial & bg", cols=cols, lwd=lwd)
	}	
	# Within sample normalization
	if (is.null(controlIndex)) {
		controlIndex <- getControlIndex(dat, noCpGWindow=ctrlNoCpGWindow, subject=subject, controlProbes=controlProbes)
	}
	if (verbose) cat("Within sample normalization: ", withinSampleNorm, "\n", sep="") 
	dat <- normalizeWithinSamples(dat, method=withinSampleNorm, useTot=useTot,
		binSize=binSize, numSegments=numSegments, cluster=cl,
		controlIndex=controlIndex, verbose=verbose)
    # Between sample normalization    
    #if (verbose) cat("Between sample normalization. 
	if(!is.null(plotDensity)) {
		plotDensity(dat, main="3. After within-sample norm", cols=cols, lwd=lwd)
	}
	if (verbose) {
		cat("Between sample normalization")
		if (is.list(betweenSampleNorm)) {
			cat(". M: ", bs$m, ", Untreated channel: ", bs$untreated, 
		        ", Methyl-depleted channel: ", bs$enriched, "\n", sep="")
		} else {
			cat (": ", betweenSampleNorm, "\n", sep="")
		}
	}
    dat <- normalizeBetweenSamples(dat, m=bs$m, 
		untreated=bs$untreated, 
        enriched=bs$enriched, controlProbes=controlProbes, 
		controlIndex=controlIndex, cluster=cl, verbose=verbose)
    if (any(class(dat)=="OfflineTilingFeatureSet2")) {
		M <- getPmM(dat)
	} else {
		M <- getM(dat)[pmindex(dat),]
	}
		
	if (msqn) {
		cat("The msqn option is not implemented\n")
		#cat("Between sample normalization: msqn\n")
		#M <- parSQN(y=M[,], ctrl.id=controlIndex, cluster=cl)
	}	
		
	if(!is.null(plotDensity)) {
		plotDensity(dat, main="4. After between-sample norm", cols=cols, lwd=lwd)
	}		
	
	if (returnM=="TRUE" | returnM=="+") {
		retval <- M
	} else if (returnM=="-") {
	    if (any(class(dat)=="OfflineTilingFeatureSet2")) {
			for (i in 1:ncol(M)) {
				M[,i] <- -M[,i]
				retval <- M
			}
		} else {
			retval <- -M
		}
	} else {
	    if (verbose) cat("Estimating percentage methylation\n")
     	retval <- methPercent(m=M, commonParams=commonMethPercentParams,
	 		ngc=countGC(dat), cluster=cl)
	}
	if(!is.null(plotDensity)) {
		if(returnM=="FALSE") plotDensity(retval, main="5. Percentage methylation", controlIndex=getControlIndex(dat), rx=c(0,1), cols=cols, lwd=lwd)
		dev.off()		
	}
   	if (all(class(cluster)!="cluster")) stopCluster(cl) 
    return(retval)
}


readCharm <- function(files, type=rep("unspecified", length(files)), path=".", ut="_532.xys", md="_635.xys", sampleNames=NULL, saveRam=FALSE, ...) {
    files <- as.character(files)
    o <- order(files)
    files <- files[o]
    type <- as.character(type[o])
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
        stop(cat("The untreated (ut) and methyl-depleted (md) file names don't match up\n"))
    if (!all(type[utIdx] == type[mdIdx])) 
        stop(cat("The untreated (ut) and methyl-depleted (md) type labels don't match up\n"))
    if (!is.null(sampleNames)) {
        sampleCheck <- sampleNames[utIdx] == sampleNames[mdIdx]
        if (!all(sampleCheck)) 
            stop(cat("The untreated (ut) and methyl-depleted (md) sample names don't match up\n Check:", 
                sampleNames[utIdx][!sampleCheck], "\n"))
        sampleNames <- sampleNames[utIdx]
    } else {
        sampleNames <- sub(ut, "", filesUt)
    }
    pd <- data.frame(arrayUT=filesUt, arrayMD=filesMd, type=type[utIdx], stringsAsFactors=FALSE)
    vpd <- data.frame(labelDescription= c("Untreated (UT)", "Methyl-depleted (MD)", "Sample type (e.g. tissue type, cancer/normal, etc.)"),
                     channel=factor(c("channel1", "channel2", NA), levels=c("channel1", "channel2", "_ALL_")))
    pdd <- new("AnnotatedDataFrame", data=pd, varMetadata=vpd)
    sampleNames(pdd) <- sampleNames  
    if (saveRam) {
		dat <- read.xysfiles2.offline(channel1=file.path(path, filesUt), 
			channel2=file.path(path, filesMd),
        	phenoData=pdd, sampleNames=sampleNames, ...)
	} else {
    	dat <- read.xysfiles2(channel1=file.path(path, filesUt), 
			channel2=file.path(path, filesMd),
        	phenoData=pdd, sampleNames=sampleNames, ...)     
	}
	return(dat)
}

charmCluster <- function(cluster=NULL, type="SOCK", verbose=FALSE) {
	if(any(class(cluster)=="cluster")) { 
		return(cluster)
	} else {
		if(is.numeric(cluster)) {
			if (verbose>1) {
				cat("Starting", cluster, "node", type, "cluster\n")
			}
			cl <- makeCluster(cluster, type=type)			
		} else {
			#cat("Starting single node", type, "cluster\n")
			cl <- makeCluster(1, type=type)			
		}
		clusterEvalQ(cl, library(charm))
		return(cl)
	}
}

plotDensity <- function(dat, rx=c(-4,6), controlIndex=NULL, 
		pdfFile=NULL, main=NULL, cols=NULL, lwd=NULL) {
	if (!is.null(pdfFile)) {
		pdf(pdfFile)
		par(mfcol=c(2,1))
	}
	if (any(class(dat)=="OfflineTilingFeatureSet2")) {
		M <- getPmM(dat)
		lab <- sampleNames(dat)
	} else 	if (any(class(dat)=="TilingFeatureSet2")) {
		M <- getM(dat)[pmindex(dat),,drop=FALSE]
		lab <- sampleNames(dat)
	} else {
		M <- dat[,]
		lab <- colnames(dat)
	}
	if (is.null(cols)) cols <- 1:ncol(M)
	if (is.null(lwd)) lwd <- rep(1, ncol(M))
	
	if (is.null(controlIndex)) controlIndex <- getControlIndex(dat)
	plotDensityMat(M[,], xlab="M", lab=lab, 
		main=paste(main,"\nAll probes"), rx=rx, cols=cols, lwd=lwd)
	plotDensityMat(M[controlIndex,], xlab="M", lab=lab, 
		main=paste(main, "\nControl probes"), rx=rx, cols=cols, lwd=lwd)
	if (!is.null(pdfFile)) dev.off()
}

normalizeLoess <- function(dat, controlIndex=NULL, 
		controlProbes=c("CONTROL_PROBES", "CONTROL_REGIONS"), span=0.3,
		by="A", approx=TRUE, breaks=1000, cluster=NULL) {
	if (is.null(controlIndex)) {
    	controlIndex <- getControlIndex(dat, controlProbes=controlProbes)
	}
	pms <- pm(dat)
	if (class(dat)=="OfflineTilingFeatureSet2") {
		cl <- charmCluster(cluster)
		M <- getPmM(dat)
		clRet <- parLapply(cl, 1:ncol(M), function(i, pms, M, by, controlIndex) {
			y <- M[,i]
			c1 <- log2(pms[,i,"channel1"])
			c2 <- log2(pms[,i,"channel2"])
			if (by=="tot") {
				x <- c1
			} else if (by=="A") {
				x <- (c1+c2)/2
			}
	 		fit <- loess(y ~ x, subset = controlIndex,
	               na.action = na.exclude, degree=1, surface = "direct",
	 			   span=span)
			adj <- predictLoess(fit, newdata=x, approx=approx, breaks=breaks)
			ret <- 2^(c2 + adj)
			pms[, i, "channel2"] <- ret
			return(NULL)
		}, pms, M, by, controlIndex)
		if (all(class(cluster)!="cluster")) stopCluster(cl) 
	} else {
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
	}
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
	


affinityBin <- function(dat, M=NULL, controlIndex=NULL, binSize=500,
		modelChannel="M", controlProbes=c("CONTROL_PROBES", "CONTROL_REGIONS"),
		numSegments=3, window=NULL, subject=NULL, useTot=TRUE, verbose=FALSE) {
	if (is.null(controlIndex)) {
		controlIndex <- getControlIndex(dat, controlProbes=controlProbes)   
    }
	affinity <- getAffinity(dat, channel=modelChannel, idx=controlIndex,
		numSegments=numSegments, window=window, subject=subject, 
		useTot=useTot, verbose=verbose)                       
	affinityCtrl <- affinity[controlIndex]
    bin <- rep(1:length(affinityCtrl), length.out=length(affinityCtrl), each=binSize)
	o <- order(affinityCtrl)
	binMedianAffinity <- tapply(affinityCtrl[o], bin, median) 
	affinityBin <- sapply(affinity, function(x) which.min(abs(x-binMedianAffinity)))
	n <- ifelse(is.null(M), dim(dat)["Samples"], ncol(M))
	#if (is.null(M) & class(dat)=="TilingFeatureSet2") {
	#	M <- getM(dat)[pmindex(dat),,drop=FALSE]
		#if (!is.matrix(M)) M <- as.matrix(M)
	#} 
	affinityBinMed <- matrix(NA, nrow=length(unique(bin)), ncol=n)
	affinityBinSd <- matrix(NA, nrow=length(unique(bin)), ncol=n)
	if (is.null(M)){
		pms <- log2(pm(dat)[controlIndex,,])
	}
	for (i in 1:n) {
		if (is.null(M)) {
			m <- pms[,i,1] - pms[,i,2]
		} else {
			m <- M[controlIndex,i]
		}
		affinityBinMed[,i] <- tapply(m[o], bin, median) 
		affinityBinSd[,i] <- tapply(m[o], bin, mad) 
	}        
	return(list(bin=affinityBin, binMed=affinityBinMed, binSd=affinityBinSd))
}

regionFilter <- function(x, region, f) {
	x<-unlist(tapply(x, region, myfilter, f))
	names(x) <- NULL
	return(x)
}

normalizeWithinSamples <- function (dat, method = "loess", 
	controlProbes = "CONTROL_REGIONS", 
	controlIndex = NULL, numSegments=3, bins=NULL, binSize=500,
	approx=TRUE, breaks=1000, useTot=TRUE,
	cluster=NULL, verbose=FALSE) {
    #if (method=="gc") {
	if (grepl("gc", method)) {
		if (class(dat)=="OfflineTilingFeatureSet2") {
			stop("normalizeWithinSamples does not yet support the method='gc' option for OfflineTilingFeatureSet2 objects. Please use the method='affinity' option.")
		}
        mAdj <- diffAmpEstGC(dat, method = method, controlProbes = controlProbes)
        dat <- diffAmpAdjustGC(dat, mAdj)
    }
	if (grepl("loess", method)) {
			dat <- normalizeLoess(dat, controlIndex=controlIndex, 
				controlProbes=controlProbes, 
				approx=approx, breaks=breaks, cluster=cluster)
	}
	if (grepl("affinity", method)) {   
        if (is.null(controlIndex)) {
            controlIndex <- getControlIndex(dat, controlProbes=controlProbes)
        }
		if (is.null(bins)) {
			bins <- affinityBin(dat, controlIndex=controlIndex, 
				binSize=binSize, useTot=useTot, verbose=verbose)
		}
		pms <- pm(dat)
		for (i in 1:ncol(pms)) {	
			mAdj <- bins$binMed[bins$bin,i]
               pms[, i, "channel2"] <- 2^(log2(pms[, i, "channel2"]) + mAdj)
           }        
		pm(dat) <- pms        
    }
	if (grepl("median", method)) {
		if (class(dat)=="OfflineTilingFeatureSet2") {
			stop("normalizeWithinSamples does not yet support the method='median' option for OfflineTilingFeatureSet2 objects. Please use the method='affinity' option.")
		}
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
    return(dat)
}

normalizeBetweenSamples <- function(dat, m="allQuantiles", untreated="none", enriched="none", 
		mapping="affinity", numSegments=3, controlProbes="CONTROL_REGIONS", 
		controlIndex=NULL, cluster=NULL, verbose=FALSE) {    
    #datPm <- log2(pm(dat)) 
	pms <- pm(dat)
    if (class(dat)=="OfflineTilingFeatureSet2") {
		M <- getPmM(dat)
	} else {
		M <- getM(dat)[pmindex(dat),,drop=FALSE]
	}
	#M <- datPm[,,"channel1"] - datPm[,,"channel2"]
	#if (!is.matrix(M)) M <- as.matrix(M)
	
	if (m!="none"){
		M <- normQuantile(M, method=m)
		for (i in 1:ncol(M)) {
			pms[,i,"channel2"] <- 2^(log2(pms[,i,"channel1"]) - M[,i])	
			#datPm[,,"channel2"] <- datPm[,,"channel1"] - M	
		}		
	}
    ## Normalize untreated (leaving M unchanged)
    if (untreated=="complete") {
		if (any(class(pms)=="ff")) {
			m <- ffrowapply(rowMedians(log2(pms[i1:i2,,"channel1"]), na.rm=TRUE), X=pms, 
				RETURN=TRUE, RETCOL=NULL, BATCHSIZE=5000)[]
		} else {
        	m <- rowMedians(log2(pms[,,"channel1"]))
		}
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
		if (any(class(pms)=="ff")) {
			pms <- SQN.ff(pms, channel=2, ctrl.id=controlIndex, cluster=cluster)
		} else {
			if (is.null(cluster)) {
				pms[,,"channel2"] <- 2^(SQN(y=log2(pms[,,"channel2"]), ctrl.id=controlIndex))
			} else {
				pms[,,"channel2"] <- 2^(parSQN(y=log2(pms[,,"channel2"]), ctrl.id=controlIndex, cluster=cluster))
			}
		}
	} else {
        pms[,,"channel2"] <- 2^(normQuantile(log2(pms[,,"channel2"]), method=enriched))
    }
    pm(dat) <- pms
    return(dat)
} 

SQN.ff <- function (y, channel=2, N.mix = 5, ctrl.id, model.weight = 0.9, cluster=NULL) {
	cl <- charmCluster(cluster)
    QE = apply(log2(y[ctrl.id, ,channel]), 2, sort)
    QN = apply(QE, 1, median)
    mix.param = Mclust(QN, G = N.mix)$parameters
    mix.param = norMix(mu = mix.param$mean, sig2 = mix.param$variance$sigmasq, 
        w = mix.param$pro)
    qq = seq(1/(2 * length(QN)), 1 - 1/(2 * length(QN)), 1/length(QN))
    qq = qnorMix(qq, mix.param)
    QN1 = QN * (1 - model.weight) + qq * model.weight

	#for (i in 1:ncol(y)) {
	parLapply(cl, 1:ncol(y), function(i, y, channel, ctrl.id, QN1, mix.param) {
		tmp <- mix.qn(log2(y[,i,channel]), ctrl.id, NQ = QN1, mix.param = mix.param, 
	        max.q = 0.95, low = quantile(QN1, 0.05))
		y[,i,channel] <- 2^tmp
	}, y, channel, ctrl.id, QN1, mix.param)
	if (all(class(cluster)!="cluster")) stopCluster(cl) 
	return(y)
}

parSQN <- function (y, N.mix = 5, ctrl.id, model.weight = 0.9, cluster=NULL) {
	cl <- charmCluster(cluster)
    QE = apply(y[ctrl.id, ], 2, sort)
    QN = apply(QE, 1, median)
    mix.param = Mclust(QN, G = N.mix)$parameters
    mix.param = norMix(mu = mix.param$mean, sig2 = mix.param$variance$sigmasq, 
        w = mix.param$pro)
    qq = seq(1/(2 * length(QN)), 1 - 1/(2 * length(QN)), 1/length(QN))
    qq = qnorMix(qq, mix.param)
    QN1 = QN * (1 - model.weight) + qq * model.weight
    ynorm = parApply(cl, y, 2, mix.qn, ctrl.id, NQ = QN1, mix.param = mix.param, 
        max.q = 0.95, low = quantile(QN1, 0.05))
	if (all(class(cluster)!="cluster")) stopCluster(cl) 
	return(ynorm)
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
    cat(" Normalizing with", length(c(mBg, mCtrl)), "control probes ")
    if (!is.null(affinity)) {
        cat("(affinity-based mapping)\n")
    } else {
        cat("(intensity-based mapping)\n")
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



makeX <- function(dat, idx=NULL, window="probe", type="pm", segmentSize=NULL, numSegments=1, letters=c("A", "C", "G"), subject=NULL) {
    if (class(dat)=="DNAStringSet"){
		seqs <- dat
	} else if (type=="bg") {
        #seqs = as.character(bgSequence(dat))
		seqs = bgSequence(dat)
    } else if (type=="pm"){
		if (window=="probe" | is.null(window)) {
        	#seqs = as.character(pmSequence(dat))
			seqs = pmSequence(dat)
		} else if (is.numeric(window)){
			if (is.null(subject)) {
				stop("You must specify a BSGenome object as 'subject' if choosing a window size")
			}
			chr <- pmChr(dat)
			pos <- pmPosition(dat)
			seqs <- rep(NA, length(chr))
			x <- split(1:length(chr), chr)
			for (curchr in (names(x))) {
				#cat(".")
				chrseq <- subject[[curchr]]
				curpos <- pos[x[[curchr]]]
				v <- suppressWarnings(Views(chrseq, start=curpos-window/2,
					 	end=(curpos+(window/2)-1)))
				seqs[x[[curchr]]] <-	as.character(v)
			}	
    	}
	} else {
        stop("You must either provide probe sequences or type must be 'pm' or 'bg'")
    }
    if (!is.null(idx)) seqs <- seqs[idx]
    d <- DNAStringSet(seqs)
    avgLength <-  round(mean(width(d)))
    if(is.null(segmentSize)) segmentSize <- floor(avgLength/numSegments) 
    s <- seq(1, avgLength, segmentSize)
    e <- pmin(s+segmentSize-1, avgLength)
    l <- e-s+1
    n <- length(l)
    if (l[n] < segmentSize) { # last segment is less than segment size
        s <- s[1:(n-1)]
        e <- e[1:(n-1)]
        e[n-1] <- avgLength
    }
    X <- matrix(nrow=length(seqs), ncol=0)
    for (i in 1:length(s)) {
        if (i==length(s)) {
            end <- width(d)
        } else {
            end <- e[i]
        }
        tmp <- alphabetFrequency(narrow(d, s[i], end))
		if (letters[1]=="GC") {
			tmp <- data.frame(gc=tmp[,"C"]+tmp[,"G"])
		} else {
			tmp <- tmp[,letters]
		}
        colnames(tmp) <- paste(colnames(tmp), i, sep="")
        X <- cbind(X, tmp)
    }
    return(as.data.frame(X))
}

addInteractions <- function(X) {
    n <- ncol(X)
    column <- n+1
    for (i in 1:(n-1)) {
        for (j in (i+1):n) {
            #cat(i, j, colnames(X)[i], colnames(X)[j], "\n")    
            X <- cbind(X, X[,i]*X[,j])
            colnames(X)[column] <- paste(colnames(X)[i], colnames(X)[j])
            column <- column+1
        }
    }    
    return(X)
}

makeSeqX <- function(dat, type="pm", idx=NULL, numSegments=3, segmentSize=NULL, MAT=FALSE, addInteractions=FALSE, addSquared=FALSE, window=NULL, subject=NULL) {
    if (MAT) {
        X1 <- makeX(dat=dat, idx=idx, type=type, segmentSize=1, letters=c("A", "C", "G"))
        X2 <- makeX(dat=dat, idx=idx, type=type, segmentSize=25, letters=c("A", "C", "G", "T"))
        colnames(X2) <- c("A.sq", "C.sq", "G.sq", "T.sq")
        X <- cbind(X1, X2^2)        
    } else {
        X <- makeX(dat=dat, idx=idx, type=type, numSegments=numSegments, segmentSize=segmentSize, letters=c("A", "C", "G"))
		if(!is.null(window)) {
			Xw <- makeX(dat=dat, idx=idx, window=window, numSegments=1, subject=subject,  
				letters=c("A", "C", "G", "T"))
			colnames(Xw) <- paste("W", colnames(Xw), sep="")
			X <- cbind(X, Xw[,1:3])
		}
		if(addInteractions) {
			X <- addInteractions(X)
		} else if (addSquared) {
			X1 <- makeX(dat=dat, idx=idx, type=type, numSegments=1, letters=c("A", "C", "G", "T"))
			Xsq <- X1^2
			colnames(Xsq) <- paste(colnames(X1), "sq", sep=".")
			X <- cbind(X, Xsq)
			if (!is.null(window)){
				Xwsq <- Xw^2
				colnames(Xwsq) <- paste(colnames(Xw), "sq", sep=".")
				X <- cbind(X, Xwsq)
			}
		}
    }
    return(as.data.frame(X))
}

getAffinity <- function(dat, type="pm", channel=2, idx=NULL, returnIdx=NULL, exclude=NULL, segmentSize=NULL, numSegments=3, window=NULL, addSquared=TRUE, subject=NULL, MAT=FALSE, useTot=useTot, verbose=FALSE) {
    if(is.null(segmentSize)) {
        d = pmSequence(dat)
        avgLength <-  round(mean(width(d)))
        segmentSize <- floor(avgLength/numSegments)
    }
	if (is.null(idx)) {
		idx <- 1:(dim(pm(dat))[1])
	}
    if (is.matrix(channel)) {
		y <- channel
		y<-y[idx,,drop=FALSE]
	} else {
		if (type=="pm") {
        	if (channel=="M") {
            	if (class(dat)=="OfflineTilingFeatureSet2") {
					y <- getPmM(dat, rows=idx)
				} else {
					y <- getM(dat)[pmindex(dat)[idx], ]
				}
        	} else {
            	y <- log2(pm(dat)[idx,,channel])
        	}
        	if (is.vector(y)) y <- as.matrix(y)
        	#if (!is.null(idx)) y<-y[idx,]
        	#if (is.vector(y)) y <- as.matrix(y)
    	} else if (type=="bg") {
			if (class(dat)=="OfflineTilingFeatureSet2") {
				stop("getAffinity does not yet implement the type='bg' option for OfflineTilingFeatureSet2\n")
			} 
        	if (channel=="M") {
            	y <- getM(dat)[bgindex(dat), ]
        	} else {
            	y <- log2(bg(dat))[,,channel]
        	}
        	if (is.vector(y)) y <- as.matrix(y)
    	}
    	if (!is.null(exclude)) {
        	exclude <- as.matrix(exclude)
        	y[exclude] <- NA
    	}
	}
	yMed <- rowMedians(y, na.rm=TRUE)
	
	X <- makeSeqX(dat, type=type, idx=idx, segmentSize=segmentSize, 
		window=window, subject=subject, addSquared=addSquared)
	#browser()
	#rhs <- paste(c("totMed", colnames(X)), collapse="+")
	rhs <- paste(colnames(X), collapse="+")
	if (useTot) {
		totMed <-  rowMedians(log2(pm(dat)[idx,,1]))
		totMed2 <- totMed^2
		rhs <- paste(c("totMed", "totMed2", rhs), collapse="+")
		#rhs <- paste(c("totMed", rhs), collapse="+")
	}
	
	formul <- formula(paste("yMed~", rhs))
	#w<- -log(1-pmin(0.99,pmQualCtrlM))
	lmfit <- lm(formul, data=X)
	#lmfit <- lm(yMed~.,X)
	#lmfit.base <- lm(formul, data=X)
	#if (addSquared) { # Don't add interactions
	#	lmfit <- step(lmfit.base, trace=0)	
	#} else { # Add interactions
	#	lmfit <- step(lmfit.base, ~ .^2, trace=0)	
	#}
	#cat("  Built control probe affinity model with", length(coef(lmfit)), "coefficients (out of", length(coef(lmfit.base)), "possible). R2 =", round(summary(lmfit)$r.squared,3), "\n")
	if(verbose>1) {
		cat("  Built control probe affinity model with", length(coef(lmfit)), 
			"coefficients. R2 =", 
			round(summary(lmfit)$r.squared,3), "\n")
	}
	#print(summary(lmfit))

    X <- makeSeqX(dat, segmentSize=segmentSize, window=window, subject=subject, addSquared=addSquared)
	if (useTot) {
		if (class(dat)=="OfflineTilingFeatureSet2") {
			pms <- pm(dat)
			totMed <- ffrowapply(rowMedians(log2(pms[i1:i2,,1])), X=pms, 
				RETURN=TRUE, RETCOL=NULL, BATCHSIZE=5000)[]	
		} else {
			totMed <-  rowMedians(log2(pm(dat)[,,1]))
		}
		totMed2 <- totMed^2
		X <- cbind(totMed, totMed2, X)
		#X <- cbind(totMed, X)
	}
    affinity <- predict(lmfit, newdata=data.frame(X))
    return(affinity)
}  


getControlIndex <- function(dat, controlProbes=c("CONTROL_PROBES", "CONTROL_REGIONS"), 
			noCpGWindow=NULL, subject, onlyGood=FALSE, matrix=TRUE) {
    if (is.null(noCpGWindow)) {
        controlIndex <- which(getContainer(dat) %in% controlProbes)
    } else {
        if (class(subject)!="BSgenome")
            stop("You must supply a BSgenome object for subject if using the noCpGWindow option. E.g. Hsapiens from the BSgenome.Hsapiens.UCSC.hg18 package.\n")
        chr <- pmChr(dat)
        pos <- pmPosition(dat)
        cpgd <- cpgdensity(subject, chr=chr, pos=pos, 
			windowSize=noCpGWindow, showProgress=FALSE)
        controlIndex <- which(cpgd==0)
    }
    if (onlyGood==FALSE) {
        ret <- controlIndex
    } else {
        datPm <- log2(pm(dat)[controlIndex,,]) 			
        M <- datPm[,,"channel1"] - datPm[,,"channel2"] 
        A <- (datPm[,,"channel1"] + datPm[,,"channel2"])/2
        samps <- sampleNames(dat)
        if (matrix) {
            ret <- sapply (samps, function(samp) { 
                m <- M[, samp]
                a <- A[, samp]
                breakpoint <- optimize(checkBreakpoint, range(a), ms=m, as=a)$minimum
                a>breakpoint
            })
        } else {
            ret <- lapply (samps, function(samp) { 
                m <- M[, samp]
                a <- A[, samp]
                breakpoint <- optimize(checkBreakpoint, range(a), ms=m, as=a)$minimum
                keep <- a>breakpoint
                controlIndex[keep]
            })
            names(ret) <- samps
        }
    }
    return(ret)
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

cpgdensity <-function(subject, chr, pos, windowSize=500, sequence="CG", showProgress=FALSE) {
    idx <- split(1:length(chr), chr)
    s <- DNAString(sequence) 
    cpgdensity <- rep(NA, length(pos))
    for (curchr in (names(idx))) {
            if (showProgress) cat(".")
            if (curchr %in% names(subject)) {
                chrseq <- subject[[curchr]]
                curpos <- pos[idx[[curchr]]]
                v <- suppressWarnings(Views(chrseq, start=curpos-windowSize/2, end=curpos+windowSize/2))
                d <- suppressWarnings(DNAStringSet(v))
                numcpg <- vcountPattern(s, d, fixed=TRUE)
                cpgdensity[idx[[curchr]]] <- numcpg/windowSize       
            }
    }
    if (showProgress) cat("\n")
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

qcReport <- function(dat, file=NULL, utRange=c(30,100), enRange=c(8,12), numProbes=5e+5, cluster=NULL) {
    # Calculate summary quality scores
	n <- nrow(pm(dat))
	if (numProbes==-1) numProbes <- n
	if (n < numProbes) numProbes <- n
	idx <- as.integer(seq(1, n, length.out=numProbes))
    pmQual <- pmQuality(dat, idx=idx, cluster=cluster) 
    X <- getX(dat, "pm")[idx]
    Y <- getY(dat, "pm")[idx]
    imgs1 <- arrayImage(X,Y, pmQual, cluster=cluster)
 	#sd1 <- unlist(lapply(imgs1, function(x) sd(as.vector(x$z), na.rm=TRUE)))

	if (class(dat)=="OfflineTilingFeatureSet2") {
		pms <- pm(dat)
		c1 <- ff(as.vector(log2(pms[idx,,"channel1"])), dim=c(length(idx), ncol(pms)), 
			finalizer="close")
		colnames(c1) <- colnames(pms)	
		tmp <- arrayImage(X,Y, c1, cluster=cluster)
 		sd1 <- unlist(lapply(tmp, function(x) sd(as.vector(x$z), na.rm=TRUE)))

		c2 <- ff(as.vector(log2(pms[idx,,"channel2"])), dim=c(length(idx), ncol(pms)), 
			finalizer="close")
		colnames(c2) <- colnames(pms)	
		imgs2 <- arrayImage(X,Y, c2, cluster=cluster)
	 	sd2 <- unlist(lapply(imgs2, function(x) sd(as.vector(x$z), na.rm=TRUE)))
	} else {
		tmp <- arrayImage(X,Y, log2(pm(dat)[idx,,"channel1"]))
 		sd1 <- unlist(lapply(tmp, function(x) sd(as.vector(x$z), na.rm=TRUE)))
		imgs2 <- arrayImage(X,Y, log2(pm(dat)[idx,,"channel2"]))
	 	sd2 <- unlist(lapply(imgs2, function(x) sd(as.vector(x$z), na.rm=TRUE)))
	}
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

arrayImage <- function(x,y,z, view="2d", blockSize=50, cluster=NULL) {
    if (ncol(z)==1) z <- as.matrix(z)
    if (view=="col") {
		if(any(class(z)=="ff")) stop("arrayImage doesn't yet support view='col' for OfflineTilingFeatureSet2\n")
		tmp <- vec2array(x,y,z)
        ret <- apply(tmp, 3, rowMeans, na.rm=TRUE)
        colnames(ret) <- colnames(z)
    }
    if (view=="row") {
		if(any(class(z)=="ff")) stop("arrayImage doesn't yet support view='row' for OfflineTilingFeatureSet2\n")
		tmp <- vec2array(x,y,z)
        ret <- apply(aperm(tmp, c(2,1,3)), 3, rowMeans, na.rm=TRUE)  
        colnames(ret) <- colnames(z)
    }
    if (view=="2d") {
		nx <- max(x)/blockSize
		ny <- max(y)/blockSize
		d <- discretize.image(cbind(y,x), m=ny, n=nx)
		if (any(class(z)=="ff")) {
			cl <- charmCluster(cluster)	
	        ret <- parLapply(cl, 1:ncol(z), function(i, d) {
				Z <- z[,i]
				as.image(Z, ind=d$index, nx=d$m, ny=d$n, na.rm=TRUE)		
	        }, d)   			
			if (all(class(cluster)!="cluster")) stopCluster(cl)     
		} else {
			ret <- apply(z, 2, function(vec) {
				as.image(vec, ind=d$index, nx=d$m, ny=d$n, na.rm=TRUE)		
	        })
		}
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

plotDensityMat <- function(x, cols=rep(1:5, 2), 
        lwd=rep(1:4, each=5), main=NULL, ylab="Density", xlab="p", lab=NULL,rx=NULL,ry=NULL,legendPos="topright") {
    d <- apply(x, 2, density, na.rm=TRUE)
    if (is.null(rx)) rx <- range(sapply(d, function(i) i$x))
    if (is.null(ry)) ry <- range(sapply(d, function(i) i$y))
    plot(rx, ry, type = "n", xlab = xlab, ylab = ylab, main=main)
    sapply(1:length(d), function(i) lines(d[[i]], col=cols[i], lwd=lwd[i]))
    legend(legendPos, legend=lab, col=cols, lwd=lwd, cex=0.6)    
}

## fn(pm) where fn is the ECDF of bg
pmQuality <- function(dat, channel="channel1", verbose=FALSE, idx=NULL, cluster=NULL) {
	if (is.null(idx)) idx <- 1:nrow(pm(dat))
    Ngc <- countGC(dat, "pm")[idx]
    bgNgc <- countGC(dat, "bg")  
	if (class(dat)=="OfflineTilingFeatureSet2") {
		pms <- pm(dat)
		bgs <- bg(dat)
		cl <- charmCluster(cluster)
		pmq <- ff(dim=c(length(idx), ncol(pms)), vmode="double", finalizer="close")
		parSapply(cl, 1:ncol(pms), function(i, pms, bgs, Ngc, bgNgc) {
		    if (verbose) cat(".")
		    fn <- tapply(bgs[,i,channel], bgNgc, ecdf)
		    ret <- rep(NA, length(Ngc))
			for (ngc in unique(Ngc)) {
		        idx2 <- Ngc==ngc
		        closestIdx <- order(abs(as.numeric(names(fn))-ngc))[1]
		        bgngc <- as.character(names(fn)[closestIdx])
		        #pmq[idx,i] <- 100 * fn[[bgngc]](pms[idx,i,channel])   
				ret[idx2] <- 100 * fn[[bgngc]](pms[idx,i,channel][idx2])   
		    }
			pmq[,i] <- ret
			return(NULL)
	    }, pms, bgs, Ngc, bgNgc)
	   	if (all(class(cluster)!="cluster")) stopCluster(cl) 
	} else {
		pms <- pm(dat)[idx,,]
		bgs <- bg(dat)
		pmq <- sapply(1:ncol(pms), function(i) {
		    if (verbose) cat(".")
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
	}
	colnames(pmq) <- sampleNames(dat)
	return(pmq)
}

pmvsbg <- function(...) {
	cat("The pmvsbg function has been renamed pmQuality. Both names work for now but please update your code soon\n")
	pmQuality(...)
}

dynamicRange <- function(dat, prob=0.8) {
    Ngc <- countGC(dat, "pm")
    bgNgc <- countGC(dat, "bg")
    c1 <- log2(oligo::pm(dat)[,,"channel1"])
    c2 <- log2(oligo::bg(dat)[,,"channel2"])    
    sapply(sampleNames(dat), function(i) {
            cat(".")
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
	if (any(class(Z)=="ff")) {
		retval <- ff(NA, dim=c(maxX*maxY, ncol(Z)), vmode="double", finalizer="close")
		for (i in 1:ncol(Z)){
			retval[tmp,i] <- Z[o,i]
		}
		dim(retval) <- c(maxY, maxX, ncol(Z))
	} else {
		Z <- as.matrix(Z[o,])
	    d <- matrix(NA, nrow=maxX*maxY, ncol=ncol(Z))
	    d[tmp,] <- Z
	    retval <- array(d, dim=c(maxY, maxX, ncol(Z)))
	}
	retval
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

findClosestGene <- function (chrom, pos, genome = "hg17", position = "txStart") {
    ## This function modified from the version in the ACME package to allow for chromosomes
    ## not present in the annotation file
    if (!exists("refflat")) {
        refflat <<- list()
        refflat[[genome]] <<- getRefflat(genome)
    }
    else if (!match(genome, names(refflat))) {
        refflat[[genome]] <<- getRefflat(genome)
    }
    rf <- refflat[[genome]]
    chromsub <- rf$chrom == chrom
    if (sum(chromsub)>0) {
        diffdist <- rf[chromsub, position] - pos
        sub <- which(abs(diffdist) == min(abs(diffdist)))
        rf <- rf[chromsub, 1:9][sub, ]
        return(data.frame(rf, Distance = diffdist[sub]))
    } else {
        return(data.frame(geneName="", name="", chrom=0, txStart=0, strand=""))
    }
    
}

findClosestGeneTxt <- function (chrom, pos, genome = "txtFile", position = "txStart") {
    if (chrom != "GROUPUN") {
        if (!exists("ann")) {
            ann <- read.delim(genome)
        }
        rf <- ann
        chromsub <- rf$chrom == chrom
        if (sum(chromsub)>0) {
            diffdist <- rf[chromsub, position] - pos
            sub <- which(abs(diffdist) == min(abs(diffdist)))
            rf <- rf[chromsub,][sub, ]
            return(data.frame(rf, Distance = diffdist[sub]))
        } else {
            return(data.frame(geneName="", name="", chrom=0, txStart=0, strand=""))
        }
    } else {
        return(data.frame(geneName="", name="", chrom=0, txStart=0, strand=""))
    }
}


closestTSS <- function(chr, start, end, genome="hg18", annotationTextFile=FALSE) {
    m <- rowMeans(cbind(start,end))
    t(sapply(1:length(chr), function(i) {
        if (!annotationTextFile) {
            gene <- findClosestGene(chrom=chr[i], pos=m[i], genome = genome, position = "txStart")[1,]
        } else {
            gene <- findClosestGeneTxt(chrom=chr[i], pos=m[i], genome = genome, position = "txStart")[1,]
        }
        if (gene$txStart>=start[i] & gene$txStart<=end[i]) {
            distance <- 0
            relationToTSS <- "overlaps"
        } else if (gene$chrom==0) {
            distance=0
            signedDist=0
            relationToTSS=""
        } else {
            distance <- min(abs(start[i]-gene$txStart), abs(end[i]-gene$txStart))
            signedDist <- (m[i]-gene$txStart) * ifelse(gene$strand=="+", 1, -1)
            relationToTSS <- ifelse (signedDist<=0, "upstream", "downstream")
        }
        c(symbol=as.character(gene$geneName), refSeq=as.character(gene$name), 
            distanceToTSS=distance, relationToTSS=relationToTSS, TSS=gene$txStart)
    }))
}

annotateDMRs <-function(tab, organism="hg18", minCpG=2, minProbes=3, minPDiff=0, minArea=0, maxRows=-1, sortCol=NULL) {
	if (!is.null(sortCol)) {
		o <- order(tab[,sortCol], decreasing=TRUE)
		tab <- tab[o,]
	}
	idx <- which((tab$start-tab$end)>0)
	if (length(idx)>0) {
		cat("Removing row(s)", idx, "because start is after end\n")
		tab <- tab[-idx,]
	}
    annotationTextFile <- FALSE
    if (organism=="hg18") {
        if (!require(BSgenome) & !require(BSgenome.Hsapiens.UCSC.hg18))
            stop("Please install the 'BSgenome' and 'BSgenome.Hsapiens.UCSC.hg18' packages")
        if (!require(org.Hs.eg.db))
            stop("Please install the 'org.Hs.eg.db' package")
        subject <- Hsapiens
        genome="hg18"
        entrez2NameMap <- org.Hs.egGENENAME 
        refseq2EntrezMap <- org.Hs.egREFSEQ2EG
    } else if (organism=="hg19") {
	        if (!require(BSgenome) & !require(BSgenome.Hsapiens.UCSC.hg19))
	            stop("Please install the 'BSgenome' and 'BSgenome.Hsapiens.UCSC.hg19' packages")
	        if (!require(org.Hs.eg.db))
	            stop("Please install the 'org.Hs.eg.db' package")
	        subject <- Hsapiens
	        genome="hg19"
	        entrez2NameMap <- org.Hs.egGENENAME 
	        refseq2EntrezMap <- org.Hs.egREFSEQ2EG
	 } else if (organism=="mouse") {
        if (!require(BSgenome) & !require(BSgenome.Mmusculus.UCSC.mm8))
            stop("Please install the 'BSgenome' and 'BSgenome.Mmusculus.UCSC.mm8' packages")
        if (!require(org.Mm.eg.db))
            stop("Please install the 'org.Mm.eg.db' package")
        subject <- Mmusculus
        genome="mm8"
        entrez2NameMap <- org.Mm.egGENENAME 
        refseq2EntrezMap <- org.Mm.egREFSEQ2EG
    } else if (organism=="bee") {
        genome <- "/home/bst/other/maryee/projects/methylation/data/feinberg/ApiMel4Genes.txt"
        subject <- NA
        entrez2NameMap <- NA
        refseq2EntrezMap <- NA
        annotationTextFile <- TRUE
    } else {
        subject <- organism$subject
        genome <- organism$genome
        entrez2NameMap <- organism$entrez2NameMap
        refseq2EntrezMap <- organism$refseq2EntrezMap
    }
	tab$chr <- as.character(tab$chr)
    tab$nprobes <- tab$indexEnd-tab$indexStart+1
    tab$length <- tab$end-tab$start+1
	
    if("p1" %in% colnames(tab)) tab$absDiff <- abs(tab$p2-tab$p1)
    if("M1" %in% colnames(tab)) tab$absDiff <- abs(tab$M2-tab$M1)
    if (class(subject)=="BSgenome") {
        tab$numcpg <- countSeq(subject, chr=tab$chr, start=tab$start, end=tab$end, seq="CG")
		if (minPDiff>0) {
        	tab <- subset(tab, numcpg>=minCpG & nprobes>=minProbes & absDiff>=minPDiff & area>=minArea) ## CHECK THE minArea !!!
		} else {
			tab <- subset(tab, numcpg>=minCpG & nprobes>=minProbes & area>=minArea) 
		}
    }
	if(nrow(tab)==0) return(tab)
    if (maxRows>0) {
        n <- min(maxRows, nrow(tab))
        tab <- tab[1:n,]
    }
    tab <- cbind(tab, closestTSS(tab$chr, tab$start, tab$end, genome=genome, annotationTextFile=annotationTextFile))
    tab$refSeq <- as.character(tab$refSeq)
    tab$TSS <- as.numeric(as.character(tab$TSS))
    tab$distanceToTSS <- as.numeric(as.character(tab$distanceToTSS))
    if (class(subject)=="BSgenome") {
        # Add Entrez IDs
        mapped <- mappedkeys(refseq2EntrezMap) 
        #xx <- as.list(refseq2EntrezMap[mapped]) 
        xx <- suppressWarnings(as.character(refseq2EntrezMap[mapped]))
        entrez <- xx[as.character(tab$refSeq)] #unlist is needed after as.character
        tab$entrez <- entrez
        #tab$entrez <- NA # these 3 lines only needed after as.list(refseq2EntrezMap[mapped]) 
        #idx <- match(tab$refSeq, names(entrez))
        #tab$entrez <- entrez[idx]
        # Add gene names
        mapped_genes <- mappedkeys(entrez2NameMap) 
        #xx <- as.list(entrez2NameMap[mapped_genes]) 
        xx <- as.character(entrez2NameMap[mapped_genes]) 
        geneName <- xx[as.character(tab$entrez)]  #unlist is needed after as.list
        tab$geneName <- geneName
        #tab$geneName <- NA # these 3 lines only needed after as.list(entrez2NameMap[mapped_genes])  
        #idx <- match(tab$entrez, names(geneName))
        #tab$geneName <- geneName[idx]
        # Convert to character
    } else {
        tab$entrez <- NA
        tab$geneName <- NA
    }
	tab$symbol <- as.character(tab$symbol)
    tab$relationToTSS <- as.character(tab$relationToTSS)
    return(tab)
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

spatialAdjust <- function(dat, cluster=NULL, blockSize=50, theta=1) {
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
	#clusterExport(cl, "spatialAdjustVec") 
	if (class(dat)=="OfflineTilingFeatureSet2") {
		cl <- charmCluster(cluster)
		tmp <- parLapply(cl, 1:NCOL, function(i, pms, bgs, x, y, bgX, bgY, nx, ny, d, dBg) {	
			tot <- log2(pms[,i, "channel1"])
		    bgTot <- log2(bgs[,i, "channel1"])
			sm <- spatialAdjustVec(tot, d)
			pms[,i, "channel1"] <- 2^sm$zAdj
			smBg <- spatialAdjustVec(bgTot, dBg, sm$ims)
			#bgs[,i, "channel1"] <- 2^bgTot
			bgs[,i, "channel1"] <- 2^smBg$zAdj
		    
			en <- log2(pms[,i, "channel2"])	
		    bgEn <- log2(bgs[,i, "channel2"] )	
			sm <- spatialAdjustVec(en, d)
			pms[,i, "channel2"] <- 2^sm$zAdj
			smBg <- spatialAdjustVec(bgEn, dBg, sm$ims)
			#bgs[,i, "channel2"] <- 2^bgTot
			bgs[,i, "channel2"] <- 2^smBg$zAdj
			return(NULL) ## the ff object is already updated on disk	
	    }, pms, bgs, x, y, bgX, bgY, nx, ny, d, dBg)
		if (all(class(cluster)!="cluster")) stopCluster(cl)
	} else {
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
	}    		
    return(dat)
}




spatialAdjust.poly <- function(dat, demedian=FALSE, cluster=NULL) {
    if (demedian) dat <- demedian(dat)
    X <- getX(dat, "pm")
    Y <- getY(dat, "pm")
    bgX <- getX(dat, "bg")  
    bgY <- getY(dat, "bg")     
    NCOL=dim(dat)["Samples"]
	pms <- pm(dat)
	bgs <- bg(dat)
	cl <- charmCluster(cluster)
	tmp <- parLapply(cl, 1:NCOL, function(i, pms, bgs, X, Y, bgX, bgY) {	
		tot <- log2(pms[,i, "channel1"])
	    en <- log2(pms[,i, "channel2"])
	    bgTot=log2(bgs[,i, "channel1"] )
	    bgEn=log2(bgs[,i, "channel2"] )
	
        # fit surface to enriched probes
        Z     <- en - median(en)
        fit   <- surf.ls(6, X, Y, Z)
        #tmp <- trmat(fit, min(X), max(X), min(Y), max(Y), 200)
        # adjust enriched probes
        Zhat  <- predict.trls(fit, X, Y)
        en <- en - Zhat
        # adjust enriched background probes
        Zhat  <- predict.trls(fit, bgX, bgY)
        bgEn <- bgEn - Zhat
		rm(fit); gc()
        # fit surface to total probes
        Z     <- tot - median(tot)
        fit   <- surf.ls(6, X, Y, Z)
        # adjust total probes
        Zhat  <- predict.trls(fit, X, Y)
        tot <- tot - Zhat
        # adjust total background probes
        Zhat  <- predict.trls(fit, bgX, bgY)
        bgTot <- bgTot - Zhat    
		rm(fit); gc()
		   		
		pms[,i, "channel1"] <- 2^tot
		pms[,i, "channel2"] <- 2^en
		bgs[,i, "channel1"] <- 2^bgTot
		bgs[,i, "channel2"] <- 2^bgEn
		
		if (any(class(pms)=="ff")) {
			return(NULL) ## the ff object is already updated on disk	
		} else {
			return(list(pm=pms[,i,], bg=bgs[,i,]))
		}
		
    }, pms, bgs, X, Y, bgX, bgY)
	if (all(class(cluster)!="cluster")) stopCluster(cl)
	if (all(class(pms)!="ff")) {
		for (i in 1:length(tmp)) {
			pms[,i,] <- tmp[[i]][["pm"]]
			bgs[,i,] <- tmp[[i]][["bg"]]
		}
		pm(dat) <- pms
		bg(dat) <- bgs
	}
    return(dat)
}


debatch <- function(dat, mod, B=5, sampSize=100000, seed=NULL) {
    lr <- dat$tot - dat$en
    mod0 <- cbind(rep(1,nrow(dat$pd)))
    if (is.vector(mod) | is.factor(mod))
        mod <- model.matrix(~as.factor(mod))
    #if (nrow(lr)>sampSize) {
    #    s <- sample(1:nrow(lr), sampSize)
    #} else {
    #    s <- 1:nrow(lr) 
    #}
    #sva <- irwsva(dat=lr[s,], mod=mod, mod0=mod0, B=B, seed=seed)
    sva <- irwsva(dat=lr, mod=mod, mod0=mod0, B=B, seed=seed)
    #cat(sva$n.sv, " surrogate batch variables identified\n")
    if (sva$n.sv>0) {
        Xs <- sva$sv
        colnames(Xs) <- paste("SV", 1:ncol(Xs), sep="")
        X <- cbind(mod, Xs)
        B <- lm(t(lr) ~ -1+X)$coefficients
        Bs <- B[-(1:ncol(mod)),]
        #lr2 <- lr - t(Xs %*% Bs) 
        adj <- t(Xs %*% Bs)
        dat$en <- dat$en + adj/2
        dat$tot <- dat$tot - adj/2
        dat$n.sv <- sva$n.sv
        dat$Xs <- Xs
    } else {
        dat$n.sv <- 0
    }
    return(dat)
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
bgAdjustBgp <- function (dat, cluster=NULL) {
	pms <- pm(dat)
	bgs <- bg(dat)
    Ngc <- countGC(dat, "pm")
    bgNgc <- countGC(dat, "bg")
    if (class(dat)=="OfflineTilingFeatureSet2") {
		cl <- charmCluster(cluster)
		parLapply(cl, 1:dim(dat)["Samples"], function(samp, pms, bgs, Ngc, bgNgc) {
	    	for (chan in 1:2) {
	            param <- bgParametersBgp(pms[,samp,chan], bgs[,samp,chan], Ngc, bgNgc, n.pts=2^14)
	            b <- param$sigma
	            pms[,samp,chan] <- pms[,samp,chan] - param$mu.gc[Ngc] - param$alpha * b^2
	            pms[,samp,chan] <- pms[,samp,chan] + b * ((1/sqrt(2 * pi)) * exp((-1/2) * ((pms[,samp,chan]/b)^2)))/pnorm(pms[,samp,chan]/b)
	            bgs[,samp,chan] <- bgs[,samp,chan] - param$mu.gc[bgNgc] - param$alpha * b^2
	            bgs[,samp,chan] <- bgs[,samp,chan] + b * ((1/sqrt(2 * pi)) * exp((-1/2) * ((bgs[,samp,chan]/b)^2)))/pnorm(bgs[,samp,chan]/b)
			}
			return(NULL) ## the ff object is already updated on disk	
	    }, pms, bgs, Ngc, bgNgc)
		if (all(class(cluster)!="cluster")) stopCluster(cl) 
	} else {
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
    	pm(dat) <- pms
	    bg(dat) <- bgs
	}
    return(dat)
}



##Assume the signal is S=X+Y where X~Normal(0, s2) and Y~Exp(alpha)
## Calculate E(Y|S)
logmethParameters <- function (pm, ngc, n.pts = 2^14) 
{
    max.density <- function(x, n.pts) {
        aux <- density(x, kernel = "epanechnikov", n = n.pts, 
            na.rm = TRUE)
        aux$x[order(-aux$y)[1]]
    }
    pmbg <- 0
    bg.data <- pm[pm < pmbg]
    #bgsd <- mad(bg.data, center=0, constant=1)
    #bgsd <- mad(bg.data, center=0)
    idx <- pm < pmbg
    bgsd <- median(tapply(pm[idx], ngc[idx], mad, center=0))
    sig.data <- pm[pm > pmbg]
    sig.data <- sig.data - pmbg
    #cat("Debug logmethParameters(): sig.data =", 
    #    100*round(length(sig.data)/length(pm), 2), "%\n")
    expmean <- max.density(sig.data, n.pts)
    alpha <- 1/expmean
    mubg <- pmbg
    list(alpha = alpha, mu = mubg, sigma = bgsd)
}


methPercent <- function(m, ngc, commonParams=TRUE, cluster=NULL) {
	param <- t(sapply(1:ncol(m), 
		function(i) logmethParameters(m[,i], ngc)))
	alpha <- unlist(param[,"alpha"])
	sigma <- unlist(param[,"sigma"])
	if(commonParams) {
		alpha[] <- median(alpha)
		sigma[] <- median(sigma)
	}

	if (any(class(m)=="ff")) {
		cl <- charmCluster(cluster)
		tmp <- ff(dim=dim(m), vmode="double", finalizer="close")
		parSapply(cl, 1:ncol(m), function(i, m, ngc, alpha, sigma, tmp) {	
			x <- m[,i]
	        #param <- logmethParameters(x, ngc, n.pts=2^14)
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
			#tmp[,i] <- 1-2^(-postM)
			tmp[,i] <- 1-2^(-postM)
			return(NULL)
	    }, m, ngc, alpha, sigma, tmp)
		ret <- tmp
		colnames(ret) <- colnames(m)
	   	if (all(class(cluster)!="cluster")) stopCluster(cl) 
	} else {
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
	}
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
    if(verbose) if(i%%1000==0) cat(i,",")
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
  require(genefilter)
  
  gIndex=split(seq(along=groups),groups)
  gIndex=gIndex[which(names(gIndex)%in%compare)]
  ng=length(gIndex)
  
  lm=matrix(0,nrow(l),ng)
  ls=matrix(0,nrow(l),ng)
  ns=sapply(gIndex,length)
  
  if(verbose) cat("Computing group medians and SDs for",ng,"groups:")
  for(i in seq(along=gIndex)){
    if(verbose) cat("\n",i)
    Index=gIndex[[i]]
    if(length(Index)>1){
        lm[,i]=rowMedians(l[,Index,drop=FALSE])
        ls[,i]=rowMads(l[,Index,drop=FALSE],center=lm[,i])
    } else{
        cat(paste(" ",names(gIndex)[i],"has only 1 array!"))
        lm[,i]=l[,Index,drop=FALSE]
        ls[,i]=NA
    }
  }
  colnames(lm)=names(gIndex)
  colnames(ls)=names(gIndex)

  nums  <- match(compare,colnames(lm))
  COMPS <- matrix(nums,ncol=2,byrow=TRUE)
  
  if(verbose) cat("\nDone.\n")
  return(list(lm=lm,ls=ls,ns=ns,COMPS=COMPS))
}

get.tt <- function(lm,ls,ns,filter,Indexes,COMPS,ws,verbose){
##this function assumes genome position of mat[Indexes[[1]] are ordered:
  if(is.null(filter)){
    Tukey = function(x) pmax(1 - x^2,0)^2
    filter= Tukey(seq(-ws,ws)/(ws+1));filter=filter/sum(filter)
  }
  if(!isTRUE(all.equal(sum(filter),1))) stop("filter must sum to 1.")
  
  if(verbose) cat("Smoothing")
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
  for(i in seq(along=Indexes)){
    if(verbose) if(i%%1000==0) cat(".")
    Index=Indexes[[i]]
    for(r in 1:nrow(COMPS)){
        j=COMPS[r,1]
        k=COMPS[r,2]
        sdm[Index,r]=myfilter(dm[Index,r],filter)
        if(ns[j]==1|ns[k]==1){ svr[Index,r] = 1} else{ svr[Index,r] = myfilterse(vr[Index,r],filter) }
    }
  }
  for(r in 1:nrow(COMPS)) tt[,r] = sdm[,r]/svr[,r]
  if(verbose) cat("Done.\n")
  return(tt)
}


dmrFinder <- function(eset=NULL,groups,p=NULL,l=NULL,chr=NULL,pos=NULL,pns=NULL,
					  sdBins=NULL, controlIndex=NULL,Indexes=NULL, 
					  filter=NULL,package=NULL,ws=7,
                      verbose=TRUE,compare="all",
					  bgSubtract=TRUE,
		      		  withinSampleNorm="loess", binSize=500, numSegments=3,
					  betweenSampleNorm="quantile",
                      minQCScore=NULL, 
					  controlProbes=c("CONTROL_PROBES", "CONTROL_REGIONS"), 
					  cluster=NULL,
					  cutoff=0.995,...){
  groups = as.character(groups)
  if(identical(compare,"all")) compare=comp(groups)
  if(length(compare)%%2!=0) stop("compare must have an even number of elements.")

  args=list(filter=filter, ws=ws, betweenSampleNorm=betweenSampleNorm, 
	withinSampleNorm=withinSampleNorm, minQCScore=minQCScore,
            controlProbes=controlProbes, cutoff=cutoff)

  # dmrFinder must be given either eset or p/l,chr,pos,pns, and controlIndex
  # If eset is supplied then all the latter will be taken from it (with any
  # that were given as arguments being ignored) (except p).
  # l=logit(p). Indexes=split(seq(along=pns),pns).
  if(is.null(eset)){
      if (is.null("p") & is.null("l")) stop("p or l must be supplied.")
      args = c("chr","pos","pns","controlIndex")
      nulls = sapply(args,function(x) is.null(get(x)))
      if(any(nulls))
        stop(paste("The following arguments are missing:", paste(args[nulls], collapse=", ")))
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
	  pdInfo=get(eset)
	  class(pdInfo)="TilingFeatureSet" # Trick oligo so that pmChr, pmPosition, probeNames work
	  chr=pmChr(pdInfo)
	  pos=pmPosition(pdInfo)
	  index=which(rowSums(is.na(p))==0)
      index=index[order(chr[index],pos[index])]
      chr=chr[index]
      pos=pos[index]
      pns=probeNames(pdInfo)[index]
      p=p[index,]
	  controlIndex <- which(getContainer(pdInfo) %in% controlProbes)
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
      stopifnot(length(groups)==length(eset))
      if(is.null(p) & is.null(l)){
		p <- methp(eset, bgSubtract=bgSubtract,
					withinSampleNorm=withinSampleNorm, binSize=binSize, numSegments=numSegments,
					betweenSampleNorm=betweenSampleNorm,
					controlProbes=controlProbes, 
					cluster=cluster, verbose=TRUE)[,]
	  }
      chr=pmChr(eset)
      pos=pmPosition(eset)
      index=which(rowSums(is.na(p))==0)
      index=index[order(chr[index],pos[index])]
      chr=chr[index]
      pos=pos[index]
      pns=probeNames(eset)[index]
      p=p[index,]

      controlIndex=getControlIndex(eset, controlProbes=controlProbes)
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
  Indexes=split(seq(along=pns),pns)

  tog = get.tog(l=l,groups=groups,compare=compare,verbose=verbose)
  lm=tog$lm
  ls=tog$ls
  ns=tog$ns
  COMPS=tog$COMPS
  tt = get.tt(lm=lm,ls=ls,ns=ns,COMPS=COMPS,Indexes=Indexes,filter=filter,ws=ws,verbose=verbose)

  res=vector("list",ncol(tt))
  names(res)=colnames(tt)
  if(verbose) cat("Finding DMRs for each pairwise comparison.")
  for(r in 1:nrow(COMPS)){
      j = COMPS[r,1]
      k = COMPS[r,2]
      if(verbose) cat("\n",colnames(tt)[r])
      DF=ifelse(ns[j]==1 & ns[k]==1, 1, ns[j]+ns[k]-2)
	
	  if (length(sdBins)==0) {
	      K=mad(tt[,r], na.rm=TRUE)*qt(cutoff,DF)	
	  }	else {
		  s <- tapply(tt[,r], sdBins, mad, na.rm=TRUE)
		  K=s[sdBins]*qt(cutoff,DF)	
	  }
      LAST=0
      segmentation=vector("numeric",nrow(tt))
      type=vector("numeric",nrow(tt))
      for(i in seq(along=Indexes)){
        if(verbose) if(i%%1000==0) cat(".")
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
      res[[r]]=data.frame(
           chr=tapply(chr[Index],segmentation[Index],function(x) x[1]),
           start=tapply(pos[Index],segmentation[Index],min),
           end=tapply(pos[Index],segmentation[Index],max),
           p1=NA,
           p2=NA,
           regionName=tapply(pns[Index],segmentation[Index],function(x) x[1]),
           indexStart=tapply(Index,segmentation[Index],min),
           indexEnd=tapply(Index,segmentation[Index],max))
	  if (is.null(p)) { #  We return log-ratios
		  colnames(res[[r]]) <- sub("p1", "m1", colnames(res[[r]]))
		  colnames(res[[r]]) <- sub("p2", "m2", colnames(res[[r]]))
	      res[[r]]$m1=tapply(lm[Index,j],segmentation[Index],mean)
          res[[r]]$m2=tapply(lm[Index,k],segmentation[Index],mean)
    	  length=res[[r]]$indexEnd-res[[r]]$indexStart+1
		  area=abs(res[[r]]$m2-res[[r]]$m1)*length		
		  res[[r]]$area=area
	  } else { # We return percentages
	      res[[r]]$p1=tapply(1/(1+exp(-lm[Index,j])),segmentation[Index],mean)
          res[[r]]$p2=tapply(1/(1+exp(-lm[Index,k])),segmentation[Index],mean)
    	  length=res[[r]]$indexEnd-res[[r]]$indexStart+1
		  area=abs(res[[r]]$p2-res[[r]]$p1)*length
		  res[[r]]$area=area
	  }
      res[[r]]=res[[r]][order(-area),]
  }
  if(verbose) cat("\nDone\n")
  return(list(tabs=res,p=p,m=l,chr=chr,pos=pos,pns=pns,
		      index=index,controlIndex=controlIndex,
              gm=lm,groups=groups,args=args,cutoff=cutoff,
			  filter=filter,ws=ws,package=package))
}

dmrFdr <- function(dmr, compare=1, numPerms=1000, seed=NULL, verbose=TRUE) {
	if (length(compare)!=1) stop("You must choose one comparison at a time when calculating FDRs. Please set dmr to be one of: ", 
	paste(names(dmr$tabs), collapse=", "), "\n")
	if (is.numeric(compare)) compare <- names(dmr$tabs)[compare]
	cat("Calculating q-values for DMRs between", compare, "\n")
	# Get probe order from TilingFeatureSet object
	pdInfo=get(dmr$package)
	class(pdInfo)="TilingFeatureSet" # Trick oligo so that pmChr, pmPosition work
	chr=pmChr(pdInfo)
	pos=pmPosition(pdInfo)
	o <- order(chr, pos)
	p <- matrix(NA, nrow=nrow(dmr$p), ncol=ncol(dmr$p))
	p[o,] <- dmr$p

	keep <- dmr$groups %in% unlist(strsplit(compare, "-"))
	p <- p[,keep]
	n <- sum(keep)
	n1 <- sum(dmr$groups==unlist(strsplit(compare, "-"))[1])
	maxPerms <- choose(n, n1)
	if (numPerms=="all") numPerms <- maxPerms
	if (numPerms>maxPerms) {
		cat("Given the sample sizes in the two groups the maximum number of permutations is ", maxPerms, ".\n")
		numPerms <- maxPerms
	} 

	## Reshuffled group label DMRs
	if (!is.null(seed)) set.seed(seed)
	s <- sample(1:maxPerms, numPerms)
	grp1 <- combinations(n,n1)[s,]

	if (verbose) cat("Finding permuted data DMRs. Estimating time remaining\n")
	areas <- lapply(1:numPerms, function(i) {
		groups <- rep("grp2", n)
		groups[grp1[i,]] <- "grp1"
		st <- system.time(dmrPerm <- dmrFinder(dmr$package, p=p, 
			groups=groups, cutoff=dmr$cutoff, 
			filter=dmr$filter, ws=dmr$ws, verbose=FALSE))[3]
		if (verbose & (i %in% round(seq(1, numPerms, length.out=10)))) {
			cat(i, "/", numPerms, " (", prettyTime((numPerms-i)*st), 
				" remaining)\n", sep="")
		}
		dmrPerm$tabs[[1]]$area
	})

	nullDist <- unlist(areas)
	fn <- ecdf(nullDist)
	pval <- 1-fn(dmr$tabs[[compare]]$area)
	pi0<-pi0.est(pval)$p0
	qval<-qvalue.cal(pval, pi0)
	dmr$tabs[[compare]] <- cbind(dmr$tabs[[compare]], pval, qval)
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

.onAttach <- function(libname, pkgname) {
 require(snow)
 message("Welcome to charm version ", packageDescription("charm", field="Version"))
}

#############################################
### Old bits of code not currently in use ###
#############################################
if (FALSE) { 

normalizeLoess.onlyGoodCtrl <- function(dat, controlIndex=NULL, onlyGoodCtrl=FALSE, controlProbes=c("CONTROL_PROBES", "CONTROL_REGIONS"), by="tot", iterations=4, span=0.3) {
		if (is.null(controlIndex)) {
	    	controlIndex <- getControlIndex(dat, controlProbes=controlProbes)
		}
		if (onlyGoodCtrl){
			goodCtrl <- getControlIndex(dat, controlProbes = controlProbes, onlyGood=TRUE)		
		} else {
			goodCtrl <- matrix(TRUE, nrow=length(controlIndex), ncol=dim(dat)["Samples"])
		}
		datPm <- log2(pm(dat))
	    M <- datPm[, , "channel1"] - datPm[, , "channel2"]
		for (i in 1:ncol(M)) {
			cat(i,  length(controlIndex[goodCtrl[,i]]), "probes\n")
			y <- M[,i]
			if (by=="tot") {
				x <- datPm[,i,"channel1"]			
			} else if (by=="A") {
				x <- (datPm[, i, "channel1"] + datPm[, i, "channel2"])/2
			}
	 	    fit <- loess(y ~ x, subset = controlIndex[goodCtrl[,i]],
	               na.action = na.exclude, degree = 1, surface = "direct",
	               family = "symmetric", trace.hat = "approximate",
	               iterations = iterations, span=span)
			datPm[, i, "channel2"] <- datPm[, i, "channel2"] + predict(fit, newdata=x)
		}
		pm(dat) <- 2^datPm
		return(dat)
	}

normalizeWithinSamples.20091111 <- function (dat, strata = "affinity", method = "median", affinityMode="ctrlMedianM", controlProbes = "CONTROL_REGIONS", bins=NULL,
	    controlIndex = NULL, binSize=250, onlyGoodCtrl=FALSE,
		numSegments=3, window=NULL, subject=NULL, verbose=FALSE) 
	{
	    if (strata=="gc") {
			if (class(dat)=="OfflineTilingFeatureSet2") {
				stop("normalizeWithinSamples does not yet support the strata='gc' option for OfflineTilingFeatureSet2 objects. Please use the strata='affinity' option.")
			}
	        mAdj <- diffAmpEstGC(dat, method = method, controlProbes = controlProbes)
	        dat <- diffAmpAdjustGC(dat, mAdj)
	    } else if (strata=="affinity") {   
	        if (is.null(controlIndex)) {
	            controlIndex <- getControlIndex(dat, controlProbes=controlProbes)
	            if (onlyGoodCtrl) {
					goodCtrl <- getControlIndex(dat, controlProbes = controlProbes, onlyGood=TRUE)
				}
	        }
	        if (affinityMode=="ctrlMedianM") {
				if (is.null(bins)) {
					bins <- affinityBin(dat, controlIndex=controlIndex, 
						binSize=binSize, verbose=verbose)
				}
				#datPm <- log2(pm(dat))
				pms <- pm(dat)
				for (i in 1:ncol(pms)) {	
					mAdj <- bins$binMed[bins$bin,i]
	                pms[, i, "channel2"] <- 2^(log2(pms[, i, "channel2"]) + mAdj)
	            }        
				pm(dat) <- pms
	        } else if (affinityMode=="predicted") {
				if (class(dat)=="OfflineTilingFeatureSet2") {
					stop("normalizeWithinSamples does not yet support the affinityMode='predicted' option for OfflineTilingFeatureSet2 objects. Please use the affinityMode='ctrlMedianM' option.")
				}
				datPm <- log2(pm(dat))
		        M <- datPm[, , "channel1"] - datPm[, , "channel2"]
				if (!is.matrix(M)) M <- as.matrix(M)
	            for (i in 1:ncol(M)) {
	                if (onlyGoodCtrl) {
	                    affinity <- getAffinity(dat[,i], channel="M", idx=controlIndex, exclude=!goodCtrl[,i], numSegments=numSegments, window=window, subject=subject, verbose=verbose)                
	                } else {
	                    affinity <- getAffinity(dat[,i], channel="M", idx=controlIndex,
							numSegments=numSegments, window=window, subject=subject, verbose=verbose)                
	                }
	                #affinity <- getAffinity(dat[,i], channel="M", idx=controlIndex)                
	                datPm[, i, "channel2"] <- datPm[, i, "channel2"] + affinity
	            }
	        	pm(dat) <- 2^datPm       
			}

	    } else if (strata=="none") {
			if (class(dat)=="OfflineTilingFeatureSet2") {
				stop("normalizeWithinSamples does not yet support the strata='none' option for OfflineTilingFeatureSet2 objects. Please use the strata='affinity' option.")
			}
	        if (is.null(controlIndex)) {
	            controlIndex <- getControlIndex(dat, controlProbes = controlProbes)
	            if (onlyGoodCtrl) goodCtrl <- getControlIndex(dat, controlProbes = controlProbes, onlyGood=TRUE)
	        }
	        datPm <- log2(pm(dat))
	        M <- datPm[, , "channel1"] - datPm[, , "channel2"]
	        mCtrl <- M[controlIndex, ]
	        if (onlyGoodCtrl) mCtrl[!goodCtrl] <- NA
	        if (method=="max.density") {
	            mAdj <- apply(mCtrl, 2, "max.density")
	        } else {
	            mAdj <- apply(mCtrl, 2, method, na.rm=TRUE)
	        }    
	        datPm[, , "channel2"] <- sweep(datPm[, , "channel2"], 
	            2, mAdj, FUN = "+")
	        pm(dat) <- 2^datPm
	    } else {
			stop("Invalid strata option. Please choose 'affinity', 'gc' or 'none'\n")
		}
	    return(dat)	
}


			makeSeqXGC <- function(dat, type="pm", idx=NULL, numSegments=3, segmentSize=NULL, MAT=FALSE, addInteractions=FALSE, addSquared=FALSE, window=NULL, subject=NULL) {
			    if (MAT) {
			        X1 <- makeX(dat=dat, idx=idx, type=type, segmentSize=1, letters=c("A", "C", "G"))
			        X2 <- makeX(dat=dat, idx=idx, type=type, segmentSize=25, letters=c("A", "C", "G", "T"))
			        colnames(X2) <- c("A.sq", "C.sq", "G.sq", "T.sq")
			        X <- cbind(X1, X2^2)        
			    } else {
			        X <- makeX(dat=dat, idx=idx, type=type, numSegments=numSegments, 
						segmentSize=segmentSize, letters="GC")
					if(!is.null(window)) {
						Xw <- makeX(dat=dat, idx=idx, window=window, numSegments=1, subject=subject,  
							letters=c("GC"))
						#colnames(Xw) <- paste("W", colnames(Xw), sep="")
						colnames(Xw) <- "Wgc"
						X <- cbind(X, Xw)
					}
					if(addInteractions) {
						X <- addInteractions(X)
					} else if (addSquared) {
						X1 <- makeX(dat=dat, idx=idx, type=type, numSegments=1, letters=c("GC"))
						Xsq <- X1^2
						colnames(Xsq) <- paste(colnames(X1), "sq", sep=".")
						X <- cbind(X, Xsq)
						if (!is.null(window)){
							Xwsq <- Xw^2
							colnames(Xwsq) <- paste(colnames(Xw), "sq", sep=".")
							X <- cbind(X, Xwsq)
						}
					}
			    }
			    return(as.data.frame(X))
			}    

		getAffinitySmoothTot <- function(dat, type="pm", channel=2, idx=NULL, exclude=NULL, segmentSize=NULL, numSegments=3, window=NULL, subject=NULL, MAT=FALSE) {
			    if(is.null(segmentSize)) {
			        seqs = as.character(pmSequence(dat))
			        d <- DNAStringSet(seqs)
			        avgLength <-  round(mean(width(d)))
			        segmentSize <- floor(avgLength/numSegments)
			    }
			    if (type=="pm") {
			        if (channel=="M") {
			            y <- getM(dat)[pmindex(dat), ]
			        } else {
			            y <- log2(pm(dat))[,,channel]
			        }
			        if (is.vector(y)) y <- as.matrix(y)
			        if (!is.null(idx)) y<-y[idx,]
			        if (is.vector(y)) y <- as.matrix(y)
			    } else if (type=="bg") {
			        if (channel=="M") {
			            y <- getM(dat)[bgindex(dat), ]
			        } else {
			            y <- log2(bg(dat))[,,channel]
			        }
			        if (is.vector(y)) y <- as.matrix(y)
			    }
			    if (!is.null(exclude)) {
			        exclude <- as.matrix(exclude)
			        y[exclude] <- NA
			    }
				chr <- pmChr(dat)
				pos <- pmPosition(dat)
				o <- order(chr[idx], pos[idx])
			    yMed <- rowMedians(y, na.rm=TRUE)
				totMed <-  rowMedians(log2(pm(dat)[idx,,1]))
				X <- makeSeqX(dat, type=type, idx=idx, segmentSize=segmentSize, window=window, subject=subject)
				f <- rep(1/50, 50)
				yMedS <- myfilter(yMed[o], f)
				totMedS <- myfilter(totMed[o], f)
				X <- X[o,]
				rhs <- paste(c("totMedS", colnames(X)), collapse="+")
				formul <- formula(paste("yMedS~", rhs))
				lmfit.base <- lm(formul, data=X)
				lmfit <- step(lmfit.base, ~ .^2, trace=0)
				cat("  Built smoothed control probe affinity model with", length(coef(lmfit)), "coefficients. R2 =", round(summary(lmfit)$r.squared,3), "\n")
			    X <- makeSeqX(dat, segmentSize=segmentSize, window=window, subject=subject)
				totMedS <-  rowMedians(log2(pm(dat)[,,1])) # need to call it totMedS to match name in model
				X <- cbind(totMedS, X)
			    affinity <- predict(lmfit, newdata=data.frame(X))
			    return(affinity)
			}  


		normalizeWithinSamples.old <- function(dat, method="median", controlProbes="CONTROL_REGIONS", controlIndex=NULL, useGC=FALSE) {
		    if (useGC) {
		        mAdj <- diffAmpEstGC(dat, method=method, controlProbes=controlProbes) 
		        dat <- diffAmpAdjustGC(dat, mAdj)
		    } else {
		        if (is.null(controlIndex)) controlIndex <- getControlIndex(dat, controlProbes=controlProbes)
		        datPm <- log2(pm(dat)) 
		        M <- datPm[,,"channel1"] - datPm[,,"channel2"]    
		        mAdj <- apply(M[controlIndex,], 2, method)
		        datPm[,,"channel2"] <- sweep(datPm[,,"channel2"], 2, mAdj, FUN="+")
		        pm(dat) <- 2^datPm
		    }
		    return(dat)
		}
		getMap.old <- function(x, m, df=5) {
		    nas <- is.na(x) | is.na(m)
		    x <- x[!nas]
		    m <- m[!nas]
		    d <- x-m
		    #l <- loess(d~x, surface = "direct")
		    #return(function(x1) predict(l, x1))
		    lmFit=lm(d~ns(x,df=df))
		    return(function(x1) x1 - predict(lmFit, data.frame(x=x1)))
		}


		normalizeBetweenSamples2 <- function(dat, controlProbes="CONTROL_REGIONS", onlyGood=FALSE) {
		    datPm <- log2(pm(dat))
		    M <- datPm[,,"channel1"] - datPm[,,"channel2"]
		    ctrl <- M[getControlIndex(dat), ]
		    if (onlyGood) {
		        good <- getControlIndex(dat, controlProbes="CONTROL_REGIONS", onlyGood=TRUE, matrix=TRUE) 
		        ctrl[!good] <- NA    
		    }
		    med <- apply(ctrl, 1, median, na.rm=TRUE)
		    for (i in 1:ncol(M)) {
		        map <- getMap(ctrl[,i], med)
		        M[,i] <- map(M[,i])   
		    }
		    datPm[,,"channel2"] <- datPm[,,"channel1"] - M
		    pm(dat) <- 2^datPm
		    return(dat)
		}

		normalizeBetweenSamples <- function(dat, untreated="allQuantiles", enriched="control") {    
		    datPm <- log2(pm(dat)) 
		    M <- datPm[,,"channel1"] - datPm[,,"channel2"]
		    ## Normalize untreated (leaving M unchanged)
		    if (untreated=="control") {
		        controlIndex <- getControlIndex(dat)
		        datBg <- log2(bg(dat)) 
		        datPm[,,"channel1"] <- normControlBg(pms=datPm[,,"channel1"], bgs=datBg[,,"channel1"], controlIndex=controlIndex)
		    } else if (untreated=="complete") {
		        m <- rowMedians(datPm[,,"channel1"])
		        datPm[,,"channel1"] <- matrix(rep(m,dim(dat)["Samples"]), nrow=length(m))
		    } else {
		        datPm[,,"channel1"] <- normQuantile(datPm[,,"channel1"], method=untreated)
		    }
		    datPm[,,"channel2"] <- datPm[,,"channel1"] - M
		    ## Normalize enriched 
		    #ngc <- countGC(dat, type="pm")
		    if (enriched=="control") {
		        controlIndex <- getControlIndex(dat)
		        datBg <- log2(bg(dat)) 
		        #bgNgc <- countGC(dat, type="bg")
		        #ctrlNgc <- countGC(dat, type="pm")[controlIndex] 
		        datPm[,,"channel2"] <- normControlBg(pms=datPm[,,"channel2"], bgs=datBg[,,"channel2"], controlIndex=controlIndex)
		        #M <- datPm[,,"channel1"] - datPm[,,"channel2"]
		    } else if (enriched=="ctrl") {
		        controlIndex <- getControlIndex(dat) 
		        datPm[,,"channel2"] <- normControlBg(pms=datPm[,,"channel2"], controlIndex=controlIndex)
		    } else if (enriched=="bg") {
		        datBg <- log2(bg(dat)) 
		        #bgNgc <- countGC(dat, type="bg")
		        datPm[,,"channel2"] <- normControlBg(pms=datPm[,,"channel2"], bgs=datBg[,,"channel2"])
		    } else {
		        datPm[,,"channel2"] <- normQuantile(datPm[,,"channel2"], method=enriched)
		    }
		    pm(dat) <- 2^datPm
		    return(dat)
		} 

		smooth <- function(dat) {
		    # Code taken from ~pmurakam/feinberg/CharmFiles/preprocessing.R
		    tmpM <- log(dat$p) - log(1-dat$p) # What is tmpM vs M?
		    INDEX <- dat$INDEX
		    pns <- dat$pns
    
		    lens <- sapply(INDEX,length)==1
		    if(any(lens)) message("AT LEAST 1 GROUP HAS ONLY 1 ARRAY!")
		    MM <- matrix(NA,nrow(tmpM),ncol=length(INDEX))
		    sigmas <- matrix(NA,nrow(tmpM),ncol=length(INDEX))
		    for(i in seq(along=INDEX)){
		      if(!lens[i]){
		          MM[,i] <- rowMeans(tmpM[,INDEX[[i]]])
		          sigmas[,i]=rowSds(tmpM[,INDEX[[i]]])^2
		      } else{
		          MM[,i] <- tmpM[,INDEX[[i]]]
		      }
		    }
       
		    NT=length(INDEX)
		    sMM <- matrix(NA,nrow(MM),NT)
		    colnames(sMM) <- names(INDEX) #i added this.
		    ssigmas <- matrix(NA,nrow(MM),NT)
		    FN <- 3
		    Tukey = function(x) pmax(1 - x^2,0)^2
		    fs= Tukey(seq(-FN,FN)/(FN+1));fs=fs/sum(fs)
		    Indexes=split(seq(along=pns),pns)
		        #seq(along=pns) indexes the rows of M, tmpM, MM, and sMM since
		        #pns[i] is the region that row i is a probe for.
		    for(i in seq(along=Indexes)){
		      #if(i%%1000==0) cat(paste(i,", ",sep=""))
		      Index=Indexes[[i]]
		      #x=pos[Index]  #don't know why this is here.
		      for(j in 1:NT){
		        sMM[Index,j] <- myfilter(MM[Index,j],fs)
		        if(!lens[j]) ssigmas[Index,j] <- myfilterse(sigmas[Index,j],fs,length(INDEX[[j]]))
		      }
		    }
		    dat$sMM <- sMM
		    dat$p.smooth.med <- exp(sMM) / (1+exp(sMM))
		    dat$ssigmas <- ssigmas
		    return(dat)
		}

		myfilter <- function(x,filter,...){
		  res=filter(x,filter,...)
		  ##now fill out the NAs
		  M=length(filter)
		  N=(M- 1)/2
		  L=length(x)
		  for(i in 1:N){
		    w=filter[(N-i+2):M]
		    y=x[1:(M-N+i-1)]
		    res[i] = sum(w*y)/sum(w)

		    w=rev(w)
		    ii=(L-(i-1))
		    y=x[(ii-N):L]
		    res[ii]<-sum(w*y)/sum(w)
		  }
		  return(res)
		}

		myfilterse <- function(x,filter,n,...){
		  ##filter must add up to 1!!!
		  res=filter(x,filter^2,...)
		  ##now fill out the NAs
		  M=length(filter)
		  N=(M- 1)/2
		  L=length(x)
		  for(i in 1:N){
		    w=filter[(N-i+2):M]
		    y=x[1:(M-N+i-1)]
		    res[i] = sum(w^2*y)/sum(w)^2

		    w=rev(w)
		    ii=(L-(i-1))
		    y=x[(ii-N):L]
		    res[ii]<-sum(w^2*y)/sum(w)^2
		  }
		  return(res/n)
		}
    
}