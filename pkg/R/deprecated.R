#############################################
### Old bits of code not currently in use ###
#############################################
if (FALSE) { 


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
		return(ynorm)
	}


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