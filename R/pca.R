dropVarNull <- function(x) {
    for(i in 1:length(x[, 1])) {
        if(var(x[i, ]) == 0) {
            x[i, ] <- NA
            print(paste("Dropped variable", rownames(x)[i],
                "because scaling is not possible (all values are equal)"))
        }
    }
    x <- na.omit(x)
}

pca <- function(eData, printDropped=TRUE, scale=TRUE, center=TRUE) {
    
    if(class(eData) == "ExpressionSet") {
        annotationData <- annotation(eData)
        phenoData <- phenoData(eData)
        exprsData <- exprs(eData)
        expressionData <- list(annotationData,phenoData)
    } else { 
        exprsData <- as.matrix(eData)
        expressionData <- NA
    }
    
    for(i in 1:length(exprsData[, 1])) {
        if(var(exprsData[i, ]) == 0) {
            exprsData[i, ] <- NA
            if(printDropped) {
                print(paste("Dropped variable", rownames(exprsData)[i],
                    "because scaling is not possible (all values are equal)"))
            }
        }
    }
    exprsData <- na.omit(exprsData)
    
    p <- prcomp(t(exprsData), scale=scale, center=center)
    
    res <- list(scores=p$x, loadings=p$rotation, pov=p$sdev^2/sum(p$sdev^2),
            expressionData=expressionData)
    class(res) <- "pca"
    res
}

print.pca <- function( x, ... ) {
	print(x$scores)	
}

getRankedProbeIds <- function( x , pc=1 , decreasing=TRUE ) {
	
	UseMethod("getRankedProbeIds")

}

getRankedProbeIds.pca <- function( x , pc=1 , decreasing=TRUE ) {
	return(rownames(x$loadings)[order(x$loadings[, pc], decreasing=decreasing)])
}

plot.pca <- function(x, groups, PCs=c(1, 2), printNames=TRUE, symbolColors=TRUE,
        plotCI=TRUE, GOtreeObjs=NA, primoObjs=NA, main, ... ) {
    
    if(is.list(x$expressionData)) {
        pData <- pData(x$expressionData[[2]])
        if(any(colnames(pData) == "groups")) {
            groups <- as.factor(pData[,"groups"])
        }
        if(any(colnames(pData) == "group")) {
            groups <- as.factor(pData[,"group"])
        }
        if(any(colnames(pData) == "class")) {
            groups <- as.factor(pData[,"class"])
        }
        p <- x
    } else {
        p <- x
    }
    
    if(missing(groups)) {
        groups <- as.factor(rep("All", NROW(p$scores)));
    } else {
        if(!is.factor(groups)) {
            stop("Variable 'groups' must be a factor.")
        }		
    } 
    
    if(NROW(p$scores) != length(groups)) {
        warning(
    "'groups' should have same length as the number of observation in 'x'")
    }
    
    plotColors <- rep(c(2, 3, 4, 5, 6), 10)
    plotColorsEllipse <- rep(c("#FFAAAA80","#AAFFAA80","#AAAAFF80","#AAFFFF80",
            "#FFAAFF80"),10)
    savePar <- par(no.readonly=TRUE)
    classColors <- plotColors[as.numeric(groups)]
    layout(c(1, 2), heights=c(8, 1))
    par(srt=0, xpd=NA, oma=c(0, 0, 0, 0), mar=c(2, 3, 4, 4),
            mgp=c(1.75, 0.5, 0))

    if(missing(main)) {
        main <- paste(sep="", "PCA plot of ", PCs[1], 
                ". and ", PCs[2], ". principal component")
    }
    if(is.list(GOtreeObjs) | is.list(primoObjs)) {
        main <- NULL		
    }
    
    ## require("ellipse")
    
    ## Find plot area 
    ellipseData <- ellipse(var(p$scores[, PCs]),
            centre=colMeans(p$scores[, PCs]))
    obsData <- p$scores[, PCs]
    xValues <- c(max(abs(c(ellipseData[, 1],obsData[, 1])))*-1,
            max(abs(c(ellipseData[, 1],obsData[, 1]))))
    yValues <- c(max(abs(c(ellipseData[, 2],obsData[, 2])))*-1,
            max(abs(c(ellipseData[, 2],obsData[, 2]))))
    
    plot(xValues, yValues, main=main,
            xlab=paste(sep="", "PC", PCs[1], "  -  Proportion of variance: ",
                    round(100*p$pov[PCs[1]])/100),
            ylab=paste(sep="", "PC", PCs[2], "  -  Proportion of variance: ",
                    round(100*p$pov[PCs[2]])/100), type="n")	

    par(xpd=FALSE)

    if(nlevels(groups) > 1 & plotCI) {
        for(i in 1:nlevels(groups)) {
            if(summary(groups)[i] > 2) {
                polygon(ellipse(var(p$scores[groups==levels(groups)[i], PCs]),
                    centre=colMeans(p$scores[groups==levels(groups)[i], PCs])),
                    col=plotColorsEllipse[i], border=FALSE)
            }
        }
    }
    
    lines(ellipse(var(p$scores[, PCs]), centre=colMeans(p$scores[, PCs])),
            col="#AAAAAA")
    
    if(symbolColors) {
        symbolColors <- classColors		
    } else {
        symbolColors <- c(1)
    }
    points(p$scores[, PCs], pch=as.numeric(groups), col=symbolColors, cex=0.75,
            xlab=paste(sep="", "PC", PCs[1], "  -  Proportion of variance: ",
                    round(100*p$pov[PCs[1]])/100),
            ylab=paste(sep="", "PC", PCs[2], "  -  Proportion of variance: ",
                    round(100*p$pov[PCs[2]])/100))

    abline(h=0, v=0, lty=3)
    box()
    
    ## Plot the observation names
    if( printNames ) {
        text(p$scores[, PCs], rownames(p$scores), col=classColors, cex=1, pos=1,
                offset=0.35)
    }	
    
    ## Print annotation on axis?
    if(is.list(GOtreeObjs) | is.list(primoObjs)) {
        printAnnotation(GOtreeObjs, primoObjs)
    }
    
    if(!is.null(levels(groups))) {
        # class data contains no levels, so no legend
        par(xpd=TRUE, mar=c(0, 0, 1.5, 0))
        plot.new(); 
        legend("top", legend=levels(groups), horiz=TRUE, x.intersp=0.5,
                pch=seq(1, nlevels(groups)), text.col=plotColors)
    }
    
    par(savePar)

}

printAnnotation <- function(GOtreeObjs, primoObjs) {
    
    if( is.list(GOtreeObjs) & is.list(primoObjs) ) {		
        text <- paste(sep="",
                format("Gene Ontology BP terms",width=33), # was 35
                format("Poss. TFBS",width=12), # was 15
                format("",width=5),
                format("Gene Ontology BP terms",width=33),
                format("Poss. TFBS",width=12)
        )
    } else if( is.list(GOtreeObjs) ) {		
        text <- paste(sep="",
                format("Gene Ontology BP terms",width=40),
                format("",width=25),
                format("Gene Ontology BP terms",width=40)
        )
    } else if( is.list(primoObjs) ) {		
        text <- paste(sep="",
                format("Poss. TFBS",width=15),
                format("",width=65),
                format("Poss. TFBS",width=15)
        )
    }
    
    ## Top
    mtext(text, outer=FALSE, side=3, line=2.5, adj=0, cex=0.5, family="mono",
            font=4)
    printAnnText(GOtreeObjs=GOtreeObjs, primoObjs=primoObjs, lstart=2,
            ladd=-0.5, side=3, negIndex=2, posIndex=1)
    
    ## Right
    mtext(text, outer=FALSE, side=4, line=-0.5, adj=0, cex=0.5, family="mono",
            font=4)
    printAnnText(GOtreeObjs=GOtreeObjs, primoObjs=primoObjs, lstart=0,
            ladd=0.5, side=4, negIndex=4, posIndex=3)

}

printAnnText <- function(GOtreeObjs, primoObjs, lstart=2, ladd=-0.5, side=3,
        negIndex=2, posIndex=1) {
    
    l=lstart
    for(i in 1:5) {
        GOtext1 <- character()
        GOtext2 <- character()
        primoText1 <- character()
        primoText2 <- character()
        
        if(is.list(GOtreeObjs) & is.list(primoObjs)) {		
            GOtreeWidth <- 33 # was 35
            primoWidth <- 12 # was 15
            spacer <- format("", width=5)
        } else if(is.list(GOtreeObjs)) {		
            GOtreeWidth <- 40
            spacer <- format("", width=25)
        } else if(is.list(primoObjs)) {		
            primoWidth <- 15
            spacer <- format("", width=65)
        }
        
        if(is.list(GOtreeObjs)) {
            GOtext1 <- GOtreeObjs[[negIndex]]$sigGOs$GOterm[i]
            if(nchar(as.character(GOtext1)) > GOtreeWidth-5) {
                GOtext1 <- paste(sep="", strtrim(GOtext1, GOtreeWidth-8), "...")				
            }
            GOtext1 <- format(GOtext1, width=GOtreeWidth)
            
            GOtext2 <- GOtreeObjs[[posIndex]]$sigGOs$GOterm[i]
            if( nchar(as.character(GOtext2)) > GOtreeWidth-5 ) {
                GOtext2 <- paste(sep="", strtrim(GOtext2, GOtreeWidth-8), "...")				
            }
            GOtext2 <- format(GOtext2, width=GOtreeWidth)
        }
        
        if(is.list(primoObjs)) {
            primoText1 <- primoObjs[[negIndex]]$overRepresented$gene[i]
            if(nchar(as.character(primoText1)) > primoWidth) {
                primoText1 <- paste(sep="", strtrim(primoText1, primoWidth-3),
                        "...")				
            }
            primoText1 <- format(primoText1, width=primoWidth)
            
            primoText2 <- primoObjs[[posIndex]]$overRepresented$gene[i]
            if(nchar(as.character(primoText2)) > primoWidth) {
                primoText2 <- paste(sep="", strtrim(primoText2, primoWidth-3),
                        "...")				
            }
            primoText2 <- format(primoText2, width=primoWidth)
            
        }	 			
        text <- paste(sep="", GOtext1, primoText1, spacer, GOtext2, primoText2)
        
        mtext(text, outer=FALSE, side=side, line=l, adj=0, cex=0.5,
                family="mono", font=2)
        l=l+ladd
    }	
}

pcaLeaveOut <- function( exprsData , leaveOut=0.1 , runs = NCOL(exprsData) ) {
    
    if(leaveOut == 1) {
        runs <- NCOL(exprsData)
        index <- lapply(1:runs, function(i) { c(1:NCOL(exprsData))[-i] })
    } else {	
        toTakeIn <- round(NCOL(exprsData)*(1-leaveOut))
        toTakeInIndex <- c(rep(TRUE, toTakeIn), rep(FALSE,
                        NCOL(exprsData)-toTakeIn))
        index <- lapply(1:runs, function(i) { sample(toTakeInIndex) })
    }
    
    func <- function(i) { pca(exprsData[, index[[i]]], printDropped=FALSE) }	
    if(require("parallel")) {
        res <- mclapply(1:runs, func, mc.set.seed=TRUE)
    } else {
        res <- lapply(1:runs, func)
    }
    return(res)	
}

multiPca <- function(exprsData, pc, noProbes, decreasing, leaveOut, runs,
        testFunc, inputType, org ) {
    
    pcaRes <- pcaLeaveOut(exprsData, leaveOut=leaveOut, runs=runs)
    
    ## Do we have a sign switching PCA result?
    pcaAll <- pca(exprsData)
    origSignScores <- pcaAll$scores[, pc] >= 0
    decreasingIndex <- logical(runs)
    for(i in 1:runs) {
        ## Get the sign vector with the observations in play
        actualSignScores <- origSignScores[rownames(pcaAll$scores) %in%
                        rownames(pcaRes[[i]]$scores) ]
        signScores <- pcaRes[[i]]$scores[, pc] >= 0
        
        noSignChange <- sum(!(xor(actualSignScores, signScores)))
        signChange <- sum(!(xor(actualSignScores, !(signScores))))
        
        if(noSignChange >= signChange) {
            decreasingIndex[i] <- TRUE
        }	else {
            decreasingIndex[i] <- FALSE
        }
    }
    
    if(decreasing == FALSE) {
        decreasingIndex <- !decreasingIndex
    }
    decreasingIndex <- as.character(factor(decreasingIndex,
                    labels=c("FALSE","TRUE")))
    
#    func <- function(i) {
#        testFunc(getRankedProbeIds(pcaRes[[i]], pc=pc,
#                 decreasing=decreasingIndex[i])[1:noProbes],
#                 inputType=inputType, org=org ) 
#    }
#    
#    ## multicore giver problemer med sqllite databaserne...	
#    if(require("parallel")) {
#        res <- mclapply( 1:runs, func , mc.set.seed = TRUE )
#    } else {
#        res <- lapply( 1:runs, func )
#    }

    res <- lapply( 1:runs, function(i) { 
                testFunc(getRankedProbeIds(pcaRes[[i]], pc=pc,
                decreasing=decreasingIndex[i])[1:noProbes],
                inputType=inputType, org=org ) } )

    return(res)
}
