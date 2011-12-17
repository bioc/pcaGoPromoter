convertFromProbeIdsToRefseqs <- function(x, chipType) {
    
    ## Dynamically create envir name. e.g. hgu133plus2SYMBOL
    REFSEQ <- get( paste( sep="" , chipType , "REFSEQ" ) )
    
    a <- unique(unlist(AnnotationDbi::mget(x, REFSEQ, ifnotfound=NA),
                    use.names=FALSE))
    a <- a[!is.na(a)]
    refseqType <- sapply(a, substr, 0, 2)
    nmIndex <- sapply(refseqType, '==', "NM")
    res <- a[nmIndex]
    
    return(res)
}

cleanUpTargetSet <- function(x, printIgnored=TRUE, primoData) {
    
    promotersRefseqs <- primoData[["promotersRefseqs"]]
    
    target <- x[ x %in% unlist(promotersRefseqs) ]
    
    diffTarget <- x[!(x %in% target)]
    if(length(diffTarget) > 0 & printIgnored) {
        cat(
"These refseqs could not be found in the background set and will be ignored:\n")
        print(diffTarget)
    }	
    
    return(target)
}

getMatricesWithHits <- function(target, pvalueCutOff=0.05, cutOff=0.9,
        p.adjust.method, primoData) {

    promotersRefseqs <- primoData[["promotersRefseqs"]]
    thresholds <- primoData[["thresholds"]]
    
    targetBackgroundIndex <- sapply(sapply(promotersRefseqs, "%in%", target),
            any)
    
    noTargets <- sum(targetBackgroundIndex)
    noPromoters <- length(promotersRefseqs)	
    
    resId <- character()
    resBaseId <- character()
    resPwmLength <- numeric()
    resGene <- character()
    resPvalue <- numeric()
    resUpDown <- logical()
    counter <- 0
    
    noPromotersCutOff <- round(noPromoters*cutOff)
    
    for(i in 1:length(thresholds)) {
        
        resAtCutOff <- sort( thresholds[[i]]$maxScores )[ noPromotersCutOff ]
        ## We got a round error somewhere,
        ## so we need to subtract some value form threshold.
        hitsIndex <- thresholds[[i]]$maxScores >= resAtCutOff - 0.000001
        
        targetHitsIndex <- hitsIndex & targetBackgroundIndex
        
        noHitsBackground <- sum(hitsIndex)
        noHitsTarget <- sum(targetHitsIndex)
        
        a <- noHitsBackground-noHitsTarget
        b <- noPromoters-noHitsBackground-(noTargets-noHitsTarget)
        c <- noHitsTarget
        d <- noTargets-c
        
        z <- matrix(c(a, b, c, d), nrow=2, ncol=2)
        pvalue <- fisher.test(z)$p.value
        
        distBackground <- a/noPromoters
        distTarget <- c/noTargets
        
        if( pvalue <= pvalueCutOff ) {
            counter <- counter+1
            resId[counter] <- names(thresholds[i])
            resBaseId[counter] <- thresholds[[i]][["baseId"]]
            resPwmLength[counter] <- thresholds[[i]][["pwmLength"]]
            resGene[counter] <- thresholds[[i]][["name"]]
            resPvalue[counter] <- pvalue
            
            if(distTarget > distBackground) {
                resUpDown[counter] <- TRUE
            } else {
                resUpDown[counter] <- FALSE
            }  
        }
    }
    
    resPvalueAdjusted <- p.adjust(resPvalue, method=p.adjust.method)
    
    indexUp <- order(resPvalueAdjusted[resUpDown])
    indexDown <- order(resPvalueAdjusted[!resUpDown])
    
    res <- list(overRepresented=data.frame(id=resId[resUpDown][indexUp],
                    baseId=resBaseId[resUpDown][indexUp],
                    pwmLength=resPwmLength[resUpDown][indexUp],
                    gene=resGene[resUpDown][indexUp],
                    pValue=resPvalueAdjusted[resUpDown][indexUp]),
                underRepresented=data.frame(id=resId[!resUpDown][indexDown],
                    baseId=resBaseId[!resUpDown][indexDown],
                    pwmLength=resPwmLength[!resUpDown][indexDown],
                    gene=resGene[!resUpDown][indexDown],
                    pValue=resPvalueAdjusted[!resUpDown][indexDown]))
    
    return(res)
}

primo <- function(input, inputType="hgu133plus2", org="Hs", PvalueCutOff=0.05,
        cutOff=0.9, p.adjust.method="fdr", printIgnored=FALSE,
        primoData=NULL) {
	
    convInput <- convertInput(input=input, inputType=inputType, org=org,
            outputType="refseq")
    
    refseqs <- convInput[["x"]]
    org <- convInput[["org"]]
    
    if(is.null(primoData)) {
        primoData <- new.env()	
        if(org == "Hs") {
            require("pcaGoPromoter.Hs.hg19")
            data("pcaGoPromoter.Hs.hg19",envir=primoData)
        } else if(org == "Mm") {
            require("pcaGoPromoter.Mm.mm9")
            data("pcaGoPromoter.Mm.mm9",envir=primoData)
        } else if(org == "Rn") {
            require("pcaGoPromoter.Rn.rn4")
            data("pcaGoPromoter.Rn.rn4",envir=primoData)
        } else {
            stop(paste("Unknown org",org))
        }		
        promotersRefseqs <- get("promotersRefseqs", pos=primoData)
        thresholds <- get("thresholds", pos=primoData)
        primoData <- list(thresholds=thresholds,
                promotersRefseqs=promotersRefseqs)
        class(primoData) <- "primoData"
    } else {
        if(class(primoData) != "primoData") {
            stop("primoData is not of class primoData")
        }
    }
    target <- cleanUpTargetSet(refseqs, printIgnored, primoData=primoData)
    
    res <- getMatricesWithHits(target, PvalueCutOff, cutOff=cutOff,
            p.adjust.method=p.adjust.method, primoData=primoData)
    return(res)	
}

primoHits <- function(input, inputType="hgu133plus2", org="Hs", id, cutOff=0.9,
        printIgnored=FALSE, primoData=NULL) {

    if (missing(id)) {
        stop("Argument 'id' must be set.")
    }
    
    convInput <- convertInput(input=input, inputType=inputType, org=org,
            outputType="refseq")
    refseqs <- convInput[["x"]]
    org <- convInput[["org"]]
    
    if(is.null(primoData)) {
        primoData <- new.env()	
        if(org == "Hs") {
            require("pcaGoPromoter.Hs.hg19")
            data("pcaGoPromoter.Hs.hg19", envir=primoData)
        } else if(org == "Mm") {
            require("pcaGoPromoter.Mm.mm9")
            data("pcaGoPromoter.Mm.mm9", envir=primoData)
        } else if(org == "Rn") {
            require("pcaGoPromoter.Rn.rn4")
            data("pcaGoPromoter.Rn.rn4", envir=primoData)
        } else {
            stop(paste("Unknown org",org))
        }		
        promotersRefseqs <- get("promotersRefseqs", pos=primoData)
        thresholds <- get("thresholds", pos=primoData)
        primoData <- list(thresholds=thresholds,
                promotersRefseqs=promotersRefseqs)
        class(primoData) <- "primoData"		
    } else {
        if(class(primoData) != "primoData") {
            stop("primoData is not of class primoData")
        }
    }
    target <- cleanUpTargetSet(refseqs, printIgnored, primoData=primoData)
    
    promotersRefseqs <- primoData[["promotersRefseqs"]]
    thresholds <- primoData[["thresholds"]]

    noPromoters <- length(promotersRefseqs)	
    noPromotersCutOff <- round(noPromoters*cutOff)
    resAtCutOff <- 
            sort(thresholds[[as.character(id)]]$maxScores)[ noPromotersCutOff ]
    ## We got a round error somewhere,
    ## so we need to subtract some form threshold.
    hitsIndex <-
            thresholds[[as.character(id)]]$maxScores >= resAtCutOff-0.000001
    targetBackgroundIndex <- sapply(sapply(promotersRefseqs, "%in%", target),
            any)
    targetHitsIndex <- hitsIndex & targetBackgroundIndex
    
    allRefseqs <- unlist(promotersRefseqs[targetHitsIndex])

    res <- target[target %in% allRefseqs] 

    return(res)
}

primoWithLeaveOut <- function(exprsData, inputType="hgu133plus2", org ="Hs",
        pc=1, decreasing=TRUE, noProbes=1000, leaveOut=1,
        runs=NCOL(exprsData)) {
    
    multiRes <-	multiPca( exprsData=exprsData, pc=pc , noProbes=noProbes,
            decreasing=decreasing, leaveOut=leaveOut, runs=runs, testFunc=primo,
            inputType=inputType, org=org)
    
    listName <- c("overRepresented")
    leaveOutRes <- leaveOutFindPvalues(res=multiRes, listName=listName,
            uniqueId="id")
    overRepRes <- data.frame(
            id=multiRes[[1]][[listName]][leaveOutRes[["index"]], "id"],
            baseId=multiRes[[1]][[listName]][leaveOutRes[["index"]], "baseId"],
            pwmLength=multiRes[[1]][[listName]][leaveOutRes[["index"]],
                    "pwmLength"],
            gene=multiRes[[1]][[listName]][leaveOutRes[["index"]], "gene"],
            pValue=leaveOutRes[["medianPvalues"]])
    
    listName <- c("underRepresented")
    leaveOutRes <- leaveOutFindPvalues(res=multiRes, listName=listName,
            uniqueId="id")
    underRepRes <- data.frame(
            id=multiRes[[1]][[listName]][leaveOutRes[["index"]], "id"],
            baseId=multiRes[[1]][[listName]][leaveOutRes[["index"]], "baseId"],
            pwmLength=multiRes[[1]][[listName]][leaveOutRes[["index"]],
                    "pwmLength"],
            gene=multiRes[[1]][[listName]][leaveOutRes[["index"]], "gene"],
            pValue=leaveOutRes[["medianPvalues"]])
    
    finalRes <- list(overRepresented=overRepRes[order(overRepRes$pValue), ],
            underRepresented=underRepRes[order(underRepRes$pValue), ])
    
    return(finalRes)
}
