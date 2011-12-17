leaveOutFindPvalues <- function(res, listName, uniqueId) {
    
	  runs <- length(res)
    commonIds <- as.character(res[[1]][[listName]][[uniqueId]])
    for(i in 2:runs) {
        commonIds <- commonIds[commonIds %in% res[[i]][[listName]][[uniqueId]] ]
    }
    commonIds <- as.character(commonIds)
    if(length(commonIds) == 0) {
        return(list(medianPvalues=numeric(),index=rep(FALSE,
                length(res[[1]][[listName]][[uniqueId]]))))
    }
    
    medianPvalues <- sapply(1:length(commonIds), function(i) { 
        median( sort( sapply(1:runs, function(j) {
            res[[j]][[listName]][["pValue"]][ as.character(res[[j]]
                    [[listName]][[uniqueId]]) == commonIds[i]]
            } ) ) )
        }
    )			
    index <- sapply(1:NROW(res[[1]][[listName]]), function(i) {
        any(commonIds %in% as.character(res[[1]][[listName]][[uniqueId]][i]))
    } )
    return(list(medianPvalues=medianPvalues, index=index))
}

convertInput <- function(input, inputType, org="Hs", outputType="entrezId") {
    
    if(inputType == "hgu133plus2" | inputType == "mouse4302" |
            inputType == "rat2302" |
            inputType == "hugene10st" | inputType == "mogene10st") {

        if( inputType == "hugene10st" | inputType == "mogene10st") {
            inputType <- paste(sep="", inputType, "transcriptcluster")
        }		
        packageName <- paste(sep="", inputType, ".db")
        require(packageName, character.only=TRUE, quietly=TRUE)

        if(outputType == "refseq") {
            REFSEQ <- get(paste(sep="", inputType, "REFSEQ"))	
            x <- unlist(AnnotationDbi::mget(input, REFSEQ, ifnotfound=NA),
                    use.names=FALSE)
            x <- unique(x[!is.na(x)])
            x <- x[grep("NM_",x)]
        } else { 
            ENTREZID <- get(paste( sep="", inputType, "ENTREZID"))	
            x <- unlist(AnnotationDbi::mget(input, ENTREZID, ifnotfound=NA),
                    use.names=FALSE)
            x <- unique(x[!is.na(x)])
        }
    } else {
        stop(paste(sep="", "Unknown inputType '", inputType, "'"))
    }
    return(list(x=x, org=org))
}


pcaInfoPlot <- function(exprsData, inputType="hgu133plus2", org="Hs", groups, 
        noProbes=1365, GOtermsAnnotation=TRUE, primoAnnotation=TRUE) {

    pcaObj <- pca(exprsData)

    probesPC1pos <- getRankedProbeIds(pcaObj, pc=1,
            decreasing=TRUE)[1:noProbes]
    probesPC1neg <- getRankedProbeIds(pcaObj, pc=1,
            decreasing=FALSE)[1:noProbes]
    probesPC2pos <- getRankedProbeIds(pcaObj, pc=2,
            decreasing=TRUE)[1:noProbes]
    probesPC2neg <- getRankedProbeIds(pcaObj, pc=2,
            decreasing=FALSE)[1:noProbes]

    if(GOtermsAnnotation) {
        GOtreePC1pos <- GOtree(probesPC1pos, inputType=inputType)
        GOtreePC1neg <- GOtree(probesPC1neg, inputType=inputType)
        GOtreePC2pos <- GOtree(probesPC2pos, inputType=inputType)
        GOtreePC2neg <- GOtree(probesPC2neg, inputType=inputType)
        GOtreeObjs <- list(GOtreePC1pos, GOtreePC1neg, GOtreePC2pos, 
                GOtreePC2neg)
    } else {
        GOtreeObjs <- NA
    }
    if(primoAnnotation) {
        primoPC1pos <- primo(probesPC1pos, inputType=inputType, org=org)
        primoPC1neg <- primo(probesPC1neg, inputType=inputType, org=org)
        primoPC2pos <- primo(probesPC2pos, inputType=inputType, org=org)
        primoPC2neg <- primo(probesPC2neg, inputType=inputType, org=org)
        primoObjs <- list(primoPC1pos, primoPC1neg, primoPC2pos, primoPC2neg)
    } else {
        primoObjs <- NA
    }
	plot(pcaObj, groups=groups, GOtreeObjs=GOtreeObjs, primoObjs=primoObjs)
}
