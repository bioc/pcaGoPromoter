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
            inputType == "hugene10st" | inputType == "hugene11st" |
            inputType == "mogene10st" | inputType == "mogene20st") {

        # Set org
        if( inputType == "hgu133plus2" | inputType == "hugene10st" |
            inputType == "hugene11st" ) {
            org <- "Hs"
        }        
        if( inputType == "mouse4302" | inputType == "mogene10st" |
	    inputType == "mogene20st") {
            org <- "Mm"
        }        
        if( inputType == "rat2302") {
            org <- "Rn"
        }        
        
        if( inputType == "hugene10st" | inputType == "hugene11st" | 
	    inputType == "mogene10st" | inputType == "mogene20st") {
            inputType <- paste(sep="", inputType, "transcriptcluster")
        }		
        packageName <- paste(sep="", inputType, ".db")
        require(packageName, character.only=TRUE, quietly=TRUE)

        if(outputType == "refseq") {
            REFSEQ <- get(paste(sep="", inputType, "REFSEQ"))	
#            REFSEQ <- mappedkeys(REFSEQ)
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


    } else if( inputType == "geneSymbol" ) {
        packageName <- paste(sep="", "org.", org, ".eg.db")
        require(packageName,character.only=TRUE,quietly=TRUE)
        ALIAS2EG <- get( paste( sep="" , "org.", org , ".egALIAS2EG" ) )  
        indexExists <- input %in% AnnotationDbi::keys(ALIAS2EG)
        if( !all(indexExists) ) {
            print("The following gene symbol cannot be mapped and are excluded:")
            print(input[!indexExists])
            input <- input[indexExists]
        }
        if( outputType == "refseq" ) {
            # TODO SYMBOL2EG eller ALIAS2EG ???
            SYMBOL2EG <- get( paste( sep="" , "org.", org , ".egALIAS2EG" ) ) 
            x <- unlist( AnnotationDbi::mget( input , SYMBOL2EG , ifnotfound=NA ) , use.names = FALSE)
            x <- x[ !is.na(x) ]
            REFSEQ <- get( paste( sep="" , "org.", org , ".egREFSEQ" ) )  
            x <- unlist( AnnotationDbi::mget( x , REFSEQ , ifnotfound=NA ) , use.names = FALSE)
            x <- unique(x[ !is.na(x) ])
            x <- x[grep("NM_",x)]
        } else {
            ALIAS2EG <- get( paste( sep="" , "org.", org , ".egALIAS2EG" ) )  
            x <- unlist( AnnotationDbi::mget( input , ALIAS2EG , ifnotfound=NA ) , use.names = FALSE)
            x <- unique(x[ !is.na(x) ])
        }
    } else if( inputType == "entrezID" ) {
        packageName <- paste(sep="", "org.", org, ".eg.db")
        require(packageName,character.only=TRUE,quietly=TRUE)
        if( outputType == "refseq" ) {
            REFSEQ <- get( paste( sep="" , "org.", org , ".egREFSEQ" ) )  
            x <- unlist( AnnotationDbi::mget( input , REFSEQ , ifnotfound=NA ) , use.names = FALSE)
            x <- unique(x[ !is.na(x) ])
            x <- x[grep("NM_",x)]
        } else {
            x <- input
        }
        
    
    } else {
        stop(paste(sep="", "Unknown inputType '", inputType, "'"))
    }
    return(list(x=x, org=org))
}


pcaInfoPlot <- function(eData, inputType="hgu133plus2", org="Hs", groups, 
        printNames=TRUE, plotCI=TRUE, noProbes=1365, GOtermsAnnotation=TRUE,
        primoAnnotation=TRUE ) {

 #   require(affy)
    
    pcaObj <- pca(eData)

    if( !is.na(pcaObj$expressionData[1]) ) {
        chipType <- pcaObj$expressionData[[1]]
 
        # Set inputType and org
        if( chipType == "hgu133plus2" ) {
            inputType <- chipType
            org <- "Hs"
        }
        if( chipType == "hugene10st" | chipType == "hugene10stv1" ) {
            inputType <- "hugene10st"
            org <- "Hs"
        }        
        if( chipType == "mouse4302" ) {
            inputType <- chipType
            org <- "Mm"
        }
        if( chipType == "mogene10st" | chipType == "mogene10stv1" ) {
            inputType <- "mogene10st"
            org <- "Mm"
        }        
        if( chipType == "mogene20st" ) {
            inputType <- "mogene20st"
            org <- "Mm"
        }        
        if( chipType == "rat2302") {
            inputType <- chipType
            org <- "Rn"
        }
    }
    
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
    plot(pcaObj, groups=groups, printNames=printNames, plotCI=plotCI,
         GOtreeObjs=GOtreeObjs, primoObjs=primoObjs)
}
