
getGOidsInOntology <- function(ontology="BP", org="Hs", noIEA=TRUE) {
    
    GO <- get(paste(sep="", "org.", org, ".egGO"))
    
    ## Get the entrez gene identifiers that are mapped to a GO ID
    mapped_genes <- mappedkeys(GO)
    # Convert to a list
    x <- as.list(GO[mapped_genes])
    
    # Get all GO terms in one list
    x <- unlist(x, recursive=FALSE, use.names=FALSE)
    
    OntologyIndex <- sapply(x, "[", i="Ontology") == ontology
    if(noIEA) {
        noIEAindex <- sapply(x, "[", i="Evidence") != "IEA"
    } else {
        noIEAindex <- TRUE
    }
    
    GOids <- unique(
            unlist(sapply(x, "[", i="GOID")[OntologyIndex & noIEAindex]))
    return(GOids)
}

getGOidsWithChildren <- function(specificGOids, org="Hs", noIEA=TRUE) {
    
    GO2ALLEGS <- get(paste(sep="", "org.", org, ".egGO2ALLEGS"))
    
    x <- as.list(GO2ALLEGS)
    
    GOids <- x[names(x) %in% specificGOids]
    
    GOidsNames <- names(GOids)
    
    if(noIEA) {
        noIEAindex <- sapply(GOids, function(x) { names(x) != "IEA" } )
        GOids <- lapply(1:length(GOids), function(x) {
                    unique(GOids[[x]][noIEAindex[[x]]]) } )
    }
    
    noDuplets <- lapply(1:length(GOids), function(x) {
                GOids[[x]][!duplicated(GOids[[x]])] } )
    
    names(noDuplets) <- GOidsNames
    
    return(noDuplets)
}

getGOsDownTree <- function(x, org="Hs", terms=TRUE) {

    GOTermsOntologyBP <- getGOidsWithChildren(getGOidsInOntology(org=org),
            org=org)
    
    noGenesInBP <- sapply(sapply(GOTermsOntologyBP, "%in%", x), sum)
    indexGOs <- noGenesInBP > 0
    
    GOid <- names(GOTermsOntologyBP)[indexGOs]
    genesInTerm <- noGenesInBP[indexGOs]
    names(genesInTerm) <- NULL	
    totalGenesInTerm <- sapply(GOTermsOntologyBP, length)[indexGOs]
    
    if (terms) {
        require("GO.db")
        GOterms <- as.list(GOTERM)
        GOIDtoGOTERM <- sapply(1:length(GOterms), function(i) {
                    Term(GOterms[[i]]) } )
        names(GOIDtoGOTERM) <- sapply(1:length(GOterms), function(i) {
                    GOID(GOterms[[i]]) } )		
        GOterm <- sapply(1:length(GOid), function(i) {
                    GOIDtoGOTERM[GOid[i]] } )
        names(GOterm) <- NULL
        return(data.frame(GOid, genesInTerm, totalGenesInTerm, GOterm,
                        row.names=NULL))
    } else {
        return(data.frame(GOid, genesInTerm, totalGenesInTerm, row.names=NULL))
    }
}

getSigGOTermsOnce <- function(a, noGenesInCut, org="Hs",
        statisticalTest="binom", binomAlpha=NA, p.adjust.method="fdr") {
    
    GO2ALLEGS <- get(paste(sep="", "org.", org, ".egGO2ALLEGS"))
    
    noAllGenes <- length(unique(unlist(as.list(GO2ALLEGS))))
    
    pValues <- numeric(NROW(a))
    
    if(statisticalTest == "fisher") {
        for(i in 1:NROW(a)) {
            genesInTerm <- a$genesInTerm[i]
            totalGenesInTerm <- a$totalGenesInTerm[i]
            genesNotInTerm <- totalGenesInTerm-genesInTerm
            z <- c(genesInTerm, noGenesInCut-genesInTerm, genesNotInTerm,
                    noAllGenes-genesNotInTerm-noGenesInCut)
            dim(z) <- c(2, 2)
            pValues[i] <- fisher.test(z)$p.value
        }
    } else if(statisticalTest == "binom") {
        ## Is a pValue given as argument binomAlpha?
        ## If not, set the pvalue to noGenesInCut/noAllGenes
        if(is.na(binomAlpha)) {
            binomAlpha <- noGenesInCut/noAllGenes
        }	
        pValues <- sapply(1:NROW(a), function(i) {
                    binom.test(a$genesInTerm[i], a$totalGenesInTerm[i],
                            p=binomAlpha, alternative="greater")$p.value } )
    } else {
        stop(paste("Unknown statistical test",statisticalTest))		
    }
    
    pValues <- p.adjust(pValues, p.adjust.method)	
    return(data.frame(GOid=a$GOid, genesInTerm=a$genesInTerm,
                    totalGenesInTerm=a$totalGenesInTerm, pValue=pValues,
                    GOterm=a$GOterm))
}

indexOfDupletsList <- function(a) {
    ## Get index of duplets and return a list duplets and position
    
    b <- list()
    n <- character()
    
    for (i in 1:length(a)) {
        c <- which(a == a[i])
        if (length(c) > 1) {
            b <- append(b, list(c))
            n <- append(n, a[i])
            a[c] <- NA
        }
    }
    
    names(b) <- n
    return(b)
}

newSearchGOtreeForTerms <- function(a, topNode="GO:0008150", limit=25) {
    ## Finds parent GO ids in GO tree. Only for visual purposes.
    
    b <- a[order(a$pValue)[1:limit], ]
    
    if(length(b$GOid) == 0 | limit == 0) {
        ## No GO terms...
        parent <- character()
        child <- character()
        genesInTerm <- numeric()
        totalGenesInTerm <- numeric()
        pValue <- numeric()
        GOterm <- character()
        return(NULL)
    }
    
    GOid <- as.character(b$GOid)
    parent <- character()
    child <- character()
    saved <- logical()
    replaced <- logical()
    
    saveEntry <- function(newParent, newChild) {
        assign("parent", append(parent, newParent), inherits=TRUE)
        assign("child", append(child, newChild), inherits=TRUE)
    }
    
    require("GO.db")
    GOBPAncestorList <- as.list(GOBPANCESTOR)
    
    ## Find parents
    for( i in 1:length(GOid) ) {
        possibleParents <- GOBPAncestorList[[GOid[i]]]		
        possibleParents <- possibleParents[possibleParents != "all"]		
        if(length(possibleParents) == 0) {
            next;
        }
        saved <- FALSE
        for(j in 1:length(GOid)) {
            if(GOid[j] %in% possibleParents) {
                ## Parent found ( GOid[j] )
                parent <- c(parent, GOid[j])
                child <- c(child, GOid[i])
                saved <- TRUE
            }
        }
        if(!saved) {
            parent <- c(parent, topNode)
            child <- c(child, GOid[i])
        }
    }
    
    x <- indexOfDupletsList(child)
    if (length(x)) {
        keep <- rep(TRUE, length(child))
        for(i in 1:length(x)) {
            for(j in 1:length(x[[i]])) {
                possibleParents <- GOBPAncestorList[parent[x[[i]][j]]][[1]]
                for(k in 1:length(x[[i]])) {
                    if(j == k | !x[[i]][k]) {
                        next
                    }
                    if(parent[x[[i]][k]] %in% possibleParents) {
                        keep[x[[i]][k]] <- FALSE
                    }
                }
            }
        }
        child <- child[keep]
        parent <- parent[keep]
    }
    
    noChildren <- length(child)
    
    genesInTerm <- numeric()
    totalGenesInTerm <- numeric()
    pValue <- numeric()
    GOterm <- character()
    
    for( i in 1:length(child)) {
        genesInTerm[i] <- a$genesInTerm[a$GOid == child[i]]
        totalGenesInTerm[i] <- a$totalGenesInTerm[a$GOid == child[i]]
        pValue[i] <- a$pValue[a$GOid == child[i]]
        GOterm[i] <- as.character(a$GOterm[a$GOid == child[i]])
    }
    
    return(data.frame(parent, child, genesInTerm, totalGenesInTerm, pValue,
                    GOterm))
}


makeGOtreeRgraphviz <- function (a,only1BoxPrTerm=TRUE,
        legendPosition="topright",...) {
    ## Makes Graphviz file and runs Graphviz using system() command
    
    if(length(a$child) == 0) {
        return()
    }
    
    if(!require("Rgraphviz")) {
        stop(
        "Package Rgraphiz is not installed, so cannot make graphical GO tree.")
    }

    nodeLabels <- character()
    nodeColors <- character()
    
    g <- new("graphNEL", edgemode="directed")
    
    noLegendValues <- 12
    legendValues <- numeric(noLegendValues)
    legendColors <- character(noLegendValues)
    color <- 192
    value <- 1
    for(n in 1:noLegendValues) {
        legendValues[n] <- value
        legendColors[n] <- rgb(192, color, color, maxColorValue=255)
        color <- color-16
        value <- value/100
    }
    
    g <- addNode("n0", g)  
    nodeLabels <- append(nodeLabels, paste(sep="", "GO:0008150", "\n",
                    "Biological Process"))
    nodeColors <- append(nodeColors, c("#ffffff"))
    
    if(only1BoxPrTerm) {
        b <- a[!duplicated(a$child), ]
    } else {
        b <- a
    }
    
    # Label all boxes
    for( i in 1:length(b$child)) {
        GOid <- b$child[i]
        genesInTerm <- b$genesInTerm[i]
        totalGenesInTerm <- b$totalGenesInTerm[i]
        GOterm <- b$GOterm[i]
        pValue <- sprintf("%.2e", b$pValue[i])
        
        pValueAndGenes <- paste(sep="", "\nP: ", pValue, " Genes: ",
                genesInTerm, "/", totalGenesInTerm)
        
        GOtermWrapped <- strwrap(GOterm, width=30)
        GOtermNice <- character()
        for(j in 1:min(length(GOtermWrapped), 2)) {
            GOtermNice <- paste(sep="", GOtermNice, "\n", GOtermWrapped[j])
        }
        if(length(GOtermWrapped) > 2) {
            GOtermNice <- paste(sep="", GOtermNice, "...")
        } 
        g <- addNode( paste(sep="", "n", i), g)
        nodeLabels <- append(nodeLabels, paste(sep="", GOid, GOtermNice,
                        pValueAndGenes))
        
        pValue <- as.numeric(pValue)
        for(n in 1:noLegendValues) {
            if(legendValues[n] < pValue) {
                break;
            }
        }
        if(n == 1) {
            ## p-value > 1
            nodeColors <- append(nodeColors, "#ffffff")
        } else {
            nodeColors <- append(nodeColors, legendColors[n-1])
        }
    }
    
    for(i in 1:length(a$parent)) {
        parent <- match(a$parent[i], b$child)
        
        if(only1BoxPrTerm) {
            child <- match(a$child[i], b$child)
        } else {
            child <- i
        }
        
        if(is.na(parent)) {
            ## Children with NA parents has top node as parent
            parent <- 0
        }
        g <- addEdge(paste(sep="", "n", parent), paste(sep="", "n", child), g)
    }
    
    names(nodeLabels) <- nodes(g)
    names(nodeColors) <- nodes(g)
    attrs <- list(graph=list(
                    rankdir = "TB",
                    nodesep = 0.2
            )
    )
    nodeRenderInfo(g) <- list(label=nodeLabels, iwidth=1.5, iheight=1,
            shape="box", fixedsize=FALSE, lwd=0.5, fill=nodeColors, fontsize=18)
    edgeRenderInfo(g) <- list(lwd=0.5)
    parRenderInfo(g) <- list(graph=list(...))
    
    savePar <- par(no.readonly=TRUE)
    
    par("mar"=c(3, 3, 3, 3), "cex"=0.8)  
    renderGraph(layoutGraph(g, attrs=attrs))
    
    if(!is.null(legendPosition)) {
        legend(legendPosition, legend=format(legendValues,scientific=TRUE),
                horiz=FALSE, fill=legendColors, title="p-value <")
    }
    par(savePar)
}

GOtree <- function(input, inputType="hgu133plus2", org="Hs",
        statisticalTest="binom", binomAlpha=NA, p.adjust.method="fdr") {
    
    a <- convertInput(input=input, inputType=inputType, org=org)
    
    x <- a[["x"]]
#	org <- a[["org"]]
    
    GOs <- getGOsDownTree(x, org=org)
    sigGOs <- getSigGOTermsOnce(GOs, noGenesInCut=length(x), org=org,
            statisticalTest=statisticalTest, binomAlpha=binomAlpha,
            p.adjust.method=p.adjust.method)

    res <- list(sigGOs=sigGOs[order(sigGOs$pValue), ])
    class(res) <- "GOtree"
    res
}

print.GOtree <- function(x, ... ) {
	print(x$sigGOs)
}

plot.GOtree <- function(x, boxes=25, legendPosition="topright",
        main="Gene Ontology tree, biological processes", ... ) {
    
    GOtermsRelations <- newSearchGOtreeForTerms(x$sigGOs, limit=boxes)	
    if(!is.null(GOtermsRelations)) {
        makeGOtreeRgraphviz(
            a=GOtermsRelations[order(as.numeric(GOtermsRelations$pValue)), ],
            legendPosition=legendPosition, main=main, ... )
    }
}

GOtreeHits <- function(input, inputType="hgu133plus2", org="Hs", GOid,
        returnGeneSymbols = TRUE ) {
    
    a <- convertInput(input, inputType)
    x <- a[["x"]]
    org <- a[["org"]]
    
    GOTermsOntologyBP <- getGOidsWithChildren(getGOidsInOntology(org=org),
            org=org)
    
    GOtermProbeIds <- GOTermsOntologyBP[names(GOTermsOntologyBP) == GOid][[1]]
    res <- GOtermProbeIds[GOtermProbeIds %in% x]
    
    if(returnGeneSymbols) {
        SYMBOL <- get(paste(sep="", "org.", org, ".egSYMBOL"))
        res <- unlist(AnnotationDbi::mget(res, SYMBOL), use.names=FALSE)
        res <- unique(res[!is.na(res)])
    }
    return(res)
}

GOtreeWithLeaveOut <- function(exprsData, inputType="hgu133plus2", org="Hs",
        pc=1 , decreasing=TRUE, noProbes=1000, leaveOut=1,
        runs=NCOL(exprsData)) {
	
    multiRes <-	multiPca(exprsData=exprsData, pc=pc, noProbes=noProbes,
            decreasing=decreasing, leaveOut=leaveOut, runs=runs,
            testFunc=GOtree, inputType=inputType, org=org)

    listName <- c("sigGOs")
    leaveOutRes <- leaveOutFindPvalues(res=multiRes, listName=listName,
            uniqueId="GOid")
    sigGOsRes <- data.frame(
            GOid=multiRes[[1]][[listName]][leaveOutRes[["index"]],"GOid"],
            pValue=leaveOutRes[["medianPvalues"]], genesInTerm=-1,
            totalGenesInTerm=multiRes[[1]][[listName]][leaveOutRes[["index"]],
                    "totalGenesInTerm"],
            GOterm=multiRes[[1]][[listName]][leaveOutRes[["index"]], "GOterm"])
    
    res <- list(sigGOs=sigGOsRes[order(sigGOsRes$pValue), ])
    class(res) <- "GOtree"
    res
}
