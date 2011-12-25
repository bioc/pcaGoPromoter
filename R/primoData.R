primoData <- function(promoters, matrices, cleanUpPromoters=TRUE) {
    
    if(cleanUpPromoters) {
        promoters <- primoData.removeBadSequences(promoters)
        promoters <- primoData.findDublets(promoters)
        promoters <- primoData.findDubletSequences(promoters)
    }
    promoters <- primoData.convertToNumbers(promoters)
    
    thresholds <- primoData.getAllThresholds(promoters, matrices)
    promotersRefseqs <- sapply(promoters, "[", i="refseqs")
    
    res <- list(thresholds=thresholds, promotersRefseqs=promotersRefseqs)
    class(res) <- "primoData"
    return( res )
}

primoData.getPromoters <- function(file) {
    
    require(Biostrings)
    
    promoters <- readFASTA(file, strip.descs=T)
    
    noPromoters <- length(promoters)
    res <- vector("list",noPromoters)
    seqs <- character(noPromoters)
    names <- character(noPromoters)
    
    for(i in 1:noPromoters) {
        res[[i]] <- list()
        
        desc <- strsplit(promoters[[i]]$desc, " ")
        name <- desc[[1]][1]
        strandString <- strsplit(desc[[1]][5], "=")
        if(strandString[[1]][1] != "strand" |
                (strandString[[1]][2] != "+" & strandString[[1]][2] != "-")) {
            print(
        "Strand not 5th element in description??? Or strand not given by +- ?")
            return(NULL)
        }
        strand <- strandString[[1]][2]
        names[i] <- strsplit(name, "_refGene_")[[1]][2]
        if(strand == "-") {
            ## Reverse sequence
            revNumSeq <- 5-rev(primoData.convertPromoterFromTextToNumbers(
                            tolower(promoters[[i]]$seq)))
            revNumSeq[revNumSeq == 0] <- 5
            seqs[i] <- paste(c("a", "c", "g", "t", "n")[revNumSeq], collapse="")
        } else { 
            seqs[i] <- tolower(promoters[[i]]$seq)
        }
    }
    return(list(names, seqs))
}


primoData.convertPromoterFromTextToNumbers <- function(textPromoter) {
    
    numPromoter <- as.numeric(factor(unlist(strsplit(textPromoter , "")),
                    levels=c("a", "c", "g", "t", "n")))
    
    mode(numPromoter) <- "integer"
    return(numPromoter)
}

primoData.convertToNumbers <- function(x) {
    
    for(i in 1:length(x)) {
        x[[i]]$seq <- primoData.convertPromoterFromTextToNumbers(x[[i]]$seq)
    }
    return(x)
}

primoData.removeBadSequences <- function(x) {

    names <- x[[1]]
    seqs <- x[[2]]
    
    res <- logical(length(seqs))
    
    n <- 0
    for(i in 1:length(seqs)) {
        ## Set to TRUE if sequence only contains a,c,g,t,n
        res[i] <- all(!is.na(factor(unlist(strsplit(seqs[[i]], "")),
                                levels=c("a", "c", "g", "t", "n"))))
        if(!res[i]) {
            n <- n+1
            print(paste(n, "Removed entry at index", i,
                            "because it contains illegal characters."))
        } else {
            ## Sequences of just n's are dropped 
            res[i] <- any(is.na(factor(unlist(strsplit(seqs[[i]], "")),
                                    levels=c("n"))))
            if(!res[i]) {
                n <- n+1
                print(paste(n, "Removed entry at index", i,
                                "because it only contains n's."))
            }
        }
    }
    if(n != 0) {
        print(paste(n,
                "entries removed, because they contained illegal characters."))
    }
    
    return(list(names[res], seqs[res]))
}

primoData.findDublets <- function(x) {

    seqs <- x[[2]]
    names <- x[[1]]
    dup <- duplicated(names)
    pos <- which(dup)
    
    n <- 0
    
    ## Only for control. We cannot have 2+ refseqs with same name
    ## and different sequence 
    for(i in 1:length(pos)) {
        hit <- match(names[pos[i]], names[1:pos[i]-1])
        if(seqs[pos[i]] != seqs[hit]) {
            n <- n+1
            print(paste(n, "Entries with name", names[hit],
                            ", but different sequence at", hit, "and", pos[i]))
        }
    }
    
    print(paste(length(pos),"sequences with the same name. Dublets removed. Of these sequences, ",n," had the same name but different sequence. Only the first occurence used."))
    return(list(names[!dup], seqs[!dup]))
}

primoData.findDubletSequences <- function(x) {
    
    names <- as.list(x[[1]])
    seqs <- x[[2]]
    
    dup <- duplicated(seqs)
    pos <- which(dup)
    
    for(i in 1:length(pos)) {
        hit <- match(seqs[pos[i]], seqs[1:pos[i]-1])
        names[[hit]] <- c(names[[hit]], names[[pos[i]]][1])
    }
    
    newNames <- names[!dup]
    newSeqs <- seqs[!dup]
    res <- list()
    for(i in 1:length(newNames)) {
        res[[i]] <- list(refseqs=newNames[[i]], seq=newSeqs[i])
    }
    return(res)
}

primoData.getPwmValuesInPromotersFast <- function(promoters, pwm, pwmLength,
        noPromoters) {
    
    res <- numeric(noPromoters)
    
    pwm <- matrix(c(as.vector(t(pwm)), rep(0, pwmLength)), nrow=5, byrow=TRUE)
    
    for(i in 1:noPromoters) {
        promoter <- promoters[[i]]$seq
        noPwmCalcs <- length(promoter)-pwmLength;
        revPromoter <- rev(5-promoter)
        revPromoter[revPromoter == 0] <- 5 
        m5res <- rep(1, noPwmCalcs+1)
        m3res <- rep(1, noPwmCalcs+1)
        for(n in 1:pwmLength) {
            pwmVec <- pwm[, n]
            noPwmCalcsEnd <- noPwmCalcs+n
            m5res <- m5res*pwmVec[promoter[n:noPwmCalcsEnd]]
            m3res <- m3res*pwmVec[revPromoter[n:noPwmCalcsEnd]]
        }
        res[i] <- max(c(m5res, m3res), na.rm=TRUE)
    }
    return(res)
}

primoData.getThresholdNew <- function( promoters, pwmListEntry ) {
		
    noPromoters <- length(promoters)
    
    pwmId <- names(pwmListEntry)
    pwmBaseId <- pwmListEntry[[1]][["baseId"]]
    pwmName <- pwmListEntry[[1]][["name"]]
    pwm <- pwmListEntry[[1]][["pwm"]]
    pwmLength <- length(pwm[1, ])
    
    ## Add factor 0.2 to all entries in matrix
    pwm <- pwm+0.2
    
    res <- primoData.getPwmValuesInPromotersFast(promoters, pwm, pwmLength,
            noPromoters)
    
    ## Sum of each line (ATGC) in matrix can differ,so we need to do it this way
    pwmSumAcrossMatrix <- 1
    for(i in 1:pwmLength) {
        pwmSumAcrossMatrix <- pwmSumAcrossMatrix*0.25*sum(pwm[, i])
    }
    res <- res/(pwmSumAcrossMatrix)
    
    result <- list(baseId=pwmBaseId[1], name=pwmName, pwmLength=pwmLength,
            maxScores=res)
    return(result)
}

primoData.getAllThresholds <- function(promoters, matrices) {

    func <- function(i) { primoData.getThresholdNew(promoters, matrices[i]) }
    
    if(require("parallel")) {
        res <- mclapply(1:length(matrices), func)
    } else {
        res <- lapply(1:length(matrices), func)
    }
    names(res) <- names(matrices)	
    return(res)
}
