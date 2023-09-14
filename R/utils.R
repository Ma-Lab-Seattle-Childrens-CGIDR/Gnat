#' Ensures that the provided list has names
#'
#' @param l list, Input list to which will be returned with names, if it is not
#'    already named,
#' @param prefix string, Prefix for added names
#'
#' @return A named list, with the same data as the input but which is named
#' @export
#'
#' @examples
#' testList <- list(1,2,3,4)
#' testList <- ensureNamed
ensureNamed <- function(in_list, prefix = "gn_") {
    in_list <- as.list(in_list)
    if (is.null(names(in_list))) {
        names(in_list) <- sapply(
            1:length(in_list),
            function(x) paste(prefix, x, sep = "")
        )
    }
    in_list
}


#' Return a matrix where each value has been replaces by that values rank
#'
#' @param filteredExpression Numeric matrix to convert to a rank matrix
#' @param margin Integer, margin to rank along, 1 for rows, 2 for columns
#'
#' @return Numeric matrix of ranks
#' @export
#'
#' @examples
#' testMat <- matrix(runif(20), ncol=4, nrow=5)
#' print(simpleRank(testMat))
#' print(simpleRank(testMat, margin=1))
simpleRank <- function(filteredExpression, margin=2){
    apply(filteredExpression, MARGIN=margin, rank, ties.method="first")
}



# Matrix Argument Handling ------------------------------------------------

.convertPhenotype <- function(phenotype, expressionMatrix){
    if(length(phenotype)==1){
        warning("Phenotype is only length 1, which can cause issues with later",
                " analyses")
    }
    if(length(phenotype)<1){
        stop("Phenotype is empty")
    }
    numCols <- dim(expressionMatrix)[2]
    if(is.logical(phenotype)){
        if(length(phenotype)==numCols){
            pIndex <- seq(numCols)[phenotype]
            return(pIndex)
        } else {
            stop("If phenotype is logical, must be same length as",
                 " columns of expression matrix")
        }
    }
    if(is.numeric(phenotype)){
        if(!(max(phenotype)<=numCols)){
            stop("Phenotype index out of range")
        }
        return(phenotype)
    }
    columnNames <- colnames(expressionMatrix)
    if(is.character(phenotype)){
        if(all(phenotype %in% columnNames)){
            pIndex <- match(phenotype, columnNames)
            return(pIndex)
        } else {
            stop(paste("Couldn't find all phenotype entries in expression ",
                       "matrix. Missing: ",
                       paste(setdiff(phenotype, columnNames), collapse=", "),
                       sep=""))
        }
    }
    stop(paste("Couldn't Parse Phenotype, should be logical, character, ",
               "or numeric vector", sep=""))
}

.convertNetwork <- function(geneNetwork, expressionMatrix){
    if(length(geneNetwork)==1){
        warning("Gene network is only length 1, which can cause issues with",
                " later analyses")
    }
    if(length(geneNetwork)<1){
        stop("Gene network is empty")
    }
    numRows <- dim(expressionMatrix)[1]
    if(is.logical(geneNetwork)){
        if(length(geneNetwork)==numRows){
            gIndex <- seq(numRows)[geneNetwork]
            return(gIndex)
        } else {
            stop(paste("If gene network is logical, must be same length as",
                       " rows of expression matrix", sep=""))
        }
    }
    if(is.numeric(geneNetwork)){
        if(!(max(geneNetwork)<=numRows)){
            stop("Gene network index out of range")
        }
        return(geneNetwork)
    }
    rowNames <- rownames(expressionMatrix)
    if(is.character(geneNetwork)){
        if(all(geneNetwork %in% rowNames)){
            gIndex <- match(geneNetwork, rowNames)
            return(gIndex)
        } else {
            stop(paste("Couldn't find all gene network entries in expression ",
                       "matrix. Missing: ",
                       paste(setdiff(geneNetwork, rowNames), collapse=", "),
                       sep=""))
        }
    }
    stop(paste("Couldn't parse gene network, should be logical, character, ",
               "or numeric vector"))
}

.convertNetworkList <- function(geneNetworkList, expressionMatrix){
    geneNetworkList <- ensureNamed(geneNetworkList)
    lapply(geneNetworkList, .convertNetwork, expressionMatrix=expressionMatrix)
}


# SummarizedExperiment Argument Handling ----------------------------------

.checkArgs <- function(seObject,phenotype1, phenotype2, geneNetworks){
    if(!methods::is(seObject, "SummarizedExperiment")){
        # This will not be called with subclasses (i.e. RangedExperiment)
        stop("First argument must be a SummarizedExperiment object or subclass")
    }
    p1Index <- .checkPhenotype(seObject, phenotype1)
    p2Index <- .checkPhenotype(seObject, phenotype2)
    seDims <- dim(seObject)
    geneNetworks <- ensureNamed(geneNetworks)
    geneNetworks <- lapply(geneNetworks, .checkNetwork, numRows = seDims[1],
                           rownames(seObject))
    list(p1Index=p1Index, p2Index=p2Index, geneNetworks=geneNetworks)
}


.checkPhenotype <- function(seObject, phenotype){
    if(length(phenotype)==1){
        warning("Phenotype is only length 1, which can cause issues with later",
                " analyses")
    }
    if(length(phenotype)<1){
        stop("Phenotype is empty")
    }
    dimSeObject <- dim(seObject)
    if(is.logical(phenotype)){
        if(length(phenotype)!=dimSeObject[2]){
            stop(paste("Logical Phenotype Vectors should be ",
                       "the same length as the number of samples, but ",
                       "phenotype has length: ", length(phenotype),
                       " not ", dimSeObject[2], sep=""))
        }
        pIndex <- seq(dimSeObject[2])[phenotype]
        return(pIndex)
    } else if(is.character(phenotype)){
        columnNames <- SummarizedExperiment::colnames(seObject)
        if(all(phenotype %in% columnNames)){
            pIndex <- match(phenotype, columnNames)
            return(pIndex)
        } else {
            stop(paste("Couldn't find all phenotype entries in seObject ",
                       "colnames. Missing: ",
                       paste(setdiff(phenotype, columnNames), collapse=", "),
                       sep=""))
        }
    } else if(is.numeric(phenotype)){
        if(!(max(phenotype)<=dimSeObject[2])){
            stop("Phenotype index out of range")
        }
        return(phenotype)
    }
    stop(paste("Couldn't Parse Phenotype, should be logical, character, ",
               "or numeric vector", sep=""))
}

.checkNetwork <- function(network, numRows, rowNames){
    if(length(network)==1){
        warning("Network is only length 1, which can cause issues with later",
                " analyses")
    }
    if(length(network)<1){
        stop("Network is empty")
    }
    if(is.logical(network)){
        if(length(network)!=numRows){
            stop(paste("Logical Gene Network Vectors should be ",
                       "the same length as the number of genes, but ",
                       "network has has length ", length(network),
                       " not ", numRows, sep=""))
        }
        gIndex <- seq(numRows)[network]
        return(gIndex)
    } else if(is.character(network)){
        if(all(network %in% rowNames)){
            gIndex <- match(network, rowNames)
            return(gIndex)
        } else {
            stop(paste("Couldn't find all gene network entries in seObject ",
                       "rownames. Missing: ",
                       paste(setdiff(network, rowNames), collapse=", "),
                       sep=""))
        }
    } else if(is.numeric(network)){
        if(!(max(network)<=numRows)){
            stop("Gene network index out of range")
        }
        return(network)
    }
    stop(paste("Couldn't Parse gene network, should be logical, character, ",
               "or numeric vector", sep=""))
}



