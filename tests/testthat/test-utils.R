# Test Helper Functions ---------------------------------------------------

testthat::test_that("Ensure named works", {
    testListUnamed <- list(1, 2, 3, 4)
    testListNamed <- ensureNamed(testListUnamed, prefix = "g_")
    expect_length(testListNamed, length(testListUnamed))
    expect_named(testListNamed, c("g_1", "g_2", "g_3", "g_4"))
    testListNamed <- list(a = 1, b = 2, c = 3, d = 4)
    testListEnsures <- ensureNamed(testListNamed, prefix = "g_")
    expect_length(testListEnsures, length(testListNamed))
    expect_named(testListEnsures, names(testListNamed))
    expect_type(testListEnsures, "list")
})

test_that("Simple Rank Works",{
    testMat <- matrix(runif(20), ncol=4, nrow=5)
    ranked <- simpleRank(testMat)
    expect_equal(nrow(testMat), nrow(ranked))
    expect_equal(ncol(testMat), ncol(ranked))
    orderedMat <- matrix(seq_len(5), ncol=4, nrow=5)
    ranked <- simpleRank(orderedMat)
    expect_equal(ranked, orderedMat)
})

# Test Matrix Argument Handling -------------------------------------------
test_that("Convert Phenotypes Works",{
    testMatrix <- matrix(runif(20), ncol=5)
    colnames(testMatrix) <- c("A", "B", "C", "D", "E")
    testNumIndex <- c(1,3,5)
    testLogIndex <- c(TRUE, FALSE, TRUE, FALSE, TRUE)
    testCharIndex <- c("A", "C", "E")
    convNumIndex <- .convertPhenotype(testNumIndex, testMatrix)
    convLogIndex <- .convertPhenotype(testLogIndex, testMatrix)
    convCharIndex <- .convertPhenotype(testCharIndex, testMatrix)
    expect_equal(convNumIndex, c(1,3,5))
    expect_equal(convLogIndex, c(1,3,5))
    expect_equal(convCharIndex, c(1,3,5))
    expect_error(.convertPhenotype(c(1,3,6), testMatrix),
                 "Phenotype index out of range")
    expect_error(.convertPhenotype(c("A", "B", "F", "G"), testMatrix),
                 paste("Couldn't find all phenotype entries in",
                       " expression matrix. Missing: F, G", sep=""))
    expect_error(.convertPhenotype(c(TRUE, FALSE), testMatrix),
                 paste("If phenotype is logical, must be same length as",
                       " columns of expression matrix", sep=""))
    expect_error(.convertPhenotype(c(), testMatrix), "Phenotype is empty")
    expect_warning(.convertPhenotype(5, testMatrix),
                   paste("Phenotype is only length 1, which can cause ",
                         "issues with later analyses",sep=""))
})

test_that("Convert Network is working",{
    testMatrix <- matrix(runif(20), ncol=5)
    rownames(testMatrix) <- c("A", "B", "C", "D")
    testNumIndex <- c(1,3,4)
    testLogIndex <- c(TRUE, FALSE, TRUE, TRUE)
    testCharIndex <- c("A", "C", "D")
    convNumIndex <- .convertNetwork(testNumIndex, testMatrix)
    convLogIndex <- .convertNetwork(testLogIndex, testMatrix)
    convCharIndex <- .convertNetwork(testCharIndex, testMatrix)
    expect_equal(convNumIndex, c(1,3,4))
    expect_equal(convLogIndex, c(1,3,4))
    expect_equal(convCharIndex, c(1,3,4))
    expect_error(.convertNetwork(c(1,3,6), testMatrix),
                 "Gene network index out of range")
    expect_error(.convertNetwork(c("A", "B", "F", "G"), testMatrix),
                 paste("Couldn't find all gene network entries in",
                       " expression matrix. Missing: F, G", sep=""))
    expect_error(.convertNetwork(c(TRUE, FALSE), testMatrix),
                 paste("If gene network is logical, must be same length as",
                       " rows of expression matrix", sep=""))
    expect_error(.convertNetwork(c(), testMatrix),
                 "Gene network is empty")
    expect_warning(.convertNetwork(4, testMatrix),
                   paste("Gene network is only length 1, which can cause ",
                         "issues with later analyses",sep=""))
})

test_that("Convert Network List works",{
    testMatrix <- matrix(runif(20), ncol=5)
    rownames(testMatrix) <- c("A", "B", "C", "D")
    testNumIndex <- c(1,3,4)
    testLogIndex <- c(TRUE, FALSE, TRUE, TRUE)
    testCharIndex <- c("A", "C", "D")
    netList <- list(testNumIndex, testLogIndex,testCharIndex)
    convNetList <- .convertNetworkList(netList, testMatrix)
    convNetListExp <- list(gn_1 = c(1,3,4), gn_2=c(1,3,4),
                           gn_3 = c(1,3,4))
    expect_equal(convNetList, convNetListExp)
})

# SummarizedExperiment Argument Handling ----------------------------------

test_that("Check Phenotype Works",{
    testSElist <- generateSummarizedExperiment(
        ngenesOrdered = 10,
        ngenesTotal = 20,
        nsamplesOrdered = 15,
        nsamplesTotal = 30,
        dist = rnorm,
        reorderGenes = FALSE,
        genePrefix = "g_",
        samplePrefix = "s_"
    )
    seObj <- testSElist$seObject
    orderedGenes <- testSElist$orderedGenes
    orderedSamples <- testSElist$orderedSamples
    orderedGeneNames <- testSElist$orderedGeneNames
    orderedSampleNames <- testSElist$orderedSampleNames
    unorderedGenes <- testSElist$unorderedGenes
    unorderedSamples <- testSElist$unorderedSamples
    unorderedGeneNames <- testSElist$unorderedGeneNames
    unorderedSampleNames <- testSElist$unorderedSampleNames
    expect_equal(.checkPhenotype(seObject = seObj,
                                 orderedSamples), orderedSamples)
    expect_equal(.checkPhenotype(seObject = seObj,
                                 orderedSampleNames), orderedSamples)
    expect_equal(.checkPhenotype(seObject = seObj,
                                 seObj$phenotypeNum==1), orderedSamples)
    expect_equal(.checkPhenotype(seObject = seObj,
                                 seObj$phenotypeStr=="one"), orderedSamples)
    expect_error(.checkPhenotype(seObj, c(1,2,3,35)),
                 "Phenotype index out of range")
    expect_error(.checkPhenotype(seObj, c("s_1","s_2","s_10","s_34","s_35")),
                 paste("Couldn't find all phenotype entries in",
                       " seObject colnames. Missing: s_34, s_35", sep=""))
    expect_error(.checkPhenotype(seObj, c(TRUE, FALSE, TRUE)),
                 paste("Logical Phenotype Vectors should be ",
                       "the same length as the number of samples, but ",
                       "phenotype has length: ", 3,
                       " not ", 30, sep=""))
    expect_error(.checkPhenotype(seObj, c()), "Phenotype is empty")
    expect_warning(.checkPhenotype(seObj, 5),
                   paste("Phenotype is only length 1, which can cause ",
                         "issues with later analyses",sep=""))

})

test_that("Check Network Works",{
    testSElist <- generateSummarizedExperiment(
        ngenesOrdered = 10,
        ngenesTotal = 20,
        nsamplesOrdered = 15,
        nsamplesTotal = 30,
        dist = rnorm,
        reorderGenes = FALSE,
        genePrefix = "g_",
        samplePrefix = "s_"
    )
    seObj <- testSElist$seObject
    orderedGenes <- testSElist$orderedGenes
    orderedSamples <- testSElist$orderedSamples
    orderedGeneNames <- testSElist$orderedGeneNames
    orderedSampleNames <- testSElist$orderedSampleNames
    unorderedGenes <- testSElist$unorderedGenes
    unorderedSamples <- testSElist$unorderedSamples
    unorderedGeneNames <- testSElist$unorderedGeneNames
    unorderedSampleNames <- testSElist$unorderedSampleNames
    expect_equal(.checkNetwork(network=orderedGenes,
                               numRows=length(seObj),
                               rowNames = rownames(seObj)), orderedGenes)
    expect_equal(.checkNetwork(network=orderedGeneNames,
                               numRows=length(seObj),
                               rowNames = rownames(seObj)), orderedGenes)
    expect_equal(
        .checkNetwork(network=
                          SummarizedExperiment::rowData(seObj)$networkNum==1,
                      numRows=length(seObj),
                      rowNames = rownames(seObj)), orderedGenes)
    expect_equal(
        .checkNetwork(
            network=
                SummarizedExperiment::rowData(seObj)$networkStr=="one",
            numRows=length(seObj),
            rowNames = rownames(seObj)), orderedGenes)
    expect_error(.checkNetwork(c(1,2,3,35),numRows=length(seObj),
                               rowNames = rownames(seObj)),
                 "Gene network index out of range")
    expect_error(.checkNetwork(c("g_1","g_2","g_10","g_34","g_35"),
                               numRows=length(seObj),
                               rowNames = rownames(seObj)),
                 paste("Couldn't find all gene network entries in",
                       " seObject rownames. Missing: g_34, g_35", sep=""))
    expect_error(.checkNetwork(c(TRUE, FALSE, TRUE),
                               numRows=length(seObj),
                               rowNames = rownames(seObj)),
                 paste("Logical Gene Network Vectors should be ",
                       "the same length as the number of genes, but ",
                       "network has has length ", 3,
                       " not ", 20, sep=""))
    expect_error(.checkNetwork(c(),numRows=length(seObj),
                               rowNames = rownames(seObj)), "Network is empty")
    expect_warning(.checkNetwork(5,numRows=length(seObj),
                                 rowNames = rownames(seObj)),
                   paste("Network is only length 1, which can cause ",
                         "issues with later analyses",sep=""))

})

test_that("Check Args Function Works", {
    testSElist <- generateSummarizedExperiment(
        ngenesOrdered = 10,
        ngenesTotal = 20,
        nsamplesOrdered = 15,
        nsamplesTotal = 30,
        dist = rnorm,
        reorderGenes = FALSE,
        genePrefix = "g_",
        samplePrefix = "s_"
    )
    seObj <- testSElist$seObject
    orderedGenes <- testSElist$orderedGenes
    orderedSamples <- testSElist$orderedSamples
    orderedGeneNames <- testSElist$orderedGeneNames
    orderedSampleNames <- testSElist$orderedSampleNames
    unorderedGenes <- testSElist$unorderedGenes
    unorderedSamples <- testSElist$unorderedSamples
    unorderedGeneNames <- testSElist$unorderedGeneNames
    unorderedSampleNames <- testSElist$unorderedSampleNames
    testArray <- matrix(runif(20), ncol=5)
    expect_error(
        .checkArgs(testArray, c(1,2,3), c(4,5,6), list(A=c(1,2,3),
                                                                B=c(4,5,6),
                                                                C=c(7,8,9))),
        "First argument must be a SummarizedExperiment object or subclass")
    expect_equal(
        .checkArgs(seObject = seObj,
                   phenotype1=orderedSamples,
                   phenotype2=unorderedSamples,
                   list(gn1=orderedGenes,
                        gn2=unorderedGenes)),
        list(
            p1Index=orderedSamples,
            p2Index=unorderedSamples,
            geneNetworks=list(gn1=orderedGenes, gn2=unorderedGenes)
        )
    )
    expect_equal(
        .checkArgs(seObject = seObj,
                   phenotype1=orderedSampleNames,
                   phenotype2=unorderedSampleNames,
                   list(gn1=orderedGenes,
                        gn2=unorderedGenes)),
        list(
            p1Index=orderedSamples,
            p2Index=unorderedSamples,
            geneNetworks=list(gn1=orderedGenes, gn2=unorderedGenes)
        )
    )
    expect_equal(
        .checkArgs(seObject = seObj,
                   phenotype1=seObj$phenotypeNum==1,
                   phenotype2=seObj$phenotypeNum==2,
                   list(
                       gn1=SummarizedExperiment::rowData(seObj)$networkNum==1,
                       gn2=SummarizedExperiment::rowData(seObj)$networkNum==2)),
        list(
            p1Index=orderedSamples,
            p2Index=unorderedSamples,
            geneNetworks=list(gn1=orderedGenes, gn2=unorderedGenes)
        )
    )
    expect_equal(
        .checkArgs(
            seObject = seObj,
            phenotype1=seObj$phenotypeStr=="one",
            phenotype2=seObj$phenotypeStr=="two",
            list(
                gn1=SummarizedExperiment::rowData(seObj)$networkStr=="one",
                gn2=SummarizedExperiment::rowData(seObj)$networkStr=="two")),
        list(
            p1Index=orderedSamples,
            p2Index=unorderedSamples,
            geneNetworks=list(gn1=orderedGenes, gn2=unorderedGenes)
        )
    )
})
