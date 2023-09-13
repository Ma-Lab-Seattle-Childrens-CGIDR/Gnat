


bySampleEntropy <- function(method=c("DIRAC", "CRANE", "RACE")){
    # These functions take in filtered expression
    methodFunctionList <- list(
        DIRAC = diracSampleScore,
        CRANE = craneSampleScore,
        RACE = raceSampleScore
    )
    method <- match.arg(method)



}
