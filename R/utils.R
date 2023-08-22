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
