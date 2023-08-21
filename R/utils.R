
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
ensure_named <- function(in_list, prefix="gn_"){
  if(is.null(names(in_list))){
    names(in_list) <- sapply(1:length(in_list),
                             function(x) paste(prefix, x, sep=""))
  }
  in_list
}
