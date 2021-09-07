#' Dimension Names for the Gamma Matrix
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param x Matrix or data frame.
#' @param sep Character string.
#'   Separator for variable names.
#'
#' @returns A vector.
#'
#' @examples
#' x <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
#' colnames(x) <- rownames(x) <- c("var1", "var2")
#'
#' gammacapnames(x)
#' @export
#' @family Gamma Matrix Functions
#' @keywords gammaMatrix
gammacapnames <- function(x,
                          sep = ".") {
  varnames <- colnames(x)
  if (is.null(varnames)) {
    varnames <- paste0("v", seq_len(dim(x)[2]))
  }
  return(
    vechnames(varnames, sep = sep)
  )
}
