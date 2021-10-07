#' Dimension Names for the Gamma Matrix
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param x Matrix or data frame.
#' @param sep Character string.
#'   Separator for variable names.
#' @param mean_structure Logical.
#'   Add names for means.
#'   If `mean_structure = TRUE`,
#'   the first `ncol(x)` elements
#'   will be names for the `ncol(x)` means.
#' @returns A vector.
#'
#' @examples
#' x <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
#' colnames(x) <- rownames(x) <- c("var1", "var2")
#'
#' gammacapnames(x)
#' gammacapnames(x, mean_structure = TRUE)
#' @export
#' @family Gamma Matrix Functions
#' @keywords gammaMatrix
gammacapnames <- function(x,
                          sep = ".",
                          mean_structure = FALSE) {
  varnames <- colnames(x)
  if (is.null(varnames)) {
    varnames <- paste0("v", seq_len(dim(x)[2]))
  }
  output <- vechnames(varnames, sep = sep)
  if (mean_structure) {
    mean_str <- paste0(
      "m",
      sep,
      varnames
    )
    return(
      c(
        mean_str,
        output
      )
    )
  } else {
    return(output)
  }
}
