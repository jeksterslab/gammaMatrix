#' Half-Vectorize
#'
#' Half-vectorize a symmetric matrix.
#'
#' Generates an \eqn{m(m + 1) / 2} vector
#' from the unique elements
#' of an \eqn{m \times m} symmetric matrix
#' by stacking the columns (column-major).
#'
#' @details
#' # Dependencies
#' * [vechnames()]
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param x Matrix.
#' @param names Logical.
#'   Add names.
#' @inheritParams vechnames
#'
#' @references
#' [Wikipedia: Half-vectorization](https://en.wikipedia.org/wiki/Vectorization_(mathematics)#Half-vectorization)
#'
#' @returns A vector.
#'
#' @examples
#' x <- matrix(
#'   data = c(
#'     1.0, 0.5, 0.4,
#'     0.5, 1.0, 0.6,
#'     0.4, 0.6, 1.0
#'   ),
#'   ncol = 3
#' )
#'
#' vech(x)
#' @export
#' @family Linear Algebra Functions
#' @keywords linearAlgebra symmetric vectorization
vech <- function(x,
                 names = FALSE,
                 sep = ".") {
  stopifnot(
    is.matrix(x),
    dim(x)[1] == dim(2)
  )
  output <- x[lower.tri(x, diag = TRUE)]
  if (names) {
    varnames <- colnames(x)
    if (is.null(varnames)) {
      varnames <- paste0(seq_len(dim(x)[2]))
    }
    varnames <- vechnames(varnames, sep = sep)
    names(output) <- varnames
  }
  return(output)
}
