#' Symmetric matrix A from vech(A)
#'
#' Symmetric matrix from its half-vectorization.
#'
#' Generates an \eqn{m \times m} symmetric matrix
#' from an \eqn{m(m + 1) / 2} vector.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param x Vector.
#'
#' @references
#' [Wikipedia: Half-vectorization](https://en.wikipedia.org/wiki/Vectorization_(mathematics)#Half-vectorization)
#'
#' @returns A vector.
#'
#' @examples
#' A <- matrix(
#'   data = c(
#'     1.0, 0.5, 0.4,
#'     0.5, 1.0, 0.6,
#'     0.4, 0.6, 1.0
#'   ),
#'   ncol = 3
#' )
#' vechA <- c(1.0, 0.5, 0.4, 1.0, 0.6, 1.0)
#'
#' sym_of_vech(vechA)
#' @export
#' @family Linear Algebra Functions
#' @keywords linearAlgebra symmetric
sym_of_vech <- function(x) {
  # m = m by m dimensions of the symmetric matrix
  # length(x) = m(m + 1) / 2 solve for m
  # output = symmetric matrix output
  stopifnot(is.vector(x))
  m <- 0.5 * (sqrt(1 + 8 * length(x)) - 1)
  output <- matrix(
    data = 0,
    nrow = m,
    ncol = m
  )
  if (nrow(output) != m) {
    stop("Length of \"x\" is not valid.")
  }
  output[lower.tri(output, diag = TRUE)] <- x
  output[upper.tri(output)] <- t(output)[upper.tri(output)]
  return(output)
}
